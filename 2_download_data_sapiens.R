# Phase 1: Fetch Tabula Sapiens Data, Preprocess, Split by Tissue, and Save

# --- 1. Load Necessary R Packages ---
library(tidyverse)        # For data manipulation (dplyr, stringr, readr::parse_number)
library(Seurat)           # For single-cell object management
library(cellxgene.census) # For fetching data from the Census API
library(SingleCellExperiment) # Required for converting SCE to Seurat if initial fetch fails

# --- 2. Define Paths and Global Variables ---
# Output directory for individual Seurat objects, aggregated by general tissue type
data_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/Tabula_sapiens")
dir.create(data_output_path, recursive = TRUE, showWarnings = FALSE)
message(paste0("Individual general tissue Seurat objects will be saved to: ", data_output_path))

# Census collection ID for Tabula Sapiens and the specific Census release version
collection_id_all_tissues <- "e5f58829-1a66-40b5-a624-9046778e74f5"
census_release_version <- "2025-01-30" # IMPORTANT: Update this to the latest stable release for consistency

# --- 3. Establish Census Connection and Identify Datasets ---
message("\nEstablishing connection to Cellxgene Census (version ", census_release_version, ")")
census <- cellxgene.census::open_soma(census_version = census_release_version)

message("Identifying Tabula Sapiens datasets within the Census...")
tabula_sapiens_census_dataset_ids <- NULL

tryCatch({
  census_datasets_metadata <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())
  tabula_sapiens_datasets <- census_datasets_metadata %>%
    dplyr::filter(collection_id == collection_id_all_tissues)
  
  if (nrow(tabula_sapiens_datasets) > 0) {
    tabula_sapiens_census_dataset_ids <- unique(tabula_sapiens_datasets$dataset_id)
    message(paste0("Found ", length(tabula_sapiens_census_dataset_ids), " Census dataset_ids for Tabula Sapiens collection. First 5: ", paste(head(tabula_sapiens_census_dataset_ids, 5), collapse = ", ")))
  } else {
    stop(paste0("No datasets found in Census 'datasets' table matching collection_id: ", collection_id_all_tissues, ". Cannot proceed."))
  }
}, error = function(e) {
  stop(paste0("Error during Census 'datasets' metadata retrieval: ", e$message))
})

message("\nStarting iterative fetching, preprocessing, and saving data dataset-by-dataset.")
message("This approach processes one dataset at a time to minimize peak memory usage on the cluster.")

# Define required metadata columns (original names from Census)
required_meta_cols <- c("development_stage", "sex", "tissue", "donor_id", "tissue_general")

# --- 4. Main Loop: Process Each Tabula Sapiens Dataset Individually ---
for (current_dataset_id in tabula_sapiens_census_dataset_ids) {
  message(paste0("\n--- Processing Census dataset ID: ", current_dataset_id, " ---"))
  
  current_dataset_filter_string <- paste0("dataset_id == '", current_dataset_id, "'")
  temp_seurat_obj <- NULL
  
  # Fetch data for the current dataset ID
  tryCatch({
    temp_seurat_obj <- cellxgene.census::get_seurat(
      census = census,
      organism = "Homo sapiens",
      obs_value_filter = current_dataset_filter_string
    )
    message(paste0("Successfully fetched Seurat object (", ncol(temp_seurat_obj), " cells)."))
  }, error = function(e) {
    message("  Error fetching Seurat object. Attempting to fetch as SingleCellExperiment and convert.")
    message("  Error details: ", e$message)
    temp_sce_obj <- cellxgene.census::get_single_cell_experiment(
      census = census,
      organism = "Homo sapiens",
      obs_value_filter = current_dataset_filter_string
    )
    temp_seurat_obj <- as.Seurat(temp_sce_obj, counts = "counts", data = "logcounts")
    message(paste0("  Successfully converted SingleCellExperiment to Seurat object (", ncol(temp_seurat_obj), " cells)."))
    rm(temp_sce_obj); gc()
  })
  
  if (is.null(temp_seurat_obj) || ncol(temp_seurat_obj) == 0) {
    message(paste0("  No cells retrieved for dataset. Skipping this dataset."))
    next
  }
  
  message("  Preprocessing metadata and normalizing data...")
  
  # Ensure all required metadata columns are present for processing
  missing_cols_initial <- required_meta_cols[!(required_meta_cols %in% colnames(temp_seurat_obj@meta.data))]
  if (length(missing_cols_initial) > 0) {
    message(paste0("  Skipping dataset due to missing essential metadata for filtering and splitting. Missing: ", paste(missing_cols_initial, collapse=", ")))
    rm(temp_seurat_obj); gc(); next
  }
  
  # Process metadata: rename, parse age, convert to factors
  metadata_df <- temp_seurat_obj@meta.data %>%
    dplyr::rename(
      age = development_stage, # Rename 'development_stage' to 'age'
      subject = donor_id       # Rename 'donor_id' to 'subject'
    ) %>%
    dplyr::mutate(
      age_original_string = as.character(age), # Convert to character for robust parsing
      age = NA_real_, # Initialize numeric age to NA
      
      # Parse numeric age from "X-year-old stage", "X years old", "X-month-old stage"
      age = ifelse(str_detect(age_original_string, "^\\d+-year-old stage|\\d+ year(s)? old"),
                   readr::parse_number(age_original_string), age),
      age = ifelse(str_detect(age_original_string, "^\\d+-month-old stage"),
                   readr::parse_number(age_original_string) / 12, age),
      
      # Impute numeric ages for common categorical stages (adjust values as needed)
      age = ifelse(age_original_string == "young adult stage", 25, age),
      age = ifelse(age_original_string == "adult stage", 40, age),
      age = ifelse(age_original_string == "elderly stage", 70, age),
      age = ifelse(age_original_string == "newborn human stage", 0, age),
      
      age = round(age, 2), # Round fractional ages (e.g., from months to years)
      
      # Use 'tissue_general' as the primary 'tissue' label for broader categories
      tissue = as.factor(tissue_general), 
      sex = as.factor(sex),
      subject = as.factor(subject)
    ) %>%
    dplyr::select(-age_original_string, -tissue_general) # Remove temporary and now redundant columns
  
  temp_seurat_obj@meta.data <- metadata_df # Assign processed metadata back to Seurat object
  rm(metadata_df); gc() # Clear temporary metadata data frame
  
  # Filter out cells with any missing essential metadata (age, sex, tissue, subject)
  initial_cells_dataset <- ncol(temp_seurat_obj)
  cells_to_keep <- !is.na(temp_seurat_obj$age) &
    !is.na(temp_seurat_obj$sex) &
    !is.na(temp_seurat_obj$tissue) & 
    !is.na(temp_seurat_obj$subject)
  
  if (sum(cells_to_keep) == 0) {
    message("  No cells remaining after metadata filtering (all cells have NA for essential metadata). Skipping this dataset.")
    rm(temp_seurat_obj); gc(); next
  }
  
  temp_seurat_obj <- temp_seurat_obj[, cells_to_keep]
  message(paste0("  Filtered out ", initial_cells_dataset - ncol(temp_seurat_obj), " cells due to missing metadata. Remaining cells: ", ncol(temp_seurat_obj)))
  
  # Normalize raw counts if the 'data' slot is empty (contains log-normalized counts)
  if (is.null(GetAssayData(temp_seurat_obj, slot = "data")) || all(GetAssayData(temp_seurat_obj, slot = "data") == 0)) {
    message("  Seurat 'data' slot is empty or zeros. Running default Seurat normalization (LogNormalize).")
    temp_seurat_obj <- NormalizeData(temp_seurat_obj)
  }
  
  # --- 5. Split Data by General Tissue and Save to Disk ---
  message("  Splitting data by GENERAL tissue type and saving individual Seurat objects...")
  unique_general_tissues_in_this_dataset <- unique(temp_seurat_obj@meta.data$tissue)
  message(paste0("  Found ", length(unique_general_tissues_in_this_dataset), " unique GENERAL tissues in this dataset. Objects will be saved if cells > 0."))
  
  for (current_general_tissue in unique_general_tissues_in_this_dataset) {
    safe_general_tissue_name <- gsub("[^[:alnum:]_]", "_", current_general_tissue)
    output_filename <- file.path(data_output_path, paste0("tabula_sapiens_", safe_general_tissue_name, ".rds"))
    
    # Check if a file for this general tissue already exists from a previous dataset.
    # If so, this script prioritizes keeping the first one encountered for simplicity.
    # To aggregate all cells for a general tissue across *all* datasets, a more complex merge/save logic would be needed here.
    if (file.exists(output_filename)) {
      message(paste0("    Skipping saving '", current_general_tissue, "': file already exists. (This general tissue might have been processed from an earlier dataset)."))
      next
    }
    
    tissue_seurat_to_save <- subset(temp_seurat_obj, subset = tissue == current_general_tissue)
    
    if (ncol(tissue_seurat_to_save) > 0) {
      saveRDS(tissue_seurat_to_save, file = output_filename)
      message(paste0("      Saved ", ncol(tissue_seurat_to_save), " cells for general tissue '", current_general_tissue, "' to ", output_filename))
    } else {
      message(paste0("      No cells found for general tissue '", current_general_tissue, "' in this dataset after subsetting. Skipping save."))
    }
    rm(tissue_seurat_to_save); gc()
  }
  
  rm(temp_seurat_obj); gc() # Clear memory for the current dataset
  message(paste0("--- Finished processing and clearing memory for dataset '", current_dataset_id, "' ---"))
  
} # End of main loop through dataset IDs

# --- 6. Final Cleanup and Completion Message ---
census$close()
message("\nCellxgene Census connection closed.")
message("\n--- Phase 1 Complete: All datasets processed, tissues split and saved (by general tissue). ---")
message("You can now run 'analyze_tissues_for_signatures.R' to perform downstream analysis on these individual general tissue files.")

