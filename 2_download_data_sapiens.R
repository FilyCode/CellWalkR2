# Phase 1: Fetch Tabula Sapiens Data, Preprocess, Split by Tissue, and Save

# Load necessary R packages
library(tidyverse)        # Data manipulation
library(Seurat)           # Single-cell object management
library(cellxgene.census) # Fetch data from Census API
library(SingleCellExperiment) # For converting to Seurat if needed

# Define Paths and Variables
# Output directory for individual tissue Seurat objects
data_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/Tabula_sapiens")
dir.create(data_output_path, recursive = TRUE, showWarnings = FALSE)
message(paste0("Individual tissue Seurat objects will be saved to: ", data_output_path))

# Census identifiers for the Tabula Sapiens collection
collection_id_all_tissues <- "e5f58829-1a66-40b5-a624-9046778e74f5"
census_release_version <- "2025-01-30" # IMPORTANT: Update this if the stable release changes

# Open Cellxgene Census SOMA Connection
message("\nOpening Cellxgene Census SOMA connection")
census <- cellxgene.census::open_soma(census_version = census_release_version)
message(paste0("Opened Census connection for version: ", census_release_version))

# Identify Tabula Sapiens datasets within the Census
message("\nIdentifying Tabula Sapiens datasets within the Census")
tabula_sapiens_census_dataset_ids <- NULL

tryCatch({
  census_datasets_metadata <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())
  tabula_sapiens_datasets <- census_datasets_metadata %>%
    dplyr::filter(collection_id == collection_id_all_tissues)
  
  if (nrow(tabula_sapiens_datasets) > 0) {
    tabula_sapiens_census_dataset_ids <- unique(tabula_sapiens_datasets$dataset_id)
    message(paste0("Found ", length(tabula_sapiens_census_dataset_ids), " Census dataset_ids corresponding to Tabula Sapiens collection."))
    message("First 5 identified Census dataset_ids: ", paste(head(tabula_sapiens_census_dataset_ids, 5), collapse = ", "))
  } else {
    stop(paste0("No datasets found in Census 'datasets' table matching collection_id: ", collection_id_all_tissues, ". Cannot proceed."))
  }
}, error = function(e) {
  stop(paste0("Error during Census 'datasets' metadata retrieval: ", e$message))
})

message("\nStarting iterative fetching, preprocessing, and saving dataset-by-dataset to minimize peak memory usage.")

# Define required metadata columns from the Census for subsequent analysis
required_meta_cols <- c("development_stage", "sex", "tissue", "donor_id")

# Main Loop: Process each Tabula Sapiens dataset ID individually
for (current_dataset_id in tabula_sapiens_census_dataset_ids) {
  message(paste0("\nProcessing data for Census dataset ID: ", current_dataset_id))
  
  current_dataset_filter_string <- paste0("dataset_id == '", current_dataset_id, "'")
  temp_seurat_obj <- NULL
  
  # Fetch data for this single dataset ID
  tryCatch({
    temp_seurat_obj <- cellxgene.census::get_seurat(
      census = census,
      organism = "Homo sapiens",
      obs_value_filter = current_dataset_filter_string
    )
    message(paste0("Successfully fetched Seurat object for dataset '", current_dataset_id, "' (", ncol(temp_seurat_obj), " cells)."))
  }, error = function(e) {
    message(paste0("Error fetching Seurat object for dataset '", current_dataset_id, "'. Attempting to fetch as SingleCellExperiment and convert."))
    message("Error details: ", e$message)
    temp_sce_obj <- cellxgene.census::get_single_cell_experiment(
      census = census,
      organism = "Homo sapiens",
      obs_value_filter = current_dataset_filter_string
    )
    temp_seurat_obj <- as.Seurat(temp_sce_obj, counts = "counts", data = "logcounts")
    message(paste0("Successfully fetched SingleCellExperiment and converted to Seurat object for dataset '", current_dataset_id, "' (", ncol(temp_seurat_obj), " cells)."))
    rm(temp_sce_obj); gc()
  })
  
  if (is.null(temp_seurat_obj) || ncol(temp_seurat_obj) == 0) {
    message(paste0("No cells retrieved for dataset '", current_dataset_id, "'. Skipping this dataset."))
    next
  }
  
  message(paste0("  Preprocessing metadata and normalizing data for dataset '", current_dataset_id, "'..."))
  
  # Check for essential metadata columns
  missing_cols_initial <- required_meta_cols[!(required_meta_cols %in% colnames(temp_seurat_obj@meta.data))]
  if (length(missing_cols_initial) > 0) {
    essential_for_subset <- c("development_stage", "sex", "tissue", "donor_id")
    if(any(essential_for_subset %in% missing_cols_initial)){
      message(paste0("  Skipping dataset '", current_dataset_id, "' due to missing essential metadata for filtering. Missing: ", paste(intersect(essential_for_subset, missing_cols_initial), collapse=", ")))
      rm(temp_seurat_obj); gc(); next
    }
  }
  
  metadata_df <- temp_seurat_obj@meta.data # Work on a copy of metadata
  
  # Robust age parsing (handles 'X-year-old stage', 'X-month-old stage', and categorical stages)
  metadata_df <- metadata_df %>%
    dplyr::rename(
      age = development_stage, # Rename original 'development_stage' to 'age'
      subject = donor_id       # Rename 'donor_id' to 'subject'
    ) %>%
    dplyr::mutate(
      age_original_string = as.character(age), # Convert to character for parsing
      age = NA_real_, # Initialize numeric age to NA
      
      age = ifelse(str_detect(age_original_string, "^\\d+-year-old stage|\\d+ year(s)? old"),
                   readr::parse_number(age_original_string), age),
      age = ifelse(str_detect(age_original_string, "^\\d+-month-old stage"),
                   readr::parse_number(age_original_string) / 12, age),
      
      # Impute numeric ages for specific categories (adjust values as needed)
      age = ifelse(age_original_string == "young adult stage", 25, age),
      age = ifelse(age_original_string == "adult stage", 40, age),
      age = ifelse(age_original_string == "elderly stage", 70, age),
      age = ifelse(age_original_string == "newborn human stage", 0, age),
      
      age = round(age, 2), # Round fractional ages
      
      # Ensure key metadata columns are factors
      sex = as.factor(sex),
      tissue = as.factor(tissue),
      subject = as.factor(subject)
    ) %>%
    dplyr::select(-age_original_string) # Remove temporary parsing column
  
  temp_seurat_obj@meta.data <- metadata_df # Assign processed metadata back
  rm(metadata_df); gc() 
  
  # Filter out cells with missing key metadata
  initial_cells_dataset <- ncol(temp_seurat_obj)
  cells_to_keep <- !is.na(temp_seurat_obj$age) &
    !is.na(temp_seurat_obj$sex) &
    !is.na(temp_seurat_obj$tissue) &
    !is.na(temp_seurat_obj$subject)
  
  if (sum(cells_to_keep) == 0) {
    message(paste0("  No cells remaining for dataset '", current_dataset_id, "' after metadata filtering (all cells have NA for age/sex/tissue/subject). Skipping this dataset."))
    rm(temp_seurat_obj); gc(); next
  }
  
  temp_seurat_obj <- temp_seurat_obj[, cells_to_keep]
  message(paste0("  Filtered out ", initial_cells_dataset - ncol(temp_seurat_obj), " cells due to missing metadata. Remaining cells: ", ncol(temp_seurat_obj)))
  
  # Normalize data if 'data' slot is empty
  if (is.null(GetAssayData(temp_seurat_obj, slot = "data")) || all(GetAssayData(temp_seurat_obj, slot = "data") == 0)) {
    message("  Seurat 'data' slot is empty or zeros. Running default Seurat normalization (LogNormalize).")
    temp_seurat_obj <- NormalizeData(temp_seurat_obj)
  }
  
  # Split this single dataset's data by tissue and save to disk
  message(paste0("  Splitting data from dataset '", current_dataset_id, "' by tissue and saving..."))
  unique_tissues_in_this_dataset <- unique(temp_seurat_obj@meta.data$tissue)
  message(paste0("  Found ", length(unique_tissues_in_this_dataset), " unique tissues in dataset '", current_dataset_id, "'. Cells will be saved if >0."))
  
  for (current_tissue_in_dataset in unique_tissues_in_this_dataset) {
    safe_tissue_name <- gsub("[^[:alnum:]_]", "_", current_tissue_in_dataset)
    output_filename <- file.path(data_output_path, paste0("tabula_sapiens_", safe_tissue_name, ".rds"))
    
    if (file.exists(output_filename)) {
      message(paste0("    Skipping '", current_tissue_in_dataset, "': file already exists at ", output_filename))
      next
    }
    
    tissue_seurat_to_save <- subset(temp_seurat_obj, subset = tissue == current_tissue_in_dataset)
    
    if (ncol(tissue_seurat_to_save) > 0) {
      saveRDS(tissue_seurat_to_save, file = output_filename)
      message(paste0("      Saved ", ncol(tissue_seurat_to_save), " cells for '", current_tissue_in_dataset, "' to ", output_filename))
    } else {
      message(paste0("      No cells found for '", current_tissue_in_dataset, "' in dataset '", current_dataset_id, "' after subsetting. Skipping save."))
    }
    rm(tissue_seurat_to_save); gc()
  }
  
  rm(temp_seurat_obj); gc() # Clear memory for this dataset
  message(paste0("Finished processing dataset '", current_dataset_id, "'"))
  
} # End of main loop through dataset IDs

# Final cleanup
census$close()
message("\nFinal Cellxgene Census connection closed.")
message("\nFinished Phase 1: All datasets processed, tissues split and saved.")
message("You can now run 'analyze_tissues_for_signatures.R' to process these individual tissue files.")

