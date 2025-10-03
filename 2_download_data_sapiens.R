
# Phase 1: Fetch All Tabula Sapiens Data, Preprocess, Split by Tissue, and Save

# 1. Setup and Load Libraries
library(tidyverse)      
library(Seurat)         
library(cellxgene.census) 

# --- Define Paths and Variables ---
data_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/sapiens/")
dir.create(data_output_path, recursive = TRUE, showWarnings = FALSE)
message(paste0("Individual tissue Seurat objects will be saved to: ", data_output_path))

# Census-related identifiers for Tabula Sapiens
collection_id_all_tissues <- "e5f58829-1a66-40b5-a624-9046778e74f5"
census_release_version <- "2025-01-30" 


# --- Open Cellxgene Census SOMA Connection ---
message("\n--- Opening Cellxgene Census SOMA connection ---")
census <- cellxgene.census::open_soma(census_version = census_release_version)
message(paste0("Opened Census connection for version: ", census_release_version))


# --- Identify Tabula Sapiens datasets within the Census ---
message("\n--- Identifying Tabula Sapiens datasets within the Census ---")
tabula_sapiens_census_dataset_ids <- NULL

tryCatch({
  # Get the 'datasets' table from census_info
  census_datasets_metadata <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())
  
  # Filter for datasets whose `collection_id` matches the Tabula Sapiens one
  tabula_sapiens_datasets <- census_datasets_metadata %>%
    dplyr::filter(collection_id == collection_id_all_tissues)
  
  if (nrow(tabula_sapiens_datasets) > 0) {
    tabula_sapiens_census_dataset_ids <- unique(tabula_sapiens_datasets$dataset_id)
    message(paste0("Found ", length(tabula_sapiens_census_dataset_ids), " Census dataset_ids corresponding to Tabula Sapiens collection (", collection_id_all_tissues, ")."))
    message("First 5 identified Census dataset_ids:")
    print(head(tabula_sapiens_census_dataset_ids, 5))
  } else {
    stop(paste0("No datasets found in Census 'datasets' table matching collection_id: ", collection_id_all_tissues, ". Cannot proceed."))
  }
}, error = function(e) {
  stop(paste0("Error during Census 'datasets' metadata retrieval: ", e$message))
})

# IMPORTANT: The 'census' object remains open and will be passed to get_seurat/get_single_cell_experiment.

message("\n--- Starting iterative fetching, preprocessing, and saving dataset-by-dataset ---")
message("This approach processes one dataset at a time to minimize peak memory usage.")

# Define required metadata columns from the Census for subsequent analysis
required_meta_cols <- c("development_stage", "sex", "tissue", "donor_id")


# --- Main Loop: Iterate through each identified Tabula Sapiens dataset ID ---
for (current_dataset_id in tabula_sapiens_census_dataset_ids) {
  message(paste0("\n--- Processing data for Census dataset ID: ", current_dataset_id, " ---"))
  
  # Construct filter for this single dataset ID
  current_dataset_filter_string <- paste0("dataset_id == '", current_dataset_id, "'")
  
  temp_seurat_obj <- NULL # Initialize for this iteration
  
  # Fetch data for *this single dataset ID*
  tryCatch({
    # --- CRITICAL CORRECTION HERE: Pass the 'census' object, REMOVE 'census_version' ---
    temp_seurat_obj <- cellxgene.census::get_seurat(
      census = census, # Pass the already opened census object
      organism = "Homo sapiens",
      obs_value_filter = current_dataset_filter_string
    )
    message(paste0("Successfully fetched Seurat object for dataset '", current_dataset_id, "' (", ncol(temp_seurat_obj), " cells)."))
  }, error = function(e) {
    message(paste0("Error fetching Seurat object for dataset '", current_dataset_id, "'. Attempting to fetch as SingleCellExperiment and convert."))
    message("Error details: ", e$message)
    # --- CRITICAL CORRECTION HERE: Pass the 'census' object, REMOVE 'census_version' ---
    temp_sce_obj <- cellxgene.census::get_single_cell_experiment(
      census = census, # Pass the already opened census object
      organism = "Homo sapiens",
      obs_value_filter = current_dataset_filter_string
    )
    # Ensure the assay name is correct if not "counts"
    temp_seurat_obj <- as.Seurat(temp_sce_obj, counts = "counts", data = "logcounts")
    message(paste0("Successfully fetched SingleCellExperiment and converted to Seurat object for dataset '", current_dataset_id, "' (", ncol(temp_seurat_obj), " cells)."))
    rm(temp_sce_obj); gc() # Clean up SCE object immediately
  })
  
  if (is.null(temp_seurat_obj) || ncol(temp_seurat_obj) == 0) {
    message(paste0("No cells retrieved for dataset '", current_dataset_id, "'. Skipping this dataset."))
    next # Move to the next dataset_id
  }
  
  # --- Preprocessing and Data Preparation (on the single dataset object) ---
  message(paste0("  Preprocessing metadata and normalizing data for dataset '", current_dataset_id, "'..."))
  
  # Check for required metadata columns
  missing_cols <- required_meta_cols[!(required_meta_cols %in% colnames(temp_seurat_obj@meta.data))]
  if (length(missing_cols) > 0) {
    message(paste0("  Warning: Missing one or more required metadata columns from dataset '", current_dataset_id, "': ", paste(missing_cols, collapse = ", ")))
    # If critical columns are missing, decide whether to skip or proceed with NAs
    # For now, we'll proceed, but NAs will be filtered later.
  }
  
  temp_seurat_obj@meta.data <- temp_seurat_obj@meta.data %>%
    dplyr::rename(
      age = development_stage, # Rename 'development_stage' to 'age'
      subject = donor_id       # Rename 'donor_id' to 'subject' for the model formula
    ) %>%
    dplyr::mutate(
      age = str_replace_all(as.character(age), " years| year|yrs|yr", ""), # Clean age string
      age = suppressWarnings(as.numeric(age)), # Convert to numeric, non-numeric become NA
      sex = as.factor(sex),
      tissue = as.factor(tissue),
      subject = as.factor(subject)
    )
  
  # Filter out cells with missing or unparseable key metadata
  initial_cells_dataset <- ncol(temp_seurat_obj)
  temp_seurat_obj <- subset(temp_seurat_obj, subset = !is.na(age) & !is.na(sex) & !is.na(tissue) & !is.na(subject))
  message(paste0("  Filtered out ", initial_cells_dataset - ncol(temp_seurat_obj), " cells from dataset '", current_dataset_id, "' due to missing or invalid metadata. Remaining cells: ", ncol(temp_seurat_obj)))
  
  if (ncol(temp_seurat_obj) == 0) {
    message(paste0("  No cells remaining for dataset '", current_dataset_id, "' after metadata filtering. Skipping this dataset."))
    rm(temp_seurat_obj); gc(); next # Clean up and move to next dataset
  }
  
  # Ensure the 'data' assay slot contains normalized (log-transformed) counts for MAST
  if (is.null(GetAssayData(temp_seurat_obj, slot = "data")) || all(GetAssayData(temp_seurat_obj, slot = "data") == 0)) {
    message("  Seurat 'data' slot is empty or all zeros. Running default Seurat normalization (LogNormalize).")
    temp_seurat_obj <- NormalizeData(temp_seurat_obj)
  }
  
  # --- Split this single dataset's data by tissue and save to disk ---
  message(paste0("  Splitting data from dataset '", current_dataset_id, "' by tissue and saving..."))
  unique_tissues_in_this_dataset <- unique(temp_seurat_obj@meta.data$tissue)
  message(paste0("  Found ", length(unique_tissues_in_this_dataset), " unique tissues in dataset '", current_dataset_id, "'."))
  
  for (current_tissue_in_dataset in unique_tissues_in_this_dataset) {
    # Sanitize tissue name for a valid filename
    safe_tissue_name <- gsub("[^[:alnum:]_]", "_", current_tissue_in_dataset)
    output_filename <- file.path(data_output_path, paste0("tabula_sapiens_", safe_tissue_name, ".rds"))
    
    # Check if file already exists to avoid re-processing if script is re-run
    # This is crucial for long runs, if a dataset fails, you don't re-save its already good tissues
    if (file.exists(output_filename)) {
      message(paste0("    Skipping '", current_tissue_in_dataset, "': file already exists at ", output_filename))
      next # Move to the next tissue
    }
    
    message(paste0("    Subsetting and saving tissue: ", current_tissue_in_dataset))
    tissue_seurat_to_save <- subset(temp_seurat_obj, subset = tissue == current_tissue_in_dataset)
    
    if (ncol(tissue_seurat_to_save) > 0) {
      saveRDS(tissue_seurat_to_save, file = output_filename)
      message(paste0("      Saved ", ncol(tissue_seurat_to_save), " cells for '", current_tissue_in_dataset, "' to ", output_filename))
    } else {
      message(paste0("      No cells found for '", current_tissue_in_dataset, "' in dataset '", current_dataset_id, "' after subsetting. Skipping save."))
    }
    rm(tissue_seurat_to_save); gc() # Clear memory for this tissue object
  }
  
  # Clean up the temporary Seurat object for this dataset to free memory
  rm(temp_seurat_obj)
  gc() # Force garbage collection
  message(paste0("--- Finished processing and clearing memory for dataset '", current_dataset_id, "' ---"))
  
} # End of main loop through dataset IDs

# --- Close the ONE global census connection after the loop completes ---
census$close()
message("Final Cellxgene Census connection closed.")

message("\n--- Finished Phase 1: All datasets processed, tissues split and saved. ---")
message("You can now run 'analyze_tissues_for_signatures.R' to process these individual tissue files.")
