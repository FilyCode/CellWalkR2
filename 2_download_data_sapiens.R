
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



# Construct the filter string to fetch all Tabula Sapiens data
ts_filter_string <- paste0("dataset_id %in% c('", paste(tabula_sapiens_census_dataset_ids, collapse = "', '"), "')")
message(paste0("\nUsing filter string to fetch all Tabula Sapiens data: ", ts_filter_string))

message("\n--- Fetching ALL Tabula Sapiens data into a single Seurat object ---")
message("This will require significant memory (currently ~250GB based on your report) and time. Monitor your cluster job's memory usage!")

seurat_obj <- NULL # Initialize outside tryCatch for scope

tryCatch({
  seurat_obj <- cellxgene.census::get_seurat(
    census = census,
    organism = "Homo sapiens",
    obs_value_filter = ts_filter_string
  )
  message("Successfully fetched Seurat object from Cellxgene Census.")
}, error = function(e) {
  message("Error fetching Seurat object. Attempting to fetch as SingleCellExperiment and convert.")
  message("Error details: ", e$message)
  sce_obj <- cellxgene.census::get_single_cell_experiment(
    census = census,
    organism = "Homo sapiens",
    obs_value_filter = ts_filter_string
  )
  seurat_obj <- as.Seurat(sce_obj, counts = "counts", data = "logcounts")
  message("Successfully fetched SingleCellExperiment and converted to Seurat object.")
})

# Close census connection to free up resources as early as possible after data is retrieved
census$close()
message("Cellxgene Census connection closed.")

if (is.null(seurat_obj) || ncol(seurat_obj) == 0) {
  stop("Failed to retrieve any cells from Tabula Sapiens. Exiting.")
}

message(paste0("Initial Seurat object contains ", ncol(seurat_obj), " cells and ", nrow(seurat_obj), " features."))





# --- Preprocessing and Data Preparation (on the large object) ---
message("\n--- Preprocessing metadata and normalizing data (on the full dataset) ---")
# Define required metadata columns from the Census for subsequent analysis
required_meta_cols <- c("development_stage", "sex", "tissue", "donor_id")
missing_cols <- required_meta_cols[!(required_meta_cols %in% colnames(seurat_obj@meta.data))]
if (length(missing_cols) > 0) {
  stop("Missing one or more required metadata columns from the fetched dataset: ", paste(missing_cols, collapse = ", "))
}

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  dplyr::rename(
    age = development_stage, # Rename 'development_stage' to 'age'
    subject = donor_id       # Rename 'donor_id' to 'subject' for the model formula
  ) %>%
  dplyr::mutate(
    age = str_replace_all(as.character(age), " years| year|yrs|yr", ""), # Clean age string
    age = suppressWarnings(as.numeric(age)), # Convert to numeric
    sex = as.factor(sex),
    tissue = as.factor(tissue),
    subject = as.factor(subject)
  )

# Filter out cells with missing or unparseable key metadata
initial_cells <- ncol(seurat_obj)
seurat_obj <- subset(seurat_obj, subset = !is.na(age) & !is.na(sex) & !is.na(tissue) & !is.na(subject))
message(paste0("Filtered out ", initial_cells - ncol(seurat_obj), " cells due to missing or invalid metadata. Remaining cells: ", ncol(seurat_obj)))

# Ensure the 'data' assay slot contains normalized (log-transformed) counts for MAST
# If the 'data' slot is empty or not suitable, run Seurat's NormalizeData()
if (is.null(GetAssayData(seurat_obj, slot = "data")) || all(GetAssayData(seurat_obj, slot = "data") == 0)) {
  message("Seurat 'data' slot is empty or all zeros. Running default Seurat normalization (LogNormalize).")
  seurat_obj <- NormalizeData(seurat_obj)
}




# --- Split by tissue and save to disk ---
message("\n--- Splitting the full dataset by tissue and saving to individual .rds files ---")
unique_tissues_all_data <- unique(seurat_obj@meta.data$tissue)
message(paste0("Found ", length(unique_tissues_all_data), " unique tissues across the dataset."))

for (current_tissue in unique_tissues_all_data) {
  # Sanitize tissue name for a valid filename (replace non-alphanumeric with underscore)
  safe_tissue_name <- gsub("[^[:alnum:]_]", "_", current_tissue)
  output_filename <- file.path(data_output_path, paste0("tabula_sapiens_", safe_tissue_name, ".rds"))
  
  # Check if file already exists to avoid re-processing if script is re-run
  if (file.exists(output_filename)) {
    message(paste0("  Skipping '", current_tissue, "': file already exists at ", output_filename))
    next # Move to the next tissue
  }
  
  message(paste0("  Subsetting and saving tissue: ", current_tissue))
  tissue_seurat_to_save <- subset(seurat_obj, subset = tissue == current_tissue)
  
  if (ncol(tissue_seurat_to_save) > 0) {
    saveRDS(tissue_seurat_to_save, file = output_filename)
    message(paste0("    Saved ", ncol(tissue_seurat_to_save), " cells for '", current_tissue, "' to ", output_filename))
  } else {
    message(paste0("    No cells found for '", current_tissue, "' after subsetting. Skipping save."))
  }
}



# Clean up the large Seurat object to free memory after all tissues are saved
rm(seurat_obj)
gc() # Force garbage collection
message("\n--- Finished Phase 1: All tissues split and saved. Large Seurat object removed from memory. ---")
message("You can now run 'analyze_tissues_for_signatures.R' to process these individual tissue files.")


