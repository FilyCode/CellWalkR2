# --- 1. Load Necessary R Packages ---
library(tidyverse)        # For data manipulation (dplyr, stringr, readr::parse_number)
library(Seurat)           # For single-cell object management
library(cellxgene.census) # For fetching data from the Census API

# --- 2. Define Paths and Global Variables ---
# Output directory for individual Seurat objects
data_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/Tabula_sapiens_test")
dir.create(data_output_path, recursive = TRUE, showWarnings = FALSE)
message(paste0("Individual Seurat objects will be saved to: ", data_output_path))

# Census collection ID for Tabula Sapiens
collection_id_all_tissues <- "e5f58829-1a66-40b5-a624-9046778e74f5"

# Specific Census release version
census_release_version <- "2025-01-30" # IMPORTANT: Update to the latest stable release 

# --- 3. Establish Census Connection and retrieving target collection ---
message("\nEstablishing connection to Cellxgene Census (version ", 
        census_release_version, ")")
census <- cellxgene.census::open_soma(census_version = census_release_version)

message("Identifying Tabula Sapiens datasets within the Census...")
tabula_sapiens_census_dataset_ids <- NULL

tryCatch({
  census_datasets_metadata <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat())
  tabula_sapiens_datasets <- census_datasets_metadata %>%
    dplyr::filter(collection_id == collection_id_all_tissues)
  
  if (nrow(tabula_sapiens_datasets) > 0) {
    tabula_sapiens_census_dataset_ids <- unique(tabula_sapiens_datasets$dataset_id)
    message(paste0("Found ", length(tabula_sapiens_census_dataset_ids), 
                   " Census dataset_ids for Tabula Sapiens collection. First 5: ", 
                   paste(head(tabula_sapiens_census_dataset_ids, 5), collapse = ", ")))
  } else {
    stop(paste0("No datasets found in Census 'datasets' table matching collection_id: ", 
                collection_id_all_tissues, ". Cannot proceed."))
  }
}, error = function(e) {
  stop(paste0("Error during Census 'datasets' metadata retrieval: ", e$message))
})

message("\nAccessing Homo sapiens metadata...")
# Access Homo sapiens obs metadata
human_obs <- census$get("census_data")$get("homo_sapiens")$obs

# Read only metadata columns 
obs_df <- human_obs$read()$concat() |>
  as.data.frame() |>
  dplyr::select(dataset_id, tissue) # Extract dataset_id and tissue

message("\nFiltering target dataset...")
obs_tabula <- obs_df |> filter(dataset_id %in% tabula_sapiens_census_dataset_ids)

message("\nExtracting unique tissues in the dataset...")
# List available general tissues
unique_tissues <- sort(unique(obs_tabula$tissue))
print(unique_tissues)

message("\nStarting iterative fetching, preprocessing, and saving data by tissue.")

for(t in unique_tissues) {
  message(paste0("Current tissue: ", t))
  
  message("Retrieving cells...")
  tryCatch({
    obs_value_string <- sprintf(
      'tissue == "%s" & is_primary_data == TRUE & dataset_id %%in%% c("%s")',
      t,
      paste(tabula_sapiens_census_dataset_ids, collapse = '","')
    )
  }, error = function(e) {
    stop(paste0("Error during cells retrieval: ", e$message))
  })
  
  message("Creating Seurat object...")
  
  seu_tissue <- cellxgene.census::get_seurat(
    census = census,
    organism = "Homo sapiens",
    measurement_name = "RNA",
    obs_value_filter = obs_value_string
  )
  
  message("Saving Seurat object...")
  output_seurat <- paste0(data_output_path, '/TabulaSapiens_', t, '.rds')
  
  message(paste0("Filename: ", output_seurat))
  saveRDS(seu_tissue, output_seurat)
}

