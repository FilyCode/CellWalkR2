# --- 1. Load Necessary R Packages ---
library(tidyverse)        # For data manipulation (dplyr, stringr, readr::parse_number)
library(Seurat)           # For single-cell object management
library(cellxgene.census) # For fetching data from the Census API
library(SingleCellExperiment) # For converting SCE to Seurat if initial fetch fails

# --- 2. Define Paths and Global Variables ---
# Output directory for individual Seurat objects
data_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/Tabula_sapiens")
dir.create(data_output_path, recursive = TRUE, showWarnings = FALSE)
message(paste0("Individual Seurat objects will be saved to: ", data_output_path))

# Census collection ID for Tabula Sapiens
collection_id_all_tissues <- "e5f58829-1a66-40b5-a624-9046778e74f5"

# Specific Census release version
census_release_version <- "2025-01-30" # IMPORTANT: Update to the latest stable release 

# --- 3. Establish Census Connection and retrieve target collection ---
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

message("\nAccessing Homo sapiens metadata to list all unique tissues...")
# Access Homo sapiens obs metadata
human_obs <- census$get("census_data")$get("homo_sapiens")$obs

# Read only metadata columns to get all unique specific tissue names within Tabula Sapiens collection
obs_df <- human_obs$read()$concat() |>
  as.data.frame() |>
  dplyr::select(dataset_id, tissue_general) # Extract dataset_id and tissue_general

message("\nFiltering target dataset for Tabula Sapiens collection...")
obs_tabula <- obs_df |> filter(dataset_id %in% tabula_sapiens_census_dataset_ids)

message("\nExtracting unique specific tissues in the Tabula Sapiens collection...")
# List all available specific tissues within the Tabula Sapiens collection
unique_tissues <- sort(unique(obs_tabula$tissue_general))
message(paste0("Found ", length(unique_tissues), " unique specific tissues in Tabula Sapiens collection:"))
print(unique_tissues)

message("\nStarting iterative fetching, preprocessing, and saving data by specific tissue.")

# Define required metadata columns for local processing (original names from Census)
required_meta_cols_for_processing <- c("development_stage", "sex", "donor_id")

for(t in unique_tissues) {
  message(paste0("\n--- Processing specific tissue: ", t, " ---"))
  
  # Construct the filter string for fetching specific cells from the Census API
  obs_value_string <- sprintf(
    'tissue_general == "%s" & is_primary_data == TRUE & dataset_id %%in%% c("%s")',
    t,
    paste(tabula_sapiens_census_dataset_ids, collapse = '","')
  )
  
  seu_tissue <- NULL # Initialize seu_tissue for each iteration
  
  message("  Retrieving cells and creating Seurat object from Census...")
  tryCatch({
    seu_tissue <- cellxgene.census::get_seurat(
      census = census,
      organism = "Homo sapiens",
      measurement_name = "RNA", 
      obs_value_filter = obs_value_string
    )
    message(paste0("  Successfully fetched Seurat object (", ncol(seu_tissue), " cells)."))
  }, error = function(e) {
    message("  Error fetching Seurat object using get_seurat. Attempting to fetch as SingleCellExperiment and convert.")
    message("  Error details: ", e$message)
    temp_sce_obj <- NULL # Initialize to NULL
    tryCatch({ # Nested tryCatch for get_single_cell_experiment to catch potential issues there too
      temp_sce_obj <- cellxgene.census::get_single_cell_experiment(
        census = census,
        organism = "Homo sapiens",
        obs_value_filter = obs_value_string
      )
      # Included memory management, make sure the objects not needed anymore are cleared from memory
      on.exit({ if (!is.null(temp_sce_obj)) { rm(temp_sce_obj); gc() } }, add = TRUE)
      
      # Check if SCE object has data before converting
      if (is.null(temp_sce_obj) || ncol(temp_sce_obj) == 0) {
        message("  No cells retrieved with SingleCellExperiment either. Skipping this tissue.")
        return(NULL) # Skip to next iteration, as this block is inside a tryCatch
      }
      seu_tissue <- as.Seurat(temp_sce_obj, counts = "counts", data = "logcounts")
      message(paste0("  Successfully converted SingleCellExperiment to Seurat object (", ncol(seu_tissue), " cells)."))
      rm(temp_sce_obj); gc() # Clear memory
    }, error = function(sce_e) {
      message(paste0("  Error fetching with SingleCellExperiment for tissue '", t, "': ", sce_e$message))
      message("  Skipping this tissue.")
      return(NULL) # Ensures seu_tissue remains NULL or invalid to be caught by next check
    })
  })
  
  # Skip if no Seurat object was successfully created or it's empty
  if (is.null(seu_tissue) || ncol(seu_tissue) == 0) {
    message(paste0("  No cells retrieved or object empty for tissue '", t, "'. Skipping this tissue."))
    next
  }
  
  
  # --- Local Metadata Processing and Filtering ---
  
  message("  Preprocessing metadata and normalizing data...")
  
  # Ensure all required metadata columns are present for processing
  missing_cols_initial <- required_meta_cols_for_processing[!(required_meta_cols_for_processing %in% colnames(seu_tissue@meta.data))]
  if (length(missing_cols_initial) > 0) {
    message(paste0("  Skipping tissue due to missing essential metadata for filtering. Missing: ", paste(missing_cols_initial, collapse=", ")))
    rm(seu_tissue); gc(); next # Memory management
  }
  
  # Process metadata: rename, parse age, convert to factors
  metadata_df <- seu_tissue@meta.data %>%
    dplyr::mutate(
      age_original_string = as.character(development_stage), # Keep development_stage, use age_original_string for parsing
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
      
      age = round(age, 2), # Round fractional ages
      
      sex = as.factor(sex),
      donor_id = as.factor(donor_id)
    ) %>%
    dplyr::select(-age_original_string) # Remove temporary column used for parsing
  
  seu_tissue@meta.data <- metadata_df # Assign processed metadata back to Seurat object
  rm(metadata_df); gc() # Clear temporary metadata data frame for memory management
  
  # Filter out cells with any missing essential metadata (age, sex, donor_id)
  initial_cells_tissue <- ncol(seu_tissue)
  cells_to_keep <- !is.na(seu_tissue$age) & !is.na(seu_tissue$sex) & !is.na(seu_tissue$donor_id)
  
  if (sum(cells_to_keep) == 0) {
    message("  No cells remaining after metadata filtering (all cells have NA for essential metadata). Skipping this tissue.")
    rm(seu_tissue); gc(); next # Memory management
  }
  
  seu_tissue <- seu_tissue[, cells_to_keep]
  message(paste0("  Filtered out ", initial_cells_tissue - ncol(seu_tissue), " cells due to missing metadata. Remaining cells: ", ncol(seu_tissue)))
  
  # --- Explicit Normalization ---
  message("  Running default Seurat normalization (LogNormalize) on full counts and updating 'data' slot.")
  seu_tissue <- NormalizeData(seu_tissue)
  
  # --- Saving Seurat object ---
  safe_tissue_name <- gsub("[^[:alnum:]_]", "_", t) # Make filename safe using the specific tissue name 't'
  output_seurat <- paste0(data_output_path, '/TabulaSapiens_', safe_tissue_name, '.rds')
  
  message(paste0("  Saving Seurat object for '", t, "' (", ncol(seu_tissue), " cells) to: ", output_seurat))
  saveRDS(seu_tissue, file = output_seurat)
  
  # --- Memory Management ---
  rm(seu_tissue); gc() # Clear memory for the current tissue
  message(paste0("--- Finished processing and clearing memory for specific tissue '", t, "' ---"))
} # End of main loop through specific tissues

# --- 6. Final Cleanup and Completion Message ---
census$close()
message("\nCellxgene Census connection closed.")
message("\n--- Script Complete: All Tabula Sapiens specific tissues fetched, preprocessed, and saved. ---")
message("You can now proceed with downstream analysis on these individual tissue files.")