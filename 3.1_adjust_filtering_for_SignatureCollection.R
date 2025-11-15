# R Script: Adjust OmicSignatureCollection Filters

# This script loads an existing OmicSignatureCollection, re-evaluates
# individual tissue signatures based on new Seurat object metadata criteria,
# and saves a new, filtered OmicSignatureCollection.

# 1. Setup and Load Libraries (only essential ones)
library(tidyverse)    # For data manipulation
library(Seurat)       # To load individual Seurat objects
library(OmicSignature) # To work with OmicSignature objects and collections
library(Biobase)      # Often a dependency for OmicSignature
library(Matrix)       # Often a dependency for Seurat

# 2. Define Paths
data_input_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/Tabula_sapiens")
omic_signature_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/results/Tabula_sapiens/MAST/consenus_regression")

# Define the file paths for the existing and new OmicSignatureCollections
input_collection_file <- file.path(omic_signature_output_path, "Tabula_Sapiens_Aging_OmicSignatureCollection.rds")
output_collection_file <- file.path(omic_signature_output_path, "Tabula_Sapiens_Aging_OmicSignatureCollection_corrected-filter.rds") 

message(paste0("Raw Seurat objects expected from: ", data_input_path))
message(paste0("Loading existing OmicSignatureCollection from: ", input_collection_file))

# Check if the input collection file exists
if (!file.exists(input_collection_file)) {
  stop("Input OmicSignatureCollection file not found. Please ensure the path is correct and the file exists.")
}


# Define old analysis parameters for filtering and significance
min_cells_per_tissue <- 100         # Minimum cells required for MAST analysis per tissue
min_expressed_gene_threshold <- 0.1 # Gene must be expressed in at least this percentage of cells
min_genes_after_filter <- 100       # Minimum number of genes to proceed with MAST
adj_p_cutoff <- 0.05                # Adjusted p-value cutoff for significant genes in signature
log2fc_abs_cutoff <- 0.25           # Absolute log2FC cutoff for significant genes in signature
max_genes_in_signature <- 500       # Max. number of significant genes saved in the signature part of the OmicSignature object

# 3. Define the new filter parameters
new_min_subjects <- 3
new_min_age_range <- 10

message("\n--- Identifying Tissues Meeting New Criteria from Raw Files ---")
message(paste0("  Criteria: Minimum subjects >= ", new_min_subjects, " AND Minimum age range >= ", new_min_age_range))

# Get list of all raw tissue files
tissue_files <- list.files(data_input_path, pattern = "TabulaSapiens_.*\\.rds$", full.names = TRUE) 

if (length(tissue_files) == 0) {
  stop("No raw tissue Seurat object files found in ", data_input_path, ". Please check the input path.")
}
message(paste0("Found ", length(tissue_files), " raw tissue files to evaluate."))

# Initialize a vector to store the names of tissues that pass the new criteria
tissues_passing_new_criteria <- c()
tissues_failing_new_criteria_summary <- list() # To store reasons for failure

# 4. Loop through each raw Seurat file to identify tissues that pass the new criteria
for (file_path in tissue_files) {
  # Derive tissue name using the same logic as the original script for consistency
  current_tissue_name_clean <- gsub("_", " ", gsub("TabulaSapiens_|_organ|\\.rds$", "", basename(file_path)))
  
  message(paste0("  Evaluating raw file: ", basename(file_path)))
  
  tryCatch({
    # Load only the Seurat object to access metadata (should be relatively fast)
    tissue_seurat <- readRDS(file_path)
    
    # Recalculate num_subjects: count unique non-NA donor IDs
    unique_donor_ids <- unique(tissue_seurat@meta.data$donor_id)
    num_subjects <- length(unique_donor_ids[!is.na(unique_donor_ids)])
    
    # Recalculate num_age_range: ensure age is numeric, then find range
    current_ages <- tissue_seurat@meta.data$age
    if (!is.numeric(current_ages)) {
      current_ages <- as.numeric(as.character(current_ages))
    }
    
    unique_ages_no_na <- unique(current_ages[!is.na(current_ages)])
    num_age_range <- 0 # Default if insufficient distinct ages
    if (length(unique_ages_no_na) >= 2) {
      num_age_range <- max(unique_ages_no_na) - min(unique_ages_no_na)
    } else {
      # This tissue might not have enough age variation for range calculation, treat as 0
      tissues_failing_new_criteria_summary[[current_tissue_name_clean]] <- paste0("Less than 2 distinct non-NA ages (count: ", length(unique_ages_no_na), ") for age range calculation, treated as 0.")
    }
    
    message(paste0("    '", current_tissue_name_clean, "': Subjects = ", num_subjects, ", Age Range = ", round(num_age_range, 2)))
    
    # Apply the new filter criteria
    if (num_subjects >= new_min_subjects && num_age_range >= new_min_age_range) {
      tissues_passing_new_criteria <- c(tissues_passing_new_criteria, current_tissue_name_clean)
      message(paste0("    -> PASS: '", current_tissue_name_clean, "' meets new raw file criteria."))
    } else {
      fail_reason <- ""
      if (num_subjects < new_min_subjects) {
        fail_reason <- paste0(fail_reason, "Subjects (", num_subjects, ") < ", new_min_subjects, ". ")
      }
      if (num_age_range < new_min_age_range) {
        fail_reason <- paste0(fail_reason, "Age Range (", round(num_age_range, 2), ") < ", new_min_age_range, ". ")
      }
      tissues_failing_new_criteria_summary[[current_tissue_name_clean]] <- trimws(fail_reason)
      message(paste0("    -> FAIL: '", current_tissue_name_clean, "' removed. Reason: ", trimws(fail_reason)))
    }
    
    rm(tissue_seurat); gc(verbose = FALSE) # Free memory after processing each Seurat object
    
  }, error = function(e) {
    warning(paste0("    ERROR: Failed to process raw Seurat object for '", current_tissue_name_clean, "'. Reason: ", e$message, ". Skipping this raw file."))
    tissues_failing_new_criteria_summary[[current_tissue_name_clean]] <- paste0("Processing error: ", e$message)
    rm(tissue_seurat); gc(verbose = FALSE) # Attempt to free memory even on error
  })
}

message(paste0("\nIdentified ", length(tissues_passing_new_criteria), " tissues from raw files that meet the new criteria."))

# 5. Load the existing OmicSignatureCollection
aging_signature_collection <- readRDS(input_collection_file)

# Extract individual OmicSignature objects from the loaded collection
initial_omic_sig_list <- aging_signature_collection$OmicSigList
message(paste0("Existing OmicSignatureCollection contains ", length(initial_omic_sig_list), " OmicSignature objects."))

# Initialize lists to store filtered results
filtered_omic_sig_list <- list()
removed_from_collection_names <- c()

message("\n--- Filtering Existing Collection Based on Raw File Evaluation ---")

# 6. Filter the existing OmicSignatureCollection
for (tissue_name_in_collection in names(initial_omic_sig_list)) {
  current_omic_sig_obj <- initial_omic_sig_list[[tissue_name_in_collection]]
  
  tissue_name_in_collection <- gsub("Aging Signature - ", "", tissue_name_in_collection)
  
  if (tissue_name_in_collection %in% tissues_passing_new_criteria) {
    filtered_omic_sig_list[[tissue_name_in_collection]] <- current_omic_sig_obj
    message(paste0("  KEEP: '", tissue_name_in_collection, "' (passed raw file evaluation)."))
  } else {
    removed_from_collection_names <- c(removed_from_collection_names, tissue_name_in_collection)
    message(paste0("  REMOVE: '", tissue_name_in_collection, "' (did not pass raw file evaluation or was not found in passing list)."))
  }
}

# 7. Update OmicSignatureCollection metadata and create new collection
if (length(filtered_omic_sig_list) > 0) {
  message(paste0("\nNew collection will contain ", length(filtered_omic_sig_list), " OmicSignature objects."))
  
  # Clone the metadata from the original collection to update
  new_collection_metadata <- list(
    collection_name = "Tabula Sapiens Human Aging Signatures - All Tissues", 
    description = paste0("Collection of aging signatures derived from Tabula Sapiens human single-cell RNA-seq data, stratified by tissue, adjusted for sex and donor_id. ",
                         "MAST analysis used, with min cells: ", min_cells_per_tissue, ", min expressed gene threshold: ", min_expressed_gene_threshold * 100, "%, adj. p-value cutoff: ", adj_p_cutoff, 
                         ", |Log2FC| cutoff: ", log2fc_abs_cutoff, ", max genes per signature: ", max_genes_in_signature, ", min. nr. of donors per tissue: ", new_min_subjects, ", min. age range per tissue: ", new_min_age_range, "."),
    organism = "Homo sapiens",
    direction_type = "bi-directional",
    phenotype = "Aging",
    assay_type = "transcriptomics",
    platform = "transcriptomics by single-cell RNA-seq",
    author = "BU_Bioinformatics_ChallengeProject2025",
    year = as.numeric(format(Sys.Date(), "%Y")),
    keywords = c("Aging", "Tabula Sapiens", "single-cell", "MAST", "human", "sapiens")
  )
  
  # Create the new OmicSignatureCollection object
  adjusted_aging_signature_collection <- OmicSignatureCollection$new(
    metadata = new_collection_metadata,
    OmicSigList = filtered_omic_sig_list
  )
  
  # 8. Save the new, adjusted collection
  saveRDS(adjusted_aging_signature_collection, file = output_collection_file)
  message(paste0("\nSuccessfully saved the adjusted OmicSignatureCollection to: ", output_collection_file))
  
} else {
  message("\nNo OmicSignature objects passed the new filtering criteria. An adjusted OmicSignatureCollection was not created.")
}

message("\n--- Filtering Summary ---")
message(paste0("Total raw tissue files found: ", length(tissue_files)))
message(paste0("Raw tissues passing new criteria: ", length(tissues_passing_new_criteria)))
message(paste0("Initial number of signatures in collection: ", length(initial_omic_sig_list)))
message(paste0("Number of signatures kept in new collection: ", length(filtered_omic_sig_list)))
if (length(removed_from_collection_names) > 0) {
  message(paste0("Signatures removed from collection: ", paste(removed_from_collection_names, collapse = ", ")))
}
if (length(tissues_failing_new_criteria_summary) > 0) {
  message("\nRaw tissues that failed filtering and their reasons:")
  print(data.frame(
    Tissue = names(tissues_failing_new_criteria_summary),
    Reason = unlist(tissues_failing_new_criteria_summary),
    row.names = NULL
  ))
}

message("\nScript finished.")
