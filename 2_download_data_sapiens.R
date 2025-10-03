# 1. Setup and Load Libraries
library(Biobase)
library(SummarizedExperiment)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(Matrix)
library(DelayedArray)
library(cellxgene.census)
library(MAST) # For differential expression analysis
library(OmicSignature) # For storing the results



# Set a path for saving the OmicSignatures
omic_signature_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/results/sapiens/")
dir.create(omic_signature_output_path, recursive = TRUE, showWarnings = FALSE)
message(paste0("OmicSignature results will be saved to: ", omic_signature_output_path))

# Data Retrieval using cellxgene.census
# This ID refers to the entire Tabula Sapiens Human dataset, allowing us to subset by tissue later
collection_id_all_tissues <- "e5f58829-1a66-40b5-a624-9046778e74f5"
dataset_id_all_tissues <- "946fa48d-a0ac-4e5b-80fc-1d96cb5083a7" # Tabula Sapiens - All Tissues (H5AD)

message("Opening Cellxgene Census SOMA connection...")
# Specify the census_version as recommended for consistency
census_release_version <- "2025-01-30" # Update this if the stable release changes
census <- cellxgene.census::open_soma(census_version = census_release_version)
message(paste0("Opened Census connection for version: ", census_release_version))

message(paste0("Fetching Tabula Sapiens 'All Tissues' dataset (ID: ", dataset_id_all_tissues, ") into a Seurat object. This may take some time and significant memory on the cluster."))





# 2. Fetch data for the specific dataset_id directly into a Seurat object.
# The `obs_value_filter` argument is used for filtering cells based on metadata.
tryCatch({
  seurat_obj <- cellxgene.census::get_seurat(
    census = census,
    organism = "Homo sapiens",
    obs_value_filter = paste0("dataset_id == '", dataset_id_all_tissues, "'") # --- !!! CORRECTED ARGUMENT !!! ---
  )
  message("Successfully fetched Seurat object from Cellxgene Census.")
}, error = function(e) {
  message("Error fetching Seurat object. Attempting to fetch as SingleCellExperiment and convert.")
  message("Error details: ", e$message)
  sce_obj <- cellxgene.census::get_single_cell_experiment(
    census = census,
    organism = "Homo sapiens",
    obs_value_filter = paste0("dataset_id == '", dataset_id_all_tissues, "'") # --- !!! CORRECTED ARGUMENT !!! ---
  )
  # Convert to Seurat if SCE was fetched successfully
  # Ensure the assay name is correct if not "counts"
  seurat_obj <- as.Seurat(sce_obj, counts = "counts", data = "logcounts") # Assuming normalized data in 'logcounts' slot
  message("Successfully fetched SingleCellExperiment and converted to Seurat object.")
})

# Close census connection to free up resources
census$close()
message("Cellxgene Census connection closed.")




# 3. Preprocessing and Data Preparation

# Check for required metadata columns and rename for consistency with the model formula
required_cols <- c("development_stage", "sex", "tissue", "donor_id")
if (!all(required_cols %in% colnames(seurat_obj@meta.data))) {
  stop("Missing one or more required metadata columns from the fetched dataset: ",
       paste(required_cols[!(required_cols %in% colnames(seurat_obj@meta.data))], collapse = ", "))
}

message("Processing metadata...")
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  dplyr::rename(
    age = development_stage, # Rename 'development_stage' to 'age'
    subject = donor_id       # Rename 'donor_id' to 'subject' for the model formula
  ) %>%
  # Robustly convert age to numeric, handling various string formats and NAs
  dplyr::mutate(
    age = str_replace_all(as.character(age), " years| year|yrs|yr", ""), # Remove common age unit strings
    age = suppressWarnings(as.numeric(age)), # Convert to numeric, non-numeric values become NA
    sex = as.factor(sex),
    tissue = as.factor(tissue),
    subject = as.factor(subject)
  )

# Filter out cells with missing or unparseable age, sex, tissue, or subject
initial_cells <- ncol(seurat_obj)
seurat_obj <- subset(seurat_obj, subset = !is.na(age) & !is.na(sex) & !is.na(tissue) & !is.na(subject))
message(paste0("Filtered out ", initial_cells - ncol(seurat_obj), " cells due to missing or invalid metadata. Remaining cells: ", ncol(seurat_obj)))

# Ensure the 'data' assay slot contains normalized (log-transformed) counts for MAST
# If the 'data' slot is empty or not suitable, you might need to run Seurat's NormalizeData()
if (is.null(GetAssayData(seurat_obj, slot = "data")) || all(GetAssayData(seurat_obj, slot = "data") == 0)) {
  message("Seurat 'data' slot is empty or all zeros. Running default Seurat normalization (LogNormalize).")
  seurat_obj <- NormalizeData(seurat_obj)
}






# 4. Loop through tissues and calculate aging signatures
unique_tissues <- unique(seurat_obj@meta.data$tissue)
message(paste0("\nFound ", length(unique_tissues), " unique tissues for analysis."))

# Initialize an OmicSignatureCollection to store all tissue-specific signatures
omicsig_collection_metadata <- OmicSignature::createMetadata(
  signature_name = "Tabula Sapiens Human Aging Signatures - All Tissues",
  organism = "Homo Sapiens",
  direction_type = "bi-directional",
  phenotype = "Aging",
  description = "Collection of aging signatures derived from Tabula Sapiens human single-cell RNA-seq data, stratified by tissue, adjusted for sex and donor_id (subject).",
  author = "ChallengeProject2025",
  year = as.numeric(format(Sys.Date(), "%Y")),
  keywords = c("Aging", "Tabula Sapiens", "single-cell", "MAST", "human", "sapiens")
)
aging_signature_collection <- OmicSignatureCollection$new(
  metadata = omicsig_collection_metadata,
  OmicSigList = list()
)

for (current_tissue in unique_tissues) {
  message(paste0("\n--- Processing tissue: ", current_tissue, " ---"))

  # Subset Seurat object for the current tissue
  tissue_seurat <- subset(seurat_obj, subset = tissue == current_tissue)

  # Check for minimum data requirements for MAST
  min_cells <- 100 # Minimum cells required for MAST per tissue
  if (ncol(tissue_seurat) < min_cells) {
    message(paste0("Skipping '", current_tissue, "' due to insufficient cells (", ncol(tissue_seurat), " < ", min_cells, ")."))
    next
  }

  num_subjects <- length(unique(tissue_seurat@meta.data$subject))
  num_sex_groups <- length(unique(tissue_seurat@meta.data$sex))
  num_distinct_ages <- length(unique(tissue_seurat@meta.data$age))

  # Heuristic checks for model fitting: need at least 2 subjects, 2 sex groups, and 2 distinct ages
  if (num_subjects < 2 || num_sex_groups < 2 || num_distinct_ages < 2) {
    message(paste0("Skipping '", current_tissue, "' due to insufficient variation for regression (Subjects: ", num_subjects, ", Sex groups: ", num_sex_groups, ", Distinct ages: ", num_distinct_ages, ")."))
    next
  }
  
  # Prepare data for MAST: Convert Seurat object to SingleCellExperiment
  # MAST works well with `SingleCellExperiment` objects
  sce_tissue <- as.SingleCellExperiment(tissue_seurat, assay = "data") # Use the normalized 'data' slot

  # Ensure factors are re-leveled after subsetting
  colData(sce_tissue)$subject <- droplevels(colData(sce_tissue)$subject)
  colData(sce_tissue)$sex <- droplevels(colData(sce_tissue)$sex)

  # Create a MAST 'SingleCellAssay' (SCA) object
  # Filter out subjects that have only one cell within this subset, as `zlm` might struggle
  subject_counts <- table(colData(sce_tissue)$subject)
  subjects_to_keep <- names(subject_counts[subject_counts > 1])
  
  if (length(subjects_to_keep) < 2 && num_subjects >= 2) { # If we started with >1 subject but filtering reduced it
    message(paste0("Warning: In '", current_tissue, "', filtering reduced subjects with >1 cell to ", length(subjects_to_keep), ". Continuing with remaining data."))
  } else if (length(subjects_to_keep) < 2 && num_subjects < 2) {
      message(paste0("Skipping '", current_tissue, "' due to insufficient subjects with more than one cell (after filtering)."))
      next
  }
  
  sce_tissue <- sce_tissue[, colData(sce_tissue)$subject %in% subjects_to_keep]
  colData(sce_tissue)$subject <- droplevels(colData(sce_tissue)$subject) # Re-droplevels after filtering

  # Re-check cell count after all metadata/subject filtering
  if (ncol(sce_tissue) < min_cells) {
    message(paste0("Skipping '", current_tissue, "' due to insufficient cells (", ncol(sce_tissue), " < ", min_cells, ") after subject filtering."))
    next
  }

  # Filter genes to include only those expressed in a certain percentage of cells
  # This helps MAST perform better by focusing on more relevant genes
  expressed_threshold <- 0.1 # Gene expressed in at least 10% of cells
  expressed_genes <- rowSums(assay(sce_tissue, "data") > 0) / ncol(sce_tissue) > expressed_threshold
  if (sum(expressed_genes) < 10) { # Arbitrary minimum number of expressed genes
    message(paste0("Skipping '", current_tissue, "' due to insufficient highly expressed genes (", sum(expressed_genes), ")."))
    next
  }
  sce_tissue_filtered <- sce_tissue[expressed_genes, ]

  message(paste0("Running MAST for '", current_tissue, "' with ", nrow(sce_tissue_filtered), " genes and ", ncol(sce_tissue_filtered), " cells."))

  # Fit MAST model: gene ~ age + sex + subject
  # 'subject' (donor_id) is treated as a fixed effect to account for individual variability.
  # 'method = 'glm'' for generalized linear model, 'ebayes = TRUE' for empirical Bayes shrinkage.
  tryCatch({
    # Convert SingleCellExperiment to a MAST `SingleCellAssay` object
    sca_mast <- SceToSingleCellAssay(sce_tissue_filtered, class = "SingleCellAssay")

    # Fit the ZLM model
    zlm_obj <- zlm(~ age + sex + subject, sca = sca_mast, method = 'glm', ebayes = TRUE)

    # Extract coefficients and p-values for the 'age' term
    # `logFC` function from MAST provides the log-fold change, p-value, and FDR for specified contrasts.
    results_table_mast <- MAST::as.data.frame(logFC(zlm_obj, contrasts = "age"))

    if (is.null(results_table_mast) || nrow(results_table_mast) == 0) {
      message(paste0("No differential expression results found for 'age' in tissue: ", current_tissue))
      next
    }

    # Prepare results for OmicSignature
    results_table_omic <- results_table_mast %>%
      dplyr::mutate(
        symbol = PrimerID, # Gene symbol
        score = `logFC`,   # logFC as the signature score
        p_value = `Pvalue`,
        adj_p = `FDR`      # Adjusted p-value (False Discovery Rate)
      ) %>%
      dplyr::select(symbol, score, p_value, adj_p) %>%
      dplyr::mutate(
        direction = ifelse(score > 0, "+", "-") # Determine direction based on logFC
      )

    
    
    
    # 5. Create OmicSignature object for the current tissue

    # Define metadata for the tissue-specific signature
    metadata_tissue_sig <- OmicSignature::createMetadata(
      signature_name = paste0("Aging Signature - ", current_tissue),
      organism = "Homo Sapiens",
      direction_type = "bi-directional",
      phenotype = paste0("Aging in ", current_tissue),
      covariates = "sex, subject (donor_id)",
      platform = "Single-cell RNA-seq (cellxgene.census/Tabula Sapiens)",
      sample_type = paste0(current_tissue, " cells"),
      logfc_cutoff = NULL, # Cutoffs applied during signature filtering, not in metadata
      p_value_cutoff = NULL,
      adj_p_cutoff = 0.05, # Default adjusted p-value cutoff for significance
      score_cutoff = 0.25, # Default absolute logFC cutoff for significance
      keywords = c("Aging", current_tissue, "Tabula Sapiens", "single-cell", "MAST"),
      author = "User Name / AI Assistant",
      PMID = NULL, # Add Tabula Sapiens publication PMID if available/desired
      year = as.numeric(format(Sys.Date(), "%Y")),
      description = paste0("Aging signature derived from Tabula Sapiens human single-cell RNA-seq data for the ", current_tissue, " tissue. Differential expression calculated with MAST, adjusting for sex and donor_id.")
    )

    # Filter significant genes for the signature based on defined cutoffs
    maxQ <- metadata_tissue_sig$adj_p_cutoff
    minScore <- metadata_tissue_sig$score_cutoff

    sig_genes <- results_table_omic %>%
      dplyr::filter(adj_p <= maxQ & abs(score) >= minScore) %>%
      dplyr::select(symbol, score, direction)

    if (nrow(sig_genes) == 0) {
      message(paste0("No significant genes found for 'age' in tissue: ", current_tissue, " with current cutoffs (adj_p <= ", maxQ, ", |logFC| >= ", minScore, ")."))
      next # Skip to next tissue if no significant genes
    }

    # Create the OmicSignature object
    omic_sig_tissue <- OmicSignature$new(
      metadata = metadata_tissue_sig,
      signature = sig_genes,
      difexp = results_table_omic # Store the full differential expression results
    )

    # Add the OmicSignature object to the collection
    aging_signature_collection$OmicSigList[[current_tissue]] <- omic_sig_tissue
    message(paste0("Successfully created and added aging signature for ", current_tissue, ". (", nrow(sig_genes), " significant genes)"))

    # Save individual OmicSignature object (recommended for easier access)
    safe_tissue_name <- gsub("[^[:alnum:]_]", "_", current_tissue) # Sanitize for filename
    saveRDS(omic_sig_tissue, file = file.path(omic_signature_output_path, paste0("aging_signature_", safe_tissue_name, "_oSig.rds")))

  }, error = function(e) {
    message(paste0("Error during MAST or OmicSignature creation for tissue '", current_tissue, "': ", e$message))
  })
}





# 6. Save the complete OmicSignatureCollection
if (length(aging_signature_collection$OmicSigList) > 0) {
  saveRDS(aging_signature_collection, file = file.path(omic_signature_output_path, "Tabula_Sapiens_Aging_OmicSignatureCollection.rds"))
  message(paste0("\nSaved OmicSignatureCollection with ", length(aging_signature_collection$OmicSigList), " tissue signatures to '", omic_signature_output_path, "'."))
} else {
  message("\nNo aging signatures were successfully generated for any tissue and added to the collection.")
}

message("\nScript finished.")