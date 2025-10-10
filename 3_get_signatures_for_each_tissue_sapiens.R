
# Load Individual Tissue Files and Perform Aging Signature Analysis

# 1. Setup and Load Libraries
library(tidyverse)        
library(Seurat)           
library(SingleCellExperiment) 
library(MAST)             
library(OmicSignature)    
library(Biobase)
library(doParallel)


# --- Define Paths and Variables ---
data_input_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/Tabula_sapiens/")
omic_signature_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/results/Tabula_sapiens/")

dir.create(omic_signature_output_path, recursive = TRUE, showWarnings = FALSE)
message(paste0("Input tissue Seurat objects expected from: ", data_input_path))
message(paste0("OmicSignature results will be saved to: ", omic_signature_output_path))


# Define analysis parameters
min_cells_per_tissue <- 100         # Minimum cells required for MAST per tissue
min_expressed_gene_threshold <- 0.1 # Gene expressed in at least x% of cells
min_genes_after_filter <- 10        # Minimum number of genes to proceed with MAST
adj_p_cutoff <- 0.05                # Adjusted p-value cutoff for significant genes in signature
score_cutoff <- 0.25                # Absolute logFC cutoff for significant genes in signature

# Determine number of CPU cores for parallel processing with MAST's zlm
num_cores <- as.numeric(Sys.getenv("NSLOTS", unset = 1)) 

if (num_cores > 1) {
  registerDoParallel(cores = num_cores)
  message(paste0("\n  Registered parallel backend for MAST zlm with ", num_cores, " cores."))
} else {
  message("\n  Running MAST zlm in serial mode (1 core).")
}

# --- Initialize an OmicSignatureCollection ---
message("\n--- Initializing OmicSignatureCollection ---")
omicsig_collection_metadata <- list(
  collection_name = "Tabula Sapiens Human Aging Signatures - All Tissues", # This is the missing required field!
  description = paste0("Collection of aging signatures derived from Tabula Sapiens human single-cell RNA-seq data, stratified by tissue, adjusted for sex and donor_id. ",
                       "MAST analysis used, with min cells: ", min_cells_per_tissue, ", min expressed gene threshold: ", min_expressed_gene_threshold * 100, "%, adj. p-value cutoff: ", adj_p_cutoff, ", |logFC| cutoff: ", score_cutoff, "."),
  organism = "Homo sapiens",
  direction_type = "bi-directional",
  phenotype = "Aging",
  assay_type = "transcriptomics",
  platform = "transcriptomics by single-cell RNA-seq",
  author = "BU_Bioinformatics_ChallengeProject2025",
  year = as.numeric(format(Sys.Date(), "%Y")),
  keywords = c("Aging", "Tabula Sapiens", "single-cell", "MAST", "human", "sapiens")
)

# Initialize an empty R list to collect individual OmicSignature objects
# The OmicSignatureCollection object itself will be created after the loop.
all_tissue_omicsigs <- list()


# --- Get list of saved tissue files ---
tissue_files <- list.files(data_input_path, pattern = "TabulaSapiens_.*\\.rds$", full.names = TRUE)
if (length(tissue_files) == 0) {
  stop("No tissue Seurat object files found in ", data_input_path, ".")
}
message(paste0("Found ", length(tissue_files), " tissue files to analyze."))

tissue_files <- list("/restricted/projectnb/agedisease/projects/challenge2025/data/Tabula_sapiens/TabulaSapiens_bladder_organ.rds")

# --- Loop through individual tissue files and perform analysis ---
for (file_path in tissue_files) {
  current_tissue_name <- sub("^tabula_sapiens_|_\\.rds$", "", basename(file_path))
  message(paste0("\n--- Analyzing tissue from file: ", basename(file_path), " (", current_tissue_name, ") ---"))
  
  # Load the tissue-specific Seurat object
  tissue_seurat <- readRDS(file_path)
  
  # Ensure the object is not empty after loading (shouldn't be if saved correctly)
  if (ncol(tissue_seurat) == 0) {
    message(paste0("  Skipping '", current_tissue_name, "': loaded object is empty."))
    rm(tissue_seurat); gc(); next
  }
  
  
  
  # --- Data Checks for MAST ---
  if (ncol(tissue_seurat) < min_cells_per_tissue) {
    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(tissue_seurat), " < ", min_cells_per_tissue, ")."))
    rm(tissue_seurat); gc(); next
  }
  
  num_subjects <- length(levels(tissue_seurat@meta.data$donor_id))
  num_sex_groups <- length(levels(tissue_seurat@meta.data$sex))
  num_distinct_ages <- length(unique(tissue_seurat@meta.data$age))
  
  if (num_subjects < 2 || num_sex_groups < 2 || num_distinct_ages < 2) {
    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient variation for regression (Subjects: ", num_subjects, ", Sex groups: ", num_sex_groups, ", Distinct ages: ", num_distinct_ages, ")."))
    rm(tissue_seurat); gc(); next
  }
  
  # Convert Seurat object to SingleCellExperiment for MAST
  sce_tissue <- as.SingleCellExperiment(tissue_seurat, assay = "RNA")
  
  # Ensure factors are re-leveled after subsetting/loading
  colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
  colData(sce_tissue)$sex <- droplevels(colData(sce_tissue)$sex)
  
  # Filter out subjects with only one cell within this subset for MAST stability
  subject_counts <- table(colData(sce_tissue)$donor_id)
  subjects_to_keep <- names(subject_counts[subject_counts > 1])
  
  if (length(subjects_to_keep) < 2) {
    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient subjects with more than one cell (after filtering)."))
    rm(tissue_seurat, sce_tissue); gc(); next
  }
  
  sce_tissue <- sce_tissue[, colData(sce_tissue)$donor_id %in% subjects_to_keep]
  colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
  
  if (ncol(sce_tissue) < min_cells_per_tissue) {
    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(sce_tissue), " < ", min_cells_per_tissue, ") after subject filtering."))
    rm(tissue_seurat, sce_tissue); gc(); next
  }
  
  # Filter genes to include only those expressed in a certain percentage of cells
  expressed_genes <- rowSums(assay(sce_tissue, "logcounts") > 0) / ncol(sce_tissue) > min_expressed_gene_threshold
  if (sum(expressed_genes) < min_genes_after_filter) {
    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient highly expressed genes (", sum(expressed_genes), " < ", min_genes_after_filter, ")."))
    rm(tissue_seurat, sce_tissue); gc(); next
  }
  sce_tissue_filtered <- sce_tissue[expressed_genes, ]
  
  # Ensure primerid (gene ID) and wellKey (cell ID) are explicitly defined for MAST, ensures meaningful IDs
  rowData(sce_tissue_filtered)$primerid <- rownames(sce_tissue_filtered)
  colData(sce_tissue_filtered)$wellKey <- colnames(sce_tissue_filtered)
  
  # --- Perform MAST Analysis and OmicSignature Creation ---
  message(paste0("  Running MAST for '", current_tissue_name, "' with ", nrow(sce_tissue_filtered), " genes and ", ncol(sce_tissue_filtered), " cells."))
  
  tryCatch({
    sca_mast <- SceToSingleCellAssay(sce_tissue_filtered, class = "SingleCellAssay")
    # Fit the ZLM model: gene ~ age + sex + donor_id
    zlm_obj <- zlm(~ age + sex + donor_id, sca = sca_mast, method = 'glm', ebayes = TRUE)
    results_table_mast <- MAST::as.data.frame(logFC(zlm_obj, contrasts = "age"))
    
    if (is.null(results_table_mast) || nrow(results_table_mast) == 0) {
      message(paste0("  No differential expression results found for 'age' in tissue: ", current_tissue_name))
      # Skip to cleanup and next file
      
    } else {
      # Prepare results for OmicSignature (difexp data frame)
      results_table_omic <- results_table_mast %>%
        dplyr::mutate(
          probe_id = PrimerID, # Use PrimerID as the unique identifier
          feature_name = PrimerID, # Assuming PrimerID is also the feature name (e.g., gene symbol/ENSG)
          score = `logFC`,
          p_value = `Pvalue`,
          adj_p = `FDR`
        ) %>%
        # Select the columns required for difexp, ensuring correct names
        dplyr::select(probe_id, feature_name, score, p_value, adj_p) %>%
        dplyr::mutate(
          # Define group_label for bi-directional signature (required by OmicSignature)
          group_label = ifelse(score > 0, "Increased_with_Age", "Decreased_with_Age")
        )
      
      # Define metadata for the tissue-specific signature
      metadata_tissue_sig <- OmicSignature::createMetadata(
        signature_name = paste0("Aging Signature - ", current_tissue_name),
        organism = "Homo Sapiens",
        direction_type = "bi-directional",
        phenotype = paste0("Aging in ", current_tissue_name),
        covariates = "sex, donor_id",
        platform = "Single-cell RNA-seq (cellxgene.census/Tabula Sapiens)",
        sample_type = paste0(current_tissue_name, " cells"),
        adj_p_cutoff = adj_p_cutoff,
        score_cutoff = score_cutoff,
        keywords = c("Aging", current_tissue_name, "Tabula Sapiens", "single-cell", "MAST"),
        author = "ChallengeProject2025",
        PMID = NULL, 
        year = as.numeric(format(Sys.Date(), "%Y")),
        description = paste0("Aging signature derived from Tabula Sapiens human single-cell RNA-seq data for the ", current_tissue_name, " tissue. Differential expression calculated with MAST, adjusting for sex and donor_id. Filters: min cells=",min_cells_per_tissue,", min gene expr=",min_expressed_gene_threshold*100,"%, adj.p<=",adj_p_cutoff,", |logFC|>=",score_cutoff,".")
      )
      
      # Filter significant genes for the signature (signature data frame)
      sig_genes <- results_table_omic %>%
        dplyr::filter(adj_p <= adj_p_cutoff & abs(score) >= score_cutoff) %>%
        # Select columns required for signature
        dplyr::select(probe_id, feature_name, score, group_label)
      
      if (nrow(sig_genes) == 0) {
        message(paste0("  No significant genes found for 'age' in tissue: ", current_tissue_name, " with current cutoffs (adj_p <= ", adj_p_cutoff, ", |logFC| >= ", score_cutoff, ")."))
        # Skip to cleanup
        
      } else {
        # Create the OmicSignature object
        omic_sig_tissue <- OmicSignature$new(
          metadata = metadata_tissue_sig,
          signature = sig_genes,
          difexp = results_table_omic # Store the full differential expression results
        )
        
        # Add the successfully created OmicSignature object to temporary list
        all_tissue_omicsigs[[current_tissue_name]] <- omic_sig_tissue
        message(paste0("  Successfully created and added aging signature for ", current_tissue_name, ". (", nrow(sig_genes), " significant genes)"))

        # Save individual OmicSignature object (for easier access)
        safe_tissue_name <- gsub("[^[:alnum:]_]", "_", current_tissue_name) # Standardize name (remove spaces etc)
        saveRDS(omic_sig_tissue, file = file.path(omic_signature_output_path, paste0("aging_signature_", safe_tissue_name, "_oSig.rds")))
      }
    }
  }, error = function(e) {
    message(paste0("  Error during MAST or OmicSignature creation for tissue '", current_tissue_name, "': ", e$message))
  }, finally = {
    # Ensure all large objects from this iteration are removed to free memory
    rm(list = c("tissue_seurat", "sce_tissue", "sce_tissue_filtered", "sca_mast", "zlm_obj",
                "results_table_mast", "results_table_omic", "omic_sig_tissue", "sig_genes",
                "metadata_tissue_sig") %>% Filter(exists, .))
    gc() # Force garbage collection
  }) # End tryCatch for MAST analysis
} # End loop for tissue files



# --- Save the complete OmicSignatureCollection ---
if (length(all_tissue_omicsigs) > 0) {
  message("\n--- Initializing OmicSignatureCollection ---")
  aging_signature_collection <- OmicSignatureCollection$new(
    metadata = omicsig_collection_metadata,
    OmicSigList = all_tissue_omicsigs # Pass the now populated list
  )
  
  saveRDS(aging_signature_collection, file = file.path(omic_signature_output_path, "Tabula_Sapiens_Aging_OmicSignatureCollection.rds"))
  message(paste0("\nSaved OmicSignatureCollection with ", length(aging_signature_collection$OmicSigList), " tissue signatures to '", omic_signature_output_path, "'."))
} else {
  message("\nNo aging signatures were successfully generated for any tissue and added to the collection. OmicSignatureCollection was not created.")
}

message("\nScript finished.")


