
# Load Individual Tissue Files and Perform Aging Signature Analysis

# 1. Setup and Load Libraries
library(tidyverse)        
library(Seurat)           
library(SingleCellExperiment) 
library(MAST)             
library(OmicSignature)    
library(Biobase)
library(doParallel)
library(Matrix)


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
  options(mc.cores = num_cores) 
  registerDoParallel(cores = num_cores)
  message(paste0("\n  Registered parallel backend for MAST zlm with ", num_cores, " cores."))
  message(paste0("  Explicitly set options(mc.cores = ", num_cores, ") for MAST::zlm."))
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
  current_tissue_name <- gsub("_", " ", gsub("TabulaSapiens_|_organ|\\.rds$", "", basename(file_path))) # get tissue name out of file name
  message(paste0("\n--- Analyzing tissue from file: ", basename(file_path), " (", current_tissue_name, ") ---"))
  
  # Load the tissue-specific Seurat object
  tissue_seurat <- readRDS(file_path)
  
  # Explicitly ensure the 'data' slot (log-normalized counts, typically used by MAST) is sparse 
  # before converting to SingleCellExperiment. This prevents large memory allocations.
  if ("RNA" %in% names(tissue_seurat@assays)) {
    # Attempt to retrieve the 'data' slot (log-normalized data) safely
    current_data_matrix <- tryCatch(
      expr = Seurat::GetAssayData(tissue_seurat, slot = "data", assay = "RNA"),
      error = function(e) {
        # If GetAssayData errors, it means the slot is likely missing or inaccessible.
        message(paste0("  Warning: Could not access 'data' slot from 'RNA' assay for '", current_tissue_name, "'. This may indicate missing normalized data. Error: ", e$message))
        return(NULL) # Return NULL to indicate failure to retrieve data
      }
    )
    
    # Check if data matrix was successfully retrieved and is not empty
    if (!is.null(current_data_matrix) && prod(dim(current_data_matrix)) > 0) {
      # Check if it's currently a dense matrix (i.e., not inheriting from a sparse matrix class)
      if (!inherits(current_data_matrix, "sparseMatrix")) {
        message(paste0("  Converting 'data' assay (logcounts) for '", current_tissue_name, "' to sparse matrix to save memory."))
        # Assign the sparse version back to the Seurat object's data slot
        tissue_seurat@assays$RNA@data <- Matrix::Matrix(current_data_matrix, sparse = TRUE)
      } else {
        message(paste0("  'data' assay (logcounts) for '", current_tissue_name, "' is already sparse. No conversion needed."))
      }
    } else {
      # This 'else' block means current_data_matrix was NULL (from tryCatch error) or empty.
      # This is the correct place to issue a warning if 'data' is truly missing or empty.
      message(paste0("  Warning: 'data' slot in 'RNA' assay is missing or empty for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure data is normalized before analysis."))
    }
  } else {
    message(paste0("  Warning: 'RNA' assay not found in Seurat object for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure 'RNA' assay exists."))
  }
  
  
  # Ensure the object is not empty after loading (shouldn't be if saved correctly)
  if (ncol(tissue_seurat) == 0) {
    message(paste0("  Skipping '", current_tissue_name, "': loaded object is empty."))
    rm(tissue_seurat); gc(); next
  }
  
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
  
  # --- DEBUGGING CHECKS FOR DIMNAMES ERROR ---
  # These checks help ensure consistency before MAST conversion
  message(paste0("  DEBUG: Dimensions of sce_tissue_filtered: ", paste(dim(sce_tissue_filtered), collapse = "x")))
  message(paste0("  DEBUG: Length of rownames(sce_tissue_filtered): ", length(rownames(sce_tissue_filtered))))
  message(paste0("  DEBUG: Length of rowData(sce_tissue_filtered)$primerid: ", length(rowData(sce_tissue_filtered)$primerid)))
  message(paste0("  DEBUG: Length of colnames(sce_tissue_filtered): ", length(colnames(sce_tissue_filtered))))
  message(paste0("  DEBUG: Length of colData(sce_tissue_filtered)$wellKey: ", length(colData(sce_tissue_filtered)$wellKey)))
  
  if (!identical(length(rownames(sce_tissue_filtered)), length(rowData(sce_tissue_filtered)$primerid))) {
    stop(paste0("DEBUG ERROR (", current_tissue_name, "): Rownames length (", length(rownames(sce_tissue_filtered)), ") does not match primerid length (", length(rowData(sce_tissue_filtered)$primerid), ")!"))
  }
  if (!identical(length(colnames(sce_tissue_filtered)), length(colData(sce_tissue_filtered)$wellKey))) {
    stop(paste0("DEBUG ERROR (", current_tissue_name, "): Colnames length (", length(colnames(sce_tissue_filtered)), ") does not match wellKey length (", length(colData(sce_tissue_filtered)$wellKey), ")!"))
  }
  if (any(nchar(rownames(sce_tissue_filtered)) == 0)) {
    stop(paste0("DEBUG ERROR (", current_tissue_name, "): Rownames contain empty strings!"))
  }
  if (any(is.na(rowData(sce_tissue_filtered)$primerid))) {
    stop(paste0("DEBUG ERROR (", current_tissue_name, "): primerid contains NA values!"))
  }
  # --- END DEBUGGING CHECKS ---
  
  # --- Perform MAST Analysis and OmicSignature Creation ---
  message(paste0("  Running MAST for '", current_tissue_name, "' with ", nrow(sce_tissue_filtered), " genes and ", ncol(sce_tissue_filtered), " cells."))
  
  tryCatch({
    sca_mast <- SceToSingleCellAssay(sce_tissue_filtered, class = "SingleCellAssay")
    # Fit the ZLM model: gene ~ age + sex + donor_id
    zlm_obj <- zlm(~ age + sex + donor_id, sca = sca_mast, method = 'glm', ebayes = TRUE, parallel = TRUE)
    
    # Get summary results for the 'age' coefficient using doLRT
    # This performs a Likelihood Ratio Test for the 'age' variable
    summary_age_results <- summary(zlm_obj, doLRT = "age")
    
    # Extract the datatable, which contains the results for each coefficient
    # This table has columns like 'primerid', 'component', 'coef', 'logFC', 'Pvalue', 'FDR'
    results_table_mast_raw <- summary_age_results$datatable
    
    # Filter for the 'age' coefficient from the continuous component ('C')
    # Then select and rename columns to match what OmicSignature expects
    # The 'coef' column in results_table_mast_raw contains the logFC for each contrast.
    # The 'Pr(>Chisq)' column is the p-value from the LRT.
    results_table_mast <- results_table_mast_raw %>%
      dplyr::filter(component == 'C' & contrast == 'age') %>% # Filter for the 'age' contrast
      dplyr::select(
        PrimerID = primerid,          # Use 'primerid' as the gene ID
        logFC_val = coef,             # 'coef' column contains the logFC value
        Pvalue_val = `Pr(>Chisq)`     # 'Pr(>Chisq)' column contains the p-value
      ) %>%
      dplyr::mutate(
        FDR_val = p.adjust(Pvalue_val, method = "fdr") # Calculate FDR from Pvalue
      ) %>%
      dplyr::select(
        PrimerID = PrimerID,
        logFC = logFC_val,
        Pvalue = Pvalue_val,
        FDR = FDR_val
      )
    
    # Now, the 'results_table_mast' will be correctly populated for the 'if' condition check
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
          group_label = as.factor(ifelse(score > 0, "Increased_with_Age", "Decreased_with_Age"))
        )
      
      
      # --- Automated BRENDA ontology lookup for sample_type ---
      found_sample_type <- NULL
      
      # Perform a broad search for the current tissue name (case-insensitive)
      brenda_results_all <- OmicSignature::searchSampleType(current_tissue_name, contain_all = FALSE)
      
      if (nrow(brenda_results_all) > 0) {
        # If a direct "<tissue_name> tissue" match is found, use it.
        found_sample_type <- brenda_results_all$Name[1]
      } else {
        # 2. If no direct match, perform a broader search for the tissue name.
        brenda_results_broad <- OmicSignature::searchSampleType(current_tissue_name, contain_all = FALSE)
        
        if (nrow(brenda_results_broad) > 0) {
          # If broad results are found, prioritize by the fewest words.
          # We'll add a word count column and sort by it.
          brenda_results_broad <- brenda_results_broad %>%
            dplyr::mutate(word_count = sapply(strsplit(Name, "\\s+"), length)) %>%
            dplyr::arrange(word_count) # Sort by least words first
          
          found_sample_type <- brenda_results_broad$Name[1] # Pick the one with the fewest words
        }
      }
      
      
      if (is.null(found_sample_type)) {
        # Fallback if no BRENDA terms were found at all
        found_sample_type <- paste0(current_tissue_name, " cells")
        message(paste0("  Warning: No suitable BRENDA ontology term found for '", current_tissue_name, "'. Using '", found_sample_type, "'."))
      }
      
      
      # Define metadata for the tissue-specific signature
      metadata_tissue_sig <- OmicSignature::createMetadata(
        signature_name = paste0("Aging Signature - ", current_tissue_name),
        organism = "Homo sapiens",
        direction_type = "bi-directional",
        phenotype = paste0("Aging in ", current_tissue_name),
        assay_type = "transcriptomics",
        covariates = "sex, donor_id",
        platform = "transcriptomics by single-cell RNA-seq",
        sample_type = found_sample_type,
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


