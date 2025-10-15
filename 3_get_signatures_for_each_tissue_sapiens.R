
# Load Individual Tissue Files and Perform Aging Signature Analysis

# 1. Setup and Load Libraries
library(tidyverse)        
library(Seurat)           
library(SingleCellExperiment) 
library(MAST)             
library(OmicSignature)    
library(Biobase)
library(doParallel) # For parallel backend registration
library(foreach)    # For parallelizing the outer loop
library(Matrix)     # Required for efficient sparse matrix operations


# --- Define Paths and Variables ---
data_input_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/Tabula_sapiens")
omic_signature_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/results/Tabula_sapiens")

dir.create(omic_signature_output_path, recursive = TRUE, showWarnings = FALSE)
message(paste0("Input tissue Seurat objects expected from: ", data_input_path))
message(paste0("OmicSignature results will be saved to: ", omic_signature_output_path))


# Define analysis parameters
min_cells_per_tissue <- 100         # Minimum cells required for MAST per tissue
min_expressed_gene_threshold <- 0.1 # Gene expressed in at least x% of cells
min_genes_after_filter <- 10        # Minimum number of genes to proceed with MAST
adj_p_cutoff <- 0.05                # Adjusted p-value cutoff for significant genes in signature
score_cutoff <- 0.25                # Absolute logFC cutoff for significant genes in signature



# --- PARALLELISM CONFIGURATION: Setting up nested parallelization ---
# 1. Outer loop: Multiple tissues processed concurrently.
# 2. Inner loop: MAST (zlm) uses multiple cores for each tissue.

sge_total_slots <- as.numeric(Sys.getenv("NSLOTS", unset = 1)) # Get total CPU slots from SGE job

# Re-enable 4 concurrent workers to test parallel behavior
n_concurrent_tissues <- 4

# Number of CPU cores for MAST zlm to use within each concurrent tissue analysis.
mast_cores_per_tissue <- max(1, floor(sge_total_slots / n_concurrent_tissues))

# Adjust if total slots are not perfectly divisible, or if too few slots are requested
if (sge_total_slots < n_concurrent_tissues) {
  warning(paste0("Total NSLOTS (", sge_total_slots, ") is less than n_concurrent_tissues (", n_concurrent_tissues, "). Running with n_concurrent_tissues = ", sge_total_slots, " and mast_cores_per_tissue = 1."))
  n_concurrent_tissues <- sge_total_slots
  mast_cores_per_tissue <- 1
} else if (sge_total_slots %% n_concurrent_tissues != 0) {
  warning(paste0("Total NSLOTS (", sge_total_slots, ") is not perfectly divisible by n_concurrent_tissues (", n_concurrent_tissues, "). Some internal MAST jobs might implicitly get fewer cores."))
}

message(paste0("\n--- Parallelization Configuration ---"))
message(paste0("  Total SGE slots detected: ", sge_total_slots))
message(paste0("  Number of concurrent tissue analyses (outer loop): ", n_concurrent_tissues))
message(paste0("  MAST zlm will use ", mast_cores_per_tissue, " cores per tissue analysis (inner loop)."))

# Set up parallel backend for the OUTER loop (foreach).
if (n_concurrent_tissues > 1) {
  cl <- makeCluster(n_concurrent_tissues, type = "FORK") 
  registerDoParallel(cl)
  message(paste0("  Registered parallel backend for external tissue loop with ", n_concurrent_tissues, " workers."))
} else {
  message("  Running external tissue loop in serial mode (1 worker).")
  registerDoSEQ() # Register a sequential backend if not parallelizing
}


# --- Initialize an OmicSignatureCollection Metadata ---
message("\n--- Initializing OmicSignatureCollection ---")
omicsig_collection_metadata <- list(
  collection_name = "Tabula Sapiens Human Aging Signatures - All Tissues", 
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


# --- Get list of saved tissue files ---
tissue_file_names <- list.files(data_input_path, pattern = "TabulaSapiens_.*\\.rds$", full.names = FALSE)
tissue_files <- file.path(data_input_path, tissue_file_names) 

if (length(tissue_files) == 0) {
  stop("No tissue Seurat object files found in ", data_input_path, ".")
}
message(paste0("Found ", length(tissue_files), " tissue files to analyze."))


# --- Loop through individual tissue files and perform analysis (Parallelized with foreach) ---
all_tissue_results <- foreach(file_path = tissue_files[1:4], 
                              # .export: Variables needed by each parallel worker from the main R session.
                              # Rely on auto-export for most, explicitly include complex ones.
                              .export = c("omic_signature_output_path", "min_cells_per_tissue", "min_expressed_gene_threshold", 
                                          "min_genes_after_filter", "adj_p_cutoff", "score_cutoff",
                                          "mast_cores_per_tissue"), 
                              .packages = c("tidyverse", "Seurat", "SingleCellExperiment", "MAST", "OmicSignature", "Biobase", "Matrix"),
                              .combine = 'c', 
                              .init = list(),
                              .verbose = TRUE) %dopar% {
                                
                                # Set mc.cores for MAST's internal parallelism for THIS specific worker/tissue task.
                                options(mc.cores = mast_cores_per_tissue) 
                                
                                current_tissue_name <- gsub("_", " ", gsub("TabulaSapiens_|_organ|\\.rds$", "", basename(file_path)))
                                safe_tissue_name <- gsub("[^[:alnum:]_]", "_", current_tissue_name) 
                                
                                # Define worker-specific log file for direct output
                                worker_log_file <- file.path(omic_signature_output_path, paste0("log_worker_", safe_tissue_name, ".txt"))
                                
                                # Helper function to write messages to the worker's log file
                                message_to_worker_log <- function(msg, append = TRUE) {
                                  cat(paste0(Sys.time(), " [WORKER] ", msg, "\n"), file = worker_log_file, append = append)
                                }
                                
                                message_to_worker_log(paste0("--- Starting analysis for tissue: ", current_tissue_name, " ---"), append = FALSE) # Clear file
                                
                                # Define base metadata for an OmicSignature object (to be filled or used for empty objects)
                                base_metadata_for_omicSig_template <- OmicSignature::createMetadata(
                                  signature_name = paste0("Aging Signature - ", current_tissue_name),
                                  organism = "Homo sapiens", direction_type = "bi-directional", phenotype = paste0("Aging in ", current_tissue_name),
                                  assay_type = "transcriptomics", covariates = "sex, donor_id", platform = "transcriptomics by single-cell RNA-seq",
                                  sample_type = paste0(current_tissue_name, " cells"), # Placeholder, updated later
                                  adj_p_cutoff = adj_p_cutoff, score_cutoff = score_cutoff,
                                  keywords = c("Aging", current_tissue_name, "Tabula Sapiens", "single-cell", "MAST", "human", "sapiens"),
                                  author = "ChallengeProject2025", PMID = NULL, year = as.numeric(format(Sys.Date(), "%Y")),
                                  description = paste0("Aging signature derived from Tabula Sapiens human single-cell RNA-seq data for the ", current_tissue_name, " tissue. Differential expression calculated with MAST, adjusting for sex and donor_id. Filters: min cells=",min_cells_per_tissue,", min gene expr=",min_expressed_gene_threshold*100,"%, adj.p<=",adj_p_cutoff,", |logFC|>=",score_cutoff,".")
                                )
                                
                                # Define empty data frames for OmicSignature object
                                empty_sig_df <- data.frame(
                                  probe_id = character(0), feature_name = character(0), 
                                  score = numeric(0), group_label = factor(levels=c("Increased_with_Age", "Decreased_with_Age"))
                                )
                                empty_difexp_df <- data.frame(
                                  probe_id = character(0), feature_name = character(0), 
                                  score = numeric(0), p_value = numeric(0), adj_p = numeric(0)
                                )
                                
                                # Placeholder for the final OmicSignature object for this worker
                                omic_sig_result_obj <- NULL
                                
                                # Use tryCatch for robust error handling.
                                tryCatch({
                                  
                                  message_to_worker_log(paste0("Analyzing tissue from file: ", basename(file_path), " (", current_tissue_name, ")"))
                                  message_to_worker_log(paste0("Loading tissue-specific Seurat object: ", basename(file_path)))
                                  
                                  # Load the tissue-specific Seurat object
                                  tissue_seurat <- readRDS(file_path)
                                  
                                  if ("RNA" %in% names(tissue_seurat@assays)) {
                                    current_data_matrix <- tryCatch(
                                      expr = Seurat::GetAssayData(tissue_seurat, layer = "data", assay = "RNA"), 
                                      error = function(e) {
                                        message_to_worker_log(paste0("  Warning: Could not access 'data' layer from 'RNA' assay for '", current_tissue_name, "'. Error: ", e$message))
                                        stop("Skipped_NoData_Access") # Use stop() to jump to error handler
                                      }
                                    )
                                    
                                    if (!is.null(current_data_matrix) && prod(dim(current_data_matrix)) > 0) {
                                      if (!inherits(current_data_matrix, "sparseMatrix")) {
                                        message_to_worker_log(paste0("  Converting 'data' assay (logcounts) for '", current_tissue_name, "' to sparse matrix to save memory."))
                                        tissue_seurat@assays$RNA@data <- Matrix::Matrix(current_data_matrix, sparse = TRUE)
                                      }
                                    } else {
                                      message_to_worker_log(paste0("  Warning: 'data' layer in 'RNA' assay is missing or empty for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure data is normalized before analysis."))
                                      stop("Skipped_NoData")
                                    }
                                  } else {
                                    message_to_worker_log(paste0("  Warning: 'RNA' assay not found in Seurat object for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure 'RNA' assay exists."))
                                    stop("Skipped_NoRNAAssay")
                                  }
                                  
                                  if (ncol(tissue_seurat) == 0) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "': loaded object is empty."))
                                    stop("Skipped_EmptyObject") 
                                  }
                                  
                                  message_to_worker_log(paste0("Performing pre-MAST checks for ", current_tissue_name, "."))
                                  if (ncol(tissue_seurat) < min_cells_per_tissue) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(tissue_seurat), " < ", min_cells_per_tissue, ")."))
                                    stop("Skipped_InsufficientCells")
                                  }
                                  
                                  num_subjects <- length(levels(tissue_seurat@meta.data$donor_id))
                                  num_sex_groups <- length(levels(tissue_seurat@meta.data$sex))
                                  num_distinct_ages <- length(unique(tissue_seurat@meta.data$age))
                                  
                                  if (num_subjects < 2 || num_sex_groups < 2 || num_distinct_ages < 2) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient variation for regression (Subjects: ", num_subjects, ", Sex groups: ", num_sex_groups, ", Distinct ages: ", num_distinct_ages, ")."))
                                    stop("Skipped_InsufficientVariation")
                                  }
                                  
                                  message_to_worker_log(paste0("Converting to SingleCellExperiment for ", current_tissue_name, "."))
                                  sce_tissue <- as.SingleCellExperiment(tissue_seurat, assay = "RNA")
                                  
                                  colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                  colData(sce_tissue)$sex <- droplevels(colData(sce_tissue)$sex)
                                  
                                  subject_counts <- table(colData(sce_tissue)$donor_id)
                                  subjects_to_keep <- names(subject_counts[subject_counts > 1])
                                  
                                  if (length(subjects_to_keep) < 2) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient subjects with more than one cell (after filtering)."))
                                    stop("Skipped_InsufficientSubjectsPostFilter")
                                  }
                                  
                                  sce_tissue <- sce_tissue[, colData(sce_tissue)$donor_id %in% subjects_to_keep]
                                  colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                  
                                  if (ncol(sce_tissue) < min_cells_per_tissue) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(sce_tissue), " < ", min_cells_per_tissue, ") after subject filtering."))
                                    stop("Skipped_InsufficientCellsPostFilter")
                                  }
                                  
                                  expressed_genes <- rowSums(assay(sce_tissue, "logcounts") > 0) / ncol(sce_tissue) > min_expressed_gene_threshold
                                  if (sum(expressed_genes) < min_genes_after_filter) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient highly expressed genes (", sum(expressed_genes), " < ", min_genes_after_filter, ")."))
                                    stop("Skipped_InsufficientGenes")
                                  }
                                  sce_tissue_filtered <- sce_tissue[expressed_genes, ]
                                  
                                  rowData(sce_tissue_filtered)$primerid <- rownames(sce_tissue_filtered)
                                  colData(sce_tissue_filtered)$wellKey <- colnames(sce_tissue_filtered)
                                  
                                  message_to_worker_log(paste0("  DEBUG: Dimensions of sce_tissue_filtered: ", paste(dim(sce_tissue_filtered), collapse = "x")))
                                  if (!identical(length(rownames(sce_tissue_filtered)), length(rowData(sce_tissue_filtered)$primerid))) stop("DEBUG ERROR: Rownames/primerid mismatch!")
                                  if (!identical(length(colnames(sce_tissue_filtered)), length(colData(sce_tissue_filtered)$wellKey))) stop("DEBUG ERROR: Colnames/wellKey mismatch!")
                                  if (any(nchar(rownames(sce_tissue_filtered)) == 0)) stop("DEBUG ERROR: Rownames contain empty strings!")
                                  if (any(is.na(rowData(sce_tissue_filtered)$primerid))) stop("DEBUG ERROR: primerid contains NA values!")
                                  
                                  message_to_worker_log(paste0("Running MAST for '", current_tissue_name, "' with ", nrow(sce_tissue_filtered), " genes and ", ncol(sce_tissue_filtered), " cells."))
                                  
                                  if (inherits(assay(sce_tissue_filtered, "logcounts"), "sparseMatrix")) {
                                    message_to_worker_log(paste0("  DEBUG: 'logcounts' assay in sce_tissue_filtered is sparse before SceToSingleCellAssay."))
                                  } else {
                                    message_to_worker_log(paste0("  DEBUG: 'logcounts' assay in sce_tissue_filtered is dense before SceToSingleCellAssay. This might trigger coercion."))
                                  }
                                  
                                  sca_mast <- SceToSingleCellAssay(sce_tissue_filtered, class = "SingleCellAssay")
                                  
                                  if (inherits(assay(sca_mast, "logcounts"), "sparseMatrix")) {
                                    message_to_worker_log(paste0("  DEBUG: 'logcounts' assay in sca_mast is sparse after SceToSingleCellAssay."))
                                  } else {
                                    message_to_worker_log(paste0("  DEBUG: 'logcounts' assay in sca_mast is dense after SceToSingleCellAssay. This is where dense conversion for MAST occurs."))
                                  }
                                  
                                  # Fit the ZLM model: gene ~ age + sex + donor_id.
                                  # Keep parallel=TRUE as in your original working code
                                  zlm_obj <- zlm(~ age + sex + donor_id, sca = sca_mast, method = 'glm', ebayes = TRUE, parallel = TRUE, exprs_value = 'logcounts') 
                                  
                                  summary_age_results <- summary(zlm_obj, doLRT = "age")
                                  results_table_mast_raw <- summary_age_results$datatable
                                  
                                  results_table_mast <- results_table_mast_raw %>%
                                    dplyr::filter(component == 'C' & contrast == 'age') %>% 
                                    dplyr::select(
                                      PrimerID = primerid, logFC_val = coef, Pvalue_val = `Pr(>Chisq)`     
                                    ) %>%
                                    dplyr::mutate(
                                      FDR_val = p.adjust(Pvalue_val, method = "fdr") 
                                    ) %>%
                                    dplyr::select(
                                      PrimerID = PrimerID, logFC = logFC_val, Pvalue = Pvalue_val, FDR = FDR_val
                                    )
                                  
                                  if (is.null(results_table_mast) || nrow(results_table_mast) == 0) {
                                    message_to_worker_log(paste0("  No differential expression results found for 'age' in tissue: ", current_tissue_name))
                                    stop("Skipped_NoDEResults") 
                                  } else {
                                    results_table_omic <- results_table_mast %>%
                                      dplyr::mutate(
                                        probe_id = PrimerID, feature_name = PrimerID, score = `logFC`,
                                        p_value = `Pvalue`, adj_p = `FDR`
                                      ) %>%
                                      dplyr::select(probe_id, feature_name, score, p_value, adj_p) %>%
                                      dplyr::mutate(
                                        group_label = as.factor(ifelse(score > 0, "Increased_with_Age", "Decreased_with_Age"))
                                      )
                                    
                                    found_sample_type <- NULL
                                    brenda_results_all <- OmicSignature::searchSampleType(current_tissue_name, contain_all = FALSE)
                                    if (nrow(brenda_results_all) > 0) {
                                      found_sample_type <- brenda_results_all$Name[1]
                                    } else {
                                      brenda_results_broad <- OmicSignature::searchSampleType(current_tissue_name, contain_all = FALSE)
                                      if (nrow(brenda_results_broad) > 0) {
                                        brenda_results_broad <- brenda_results_broad %>%
                                          dplyr::mutate(word_count = sapply(strsplit(Name, "\\s+"), length)) %>%
                                          dplyr::arrange(word_count)
                                        found_sample_type <- brenda_results_broad$Name[1]
                                      }
                                    }
                                    
                                    # Create metadata for the tissue-specific OmicSignature.
                                    # Always start with template, then customize
                                    metadata_for_current_omicSig <- base_metadata_for_omicSig_template
                                    if (!is.null(found_sample_type)) {
                                      metadata_for_current_omicSig$sample_type <- found_sample_type
                                    } else {
                                      message_to_worker_log(paste0("  Warning: No suitable BRENDA ontology term found for '", current_tissue_name, "'. Using '", metadata_for_current_omicSig$sample_type, "'."))
                                    }
                                    
                                    # Filter for significant genes based on defined cutoffs.
                                    sig_genes <- results_table_omic %>%
                                      dplyr::filter(adj_p <= adj_p_cutoff & abs(score) >= score_cutoff) %>%
                                      dplyr::select(probe_id, feature_name, score, group_label)
                                    
                                    if (nrow(sig_genes) == 0) {
                                      message_to_worker_log(paste0("  No significant genes found for 'age' in tissue: ", current_tissue_name, " with current cutoffs. Creating empty OmicSignature object."))
                                      
                                      # Create empty OmicSignature object and save it
                                      omic_sig_result_obj <- OmicSignature$new(
                                        metadata = metadata_for_current_omicSig, 
                                        signature = empty_sig_df, # Empty signature data frame
                                        difexp = results_table_omic # Still store the full diff exp results
                                      )
                                      saveRDS(omic_sig_result_obj, file = file.path(omic_signature_output_path, paste0("aging_signature_", safe_tissue_name, "_oSig.rds")))
                                      message_to_worker_log(paste0("  Saved empty aging signature for ", current_tissue_name, "."))
                                      
                                      return(list(setNames(list(omic_sig_result_obj), current_tissue_name))) # Return empty OmicSig object
                                    } else {
                                      # Create the OmicSignature object for the current tissue.
                                      omic_sig_result_obj <- OmicSignature$new(
                                        metadata = metadata_for_current_omicSig,
                                        signature = sig_genes,
                                        difexp = results_table_omic # Store the full differential expression results
                                      )
                                      
                                      # Save individual OmicSignature object (for easier access) within each worker.
                                      saveRDS(omic_sig_result_obj, file = file.path(omic_signature_output_path, paste0("aging_signature_", safe_tissue_name, "_oSig.rds")))
                                      message_to_worker_log(paste0("  Successfully created and added aging signature for ", current_tissue_name, ". (", nrow(sig_genes), " significant genes)"))
                                      
                                      return(list(setNames(list(omic_sig_result_obj), current_tissue_name))) # Return full OmicSig object
                                    }
                                  }
                                }, error = function(e) {
                                  error_message_text <- e$message
                                  
                                  # If it's a controlled exit, parse the status
                                  status_detail <- "Error"
                                  if (grepl("^Skipped_", error_message_text)) {
                                    status_detail <- sub("^Skipped_", "", error_message_text)
                                    message_to_worker_log(paste0("  Processing skipped for '", current_tissue_name, "' due to: ", status_detail))
                                  } else {
                                    message_to_worker_log(paste0("  ERROR: An unhandled error occurred for tissue '", current_tissue_name, "': ", error_message_text))
                                  }
                                  
                                  # Create an empty OmicSignature object with error/skipped status in description
                                  error_metadata <- base_metadata_for_omicSig_template
                                  error_metadata$description <- paste0(base_metadata_for_omicSig_template$description, " Processing ", status_detail, " due to: ", error_message_text)
                                  
                                  omic_sig_result_obj <- OmicSignature$new( 
                                    metadata = error_metadata, 
                                    signature = empty_sig_df, 
                                    difexp = empty_difexp_df
                                  )
                                  
                                  return(list(setNames(list(omic_sig_result_obj), current_tissue_name))) # Return empty OmicSig object for error
                                }, warning = function(w) {
                                  message_to_worker_log(paste0("  WARNING: for tissue '", current_tissue_name, "': ", w$message))
                                  # Warnings don't block further execution or return, they are just logged.
                                }) 
                                
                                # Final cleanup for the worker's environment.
                                rm(list=ls(all.names=TRUE)) 
                                gc(verbose = FALSE) 
                                
                                # This return should only be reached if tryCatch somehow failed to set omic_sig_result_obj
                                # and didn't jump to an error handler with stop(). It's a fallback for robustness.
                                if (is.null(omic_sig_result_obj)) {
                                  fallback_metadata <- base_metadata_for_omicSig_template
                                  fallback_metadata$description <- paste0(base_metadata_for_omicSig_template$description, " Processing status unknown due to unexpected worker termination.")
                                  omic_sig_result_obj <- OmicSignature$new(
                                    metadata = fallback_metadata,
                                    signature = empty_sig_df,
                                    difexp = empty_difexp_df
                                  )
                                  message_to_worker_log(paste0("  WARNING: Unexpected path to final return, creating fallback OmicSignature object."))
                                }
                                
                                # Ensure we log the final status based on the object's properties or metadata
                                final_log_status <- "Unknown_Status_Final"
                                if (grepl("Processing failed with error:", omic_sig_result_obj$metadata$description)) {
                                  final_log_status <- "Error"
                                } else if (grepl("Processing skipped due to:", omic_sig_result_obj$metadata$description)) {
                                  final_log_status <- "Skipped"
                                } else if (grepl("No significant genes found", omic_sig_result_obj$metadata$description) && nrow(omic_sig_result_obj$signature) == 0) {
                                  final_log_status <- "No_Significant_Genes_Found"
                                } else if (nrow(omic_sig_result_obj$signature) > 0) {
                                  final_log_status <- "Success"
                                } else if (grepl("unexpected worker termination", omic_sig_result_obj$metadata$description)) {
                                  final_log_status <- "Worker_Fallback_Error"
                                }
                                message_to_worker_log(paste0("--- Worker finished for ", current_tissue_name, " (", final_log_status, ") ---"))
                                
                                return(list(setNames(list(omic_sig_result_obj), current_tissue_name)))
                              } # End foreach loop


# --- Post-Processing of Results ---

all_tissue_omicsigs <- list()
tissue_processing_summary <- data.frame(Tissue = character(), Status = character(), stringsAsFactors = FALSE)

# Iterate through each worker's distinct result item (now each is a list(named_omicSig_object))
for (worker_res_wrapper in all_tissue_results) { 
  # Extract the named OmicSignature object
  omic_sig_obj <- worker_res_wrapper[[1]] 
  tissue_name <- names(worker_res_wrapper)[1] 
  
  # Add the OmicSignature object to the main collection list
  all_tissue_omicsigs[[tissue_name]] <- omic_sig_obj
  
  # Determine status from metadata description
  status_from_obj <- "Unknown_Status"
  if (grepl("Processing failed with error:", omic_sig_obj$metadata$description)) {
    status_from_obj <- "Error"
  } else if (grepl("Processing skipped due to:", omic_sig_obj$metadata$description)) {
    status_from_obj <- "Skipped"
  } else if (grepl("No significant genes found", omic_sig_obj$metadata$description) && nrow(omic_sig_obj$signature) == 0) {
    status_from_obj <- "No_Significant_Genes_Found"
  } else if (nrow(omic_sig_obj$signature) > 0) {
    status_from_obj <- "Success"
  } else if (grepl("unexpected worker termination", omic_sig_obj$metadata$description)) {
    status_from_obj <- "Worker_Fallback_Error"
  }
  
  # Add summary to data frame
  tissue_processing_summary <- rbind(tissue_processing_summary, 
                                     data.frame(Tissue = tissue_name, 
                                                Status = status_from_obj, 
                                                stringsAsFactors = FALSE))
}


message("\n--- Summary of Tissue Processing ---")
print(tissue_processing_summary)
message("------------------------------------\n")


if (exists("cl") && inherits(cl, "cluster")) {
  stopCluster(cl)
  message("\nStopped parallel cluster for external loop.")
}


if (length(all_tissue_omicsigs) > 0) {
  message("\n--- Creating and Saving OmicSignatureCollection ---") 
  aging_signature_collection <- OmicSignatureCollection$new(
    metadata = omicsig_collection_metadata,
    OmicSigList = all_tissue_omicsigs 
  )
  
  saveRDS(aging_signature_collection, file = file.path(omic_signature_output_path, "Tabula_Sapiens_Aging_OmicSignatureCollection.rds"))
  message(paste0("\nSaved OmicSignatureCollection with ", length(aging_signature_collection$OmicSigList), " tissue signatures to '", omic_signature_output_path, "'."))
} else {
  message("\nNo aging signatures were successfully generated for any tissue and added to the collection. OmicSignatureCollection was not created.")
}

message("\nScript finished.")

