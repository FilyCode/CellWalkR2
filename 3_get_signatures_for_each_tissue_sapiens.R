
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
omic_signature_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/results/Tabula_sapiens_")

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
all_tissue_results <- foreach(file_path = tissue_files, 
                              # .export: Variables needed by each parallel worker from the main R session.
                              # Rely on auto-export for most, explicitly include complex ones.
                              .export = c("omic_signature_output_path", "min_cells_per_tissue", "min_expressed_gene_threshold", 
                                          "min_genes_after_filter", "adj_p_cutoff", "score_cutoff",
                                          "mast_cores_per_tissue"), 
                              .packages = c("tidyverse", "Seurat", "SingleCellExperiment", "MAST", "OmicSignature", "Biobase", "Matrix"),
                              .combine = 'list',
                              .init = list(),
                              .verbose = TRUE) %dopar% {
                                
                                options(mc.cores = mast_cores_per_tissue) 
                                
                                current_tissue_name <- gsub("_", " ", gsub("TabulaSapiens_|_organ|\\.rds$", "", basename(file_path)))
                                safe_tissue_name <- gsub("[^[:alnum:]_]", "_", current_tissue_name) 
                                
                                worker_log_file <- file.path(omic_signature_output_path, paste0("log_worker_", safe_tissue_name, ".txt"))
                                
                                message_to_worker_log <- function(msg, append = TRUE) {
                                  cat(paste0(Sys.time(), " [WORKER] ", msg, "\n"), file = worker_log_file, append = append)
                                }
                                
                                message_to_worker_log(paste0("--- Starting analysis for tissue: ", current_tissue_name, " ---"), append = FALSE) 
                                
                                # Initialize the result for this worker with a consistent structure
                                worker_result <- list(
                                  omicSig = NULL, 
                                  tissueName = current_tissue_name, 
                                  status = "Processing_Failed_Unspecified", 
                                  messages = character(0) # This will store worker's output for later collection
                                )
                                
                                # --- Set up temporary textConnection to capture ALL messages/output for return ---
                                captured_output_temp <- character(0) # Variable to store output
                                current_conn <- textConnection("captured_output_temp", "w", local = TRUE)
                                sink(current_conn, type = "output")
                                sink(current_conn, type = "message")
                                
                                tryCatch({
                                  
                                  message_to_worker_log(paste0("Analyzing tissue from file: ", basename(file_path), " (", current_tissue_name, ")"))
                                  message_to_worker_log(paste0("Loading tissue-specific Seurat object: ", basename(file_path)))
                                  
                                  tissue_seurat <- readRDS(file_path)
                                  
                                  if ("RNA" %in% names(tissue_seurat@assays)) {
                                    current_data_matrix <- tryCatch(
                                      expr = Seurat::GetAssayData(tissue_seurat, layer = "data", assay = "RNA"), 
                                      error = function(e) {
                                        message_to_worker_log(paste0("  Warning: Could not access 'data' layer from 'RNA' assay for '", current_tissue_name, "'. Error: ", e$message))
                                        worker_result$status <<- "Skipped_NoData_Access"
                                        stop("ControlledExit") 
                                      }
                                    )
                                    
                                    if (!is.null(current_data_matrix) && prod(dim(current_data_matrix)) > 0) {
                                      if (!inherits(current_data_matrix, "sparseMatrix")) {
                                        message_to_worker_log(paste0("  Converting 'data' assay (logcounts) for '", current_tissue_name, "' to sparse matrix to save memory."))
                                        tissue_seurat@assays$RNA@data <- Matrix::Matrix(current_data_matrix, sparse = TRUE)
                                      } else {
                                        message_to_worker_log(paste0("  'data' assay (logcounts) for '", current_tissue_name, "' is already sparse. No conversion needed."))
                                      }
                                    } else {
                                      message_to_worker_log(paste0("  Warning: 'data' layer in 'RNA' assay is missing or empty for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure data is normalized before analysis."))
                                      worker_result$status <<- "Skipped_NoData"
                                      stop("ControlledExit") 
                                    }
                                  } else {
                                    message_to_worker_log(paste0("  Warning: 'RNA' assay not found in Seurat object for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure 'RNA' assay exists."))
                                    worker_result$status <<- "Skipped_NoRNAAssay"
                                    stop("ControlledExit") 
                                  }
                                  
                                  if (ncol(tissue_seurat) == 0) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "': loaded object is empty."))
                                    worker_result$status <<- "Skipped_EmptyObject"
                                    stop("ControlledExit") 
                                  }
                                  
                                  message_to_worker_log(paste0("Performing pre-MAST checks for ", current_tissue_name, "."))
                                  if (ncol(tissue_seurat) < min_cells_per_tissue) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(tissue_seurat), " < ", min_cells_per_tissue, ")."))
                                    worker_result$status <<- "Skipped_InsufficientCells"
                                    stop("ControlledExit") 
                                  }
                                  
                                  num_subjects <- length(levels(tissue_seurat@meta.data$donor_id))
                                  num_sex_groups <- length(levels(tissue_seurat@meta.data$sex))
                                  num_distinct_ages <- length(unique(tissue_seurat@meta.data$age))
                                  
                                  if (num_subjects < 2 || num_sex_groups < 2 || num_distinct_ages < 2) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient variation for regression (Subjects: ", num_subjects, ", Sex groups: ", num_sex_groups, ", Distinct ages: ", num_distinct_ages, ")."))
                                    worker_result$status <<- "Skipped_InsufficientVariation"
                                    stop("ControlledExit") 
                                  }
                                  
                                  message_to_worker_log(paste0("Converting to SingleCellExperiment for ", current_tissue_name, "."))
                                  sce_tissue <- as.SingleCellExperiment(tissue_seurat, assay = "RNA")
                                  
                                  colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                  colData(sce_tissue)$sex <- droplevels(colData(sce_tissue)$sex)
                                  
                                  subject_counts <- table(colData(sce_tissue)$donor_id)
                                  subjects_to_keep <- names(subject_counts[subject_counts > 1])
                                  
                                  if (length(subjects_to_keep) < 2) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient subjects with more than one cell (after filtering)."))
                                    worker_result$status <<- "Skipped_InsufficientSubjectsPostFilter"
                                    stop("ControlledExit") 
                                  }
                                  
                                  sce_tissue <- sce_tissue[, colData(sce_tissue)$donor_id %in% subjects_to_keep]
                                  colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                  
                                  if (ncol(sce_tissue) < min_cells_per_tissue) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(sce_tissue), " < ", min_cells_per_tissue, ") after subject filtering."))
                                    worker_result$status <<- "Skipped_InsufficientCellsPostFilter"
                                    stop("ControlledExit") 
                                  }
                                  
                                  expressed_genes <- rowSums(assay(sce_tissue, "logcounts") > 0) / ncol(sce_tissue) > min_expressed_gene_threshold
                                  if (sum(expressed_genes) < min_genes_after_filter) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient highly expressed genes (", sum(expressed_genes), " < ", min_genes_after_filter, ")."))
                                    worker_result$status <<- "Skipped_InsufficientGenes"
                                    stop("ControlledExit") 
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
                                    worker_result$status <<- "Skipped_NoDEResults"
                                    stop("ControlledExit") 
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
                                    if (is.null(found_sample_type)) {
                                      found_sample_type <- paste0(current_tissue_name, " cells")
                                      message_to_worker_log(paste0("  Warning: No suitable BRENDA ontology term found for '", current_tissue_name, "'. Using '", found_sample_type, "'."))
                                    }
                                    
                                    metadata_tissue_sig <- OmicSignature::createMetadata(
                                      signature_name = paste0("Aging Signature - ", current_tissue_name),
                                      organism = "Homo sapiens", direction_type = "bi-directional", phenotype = paste0("Aging in ", current_tissue_name),
                                      assay_type = "transcriptomics", covariates = "sex, donor_id", platform = "transcriptomics by single-cell RNA-seq",
                                      sample_type = found_sample_type, adj_p_cutoff = adj_p_cutoff, score_cutoff = score_cutoff,
                                      keywords = c("Aging", current_tissue_name, "Tabula Sapiens", "single-cell", "MAST"),
                                      author = "ChallengeProject2025", PMID = NULL, year = as.numeric(format(Sys.Date(), "%Y")),
                                      description = paste0("Aging signature derived from Tabula Sapiens human single-cell RNA-seq data for the ", current_tissue_name, " tissue. Differential expression calculated with MAST, adjusting for sex and donor_id. Filters: min cells=",min_cells_per_tissue,", min gene expr=",min_expressed_gene_threshold*100,"%, adj.p<=",adj_p_cutoff,", |logFC|>=",score_cutoff,".")
                                    )
                                    
                                    sig_genes <- results_table_omic %>%
                                      dplyr::filter(adj_p <= adj_p_cutoff & abs(score) >= score_cutoff) %>%
                                      dplyr::select(probe_id, feature_name, score, group_label)
                                    
                                    if (nrow(sig_genes) == 0) {
                                      message_to_worker_log(paste0("  No significant genes found for 'age' in tissue: ", current_tissue_name, " with current cutoffs (adj_p <= ", adj_p_cutoff, ", |logFC| >= ", score_cutoff, ")."))
                                      worker_result$status <<- "Skipped_NoSignificantGenes"
                                      stop("ControlledExit") 
                                    } else {
                                      omic_sig_object_local <- OmicSignature$new( 
                                        metadata = metadata_tissue_sig,
                                        signature = sig_genes,
                                        difexp = results_table_omic 
                                      )
                                      
                                      saveRDS(omic_sig_object_local, file = file.path(omic_signature_output_path, paste0("aging_signature_", safe_tissue_name, "_oSig.rds")))
                                      message_to_worker_log(paste0("[SAVED] Successfully created and added aging signature for ", current_tissue_name, ". (", nrow(omic_sig_object_local$signature), " significant genes)"))
                                      
                                      worker_result$omicSig <<- omic_sig_object_local 
                                      worker_result$status <<- "Success"
                                    }
                                  } 
                                }, error = function(e) {
                                  if (grepl("ControlledExit", e$message)) {
                                    # This is our controlled early exit. Status is already set in worker_result.
                                  } else {
                                    error_message <- paste0("  ERROR: An unhandled error occurred for tissue '", current_tissue_name, "': ", e$message)
                                    message_to_worker_log(error_message) 
                                    worker_result$status <<- "Error"
                                  }
                                }, warning = function(w) {
                                  warning_message <- paste0("  WARNING: for tissue '", current_tissue_name, "': ", w$message)
                                  message_to_worker_log(warning_message) 
                                }) # End of tryCatch
                                
                                # Close the captured output connection and store messages
                                sink(type = "message")
                                sink(type = "output")
                                close(current_conn)
                                worker_result$messages <- captured_output_temp # Assign the collected output
                                
                                message_to_worker_log(paste0("--- Worker finished for ", current_tissue_name, " (", worker_result$status, ") ---"))
                                
                                rm(list=ls(all.names=TRUE)) 
                                gc(verbose = FALSE) 
                                
                                return(worker_result) # Always returns a valid, structured list
                              } # End foreach loop


# --- Post-Processing of Results ---

all_tissue_omicsigs <- list()
all_captured_messages <- list()
tissue_processing_summary <- data.frame(Tissue = character(), Status = character(), stringsAsFactors = FALSE)

# Iterate through each worker's distinct result item (now guaranteed to be a consistent list)
for (worker_res in all_tissue_results) { # worker_res is now the structured list returned by the worker
  # Collect messages from this worker
  if (!is.null(worker_res$messages) && length(worker_res$messages) > 0) {
    all_captured_messages <- c(all_captured_messages, worker_res$messages)
  }
  
  # If an OmicSignature object was successfully generated, add it to the main collection list
  if (!is.null(worker_res$omicSig) && !is.null(worker_res$tissueName)) {
    all_tissue_omicsigs[[worker_res$tissueName]] <- worker_res$omicSig
  }
  
  # Add summary to data frame
  if (!is.null(worker_res$tissueName) && !is.null(worker_res$status)) {
    tissue_processing_summary <- rbind(tissue_processing_summary, 
                                       data.frame(Tissue = worker_res$tissueName, 
                                                  Status = worker_res$status, 
                                                  stringsAsFactors = FALSE))
  }
}

# Print all captured messages from workers to the main log file.
if (length(all_captured_messages) > 0) {
  message("\n--- Captured messages from parallel workers ---")
  cat(paste(all_captured_messages, collapse = "\n"), "\n") 
  message("--- End captured messages ---\n")
} else {
  message("\n--- No captured messages from parallel workers to print ---")
}


# Print a summary of tissue processing
message("\n--- Summary of Tissue Processing ---")
print(tissue_processing_summary)
message("------------------------------------\n")


# Stop the parallel cluster for the external loop to release resources.
if (exists("cl") && inherits(cl, "cluster")) {
  stopCluster(cl)
  message("\nStopped parallel cluster for external loop.")
}

# --- Save the complete OmicSignatureCollection ---
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

