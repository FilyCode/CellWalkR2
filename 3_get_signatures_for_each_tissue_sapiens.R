
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
omic_signature_output_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/results/Tabula_sapiens_test")

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
n_concurrent_tissues <- 5

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

# --- Helper function to build a consistent result structure for each foreach worker ---
# This function prepares the return object for foreach, it does NOT contain 'return' itself.
build_structured_result_for_foreach <- function(current_tissue_name, status_val, reason_val, 
                                                initial_cells, final_cells_before_mast, 
                                                initial_genes, final_genes_before_mast, 
                                                significant_genes, omic_signature_file, 
                                                omic_sig_object) {
  
  current_tissue_summary <- data.frame(
    TissueName = current_tissue_name,
    Status = status_val,
    Reason = reason_val,
    InitialCells = initial_cells,
    FinalCellsBeforeMAST = final_cells_before_mast,
    InitialGenes = initial_genes,
    FinalGenesBeforeMAST = final_genes_before_mast,
    SignificantGenes = significant_genes,
    OmicSignatureFile = omic_signature_file,
    stringsAsFactors = FALSE
  )
  
  # Ensure omic_sig part is always a named list, even if NULL
  named_omic_sig_list <- setNames(list(omic_sig_object), current_tissue_name)
  
  return(list(omic_sig = named_omic_sig_list, summary = current_tissue_summary))
}


# --- Loop through individual tissue files and perform analysis (Parallelized with foreach) ---
all_tissue_results <- foreach(file_path = tissue_files[3:5],
                              .export = c("omic_signature_output_path", "min_cells_per_tissue", "min_expressed_gene_threshold", 
                                          "min_genes_after_filter", "adj_p_cutoff", "score_cutoff",
                                          "mast_cores_per_tissue", "build_structured_result_for_foreach"), # Export the helper
                              .packages = c("tidyverse", "Seurat", "SingleCellExperiment", "MAST", "OmicSignature", "Biobase", "Matrix"),
                              .combine = 'list', # Use 'list' to collect all results, including NULLs
                              .init = list(),
                              .verbose = TRUE) %dopar% {
                                
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
                                
                                # Initialize ALL variables that need to be captured for the final summary or OmicSignature object
                                omic_sig_object_for_return <- NULL
                                status_for_return <- "Processing" # Default to processing, will change on skips/errors/success
                                reason_for_return <- "Initialized"
                                initial_cell_count_for_return <- 0
                                final_cell_count_for_return <- 0
                                initial_gene_count_for_return <- 0
                                final_gene_count_for_return <- 0
                                significant_genes_count_for_return <- 0
                                output_file_path_for_return <- "N/A"
                                
                                # Use a flag to control whether processing continues (for early skips)
                                continue_processing <- TRUE
                                
                                # Wrap the entire processing logic in a tryCatch to handle unexpected errors gracefully
                                final_worker_output <- tryCatch({
                                  
                                  message_to_worker_log(paste0("Analyzing tissue from file: ", basename(file_path), " (", current_tissue_name, ")"))
                                  message_to_worker_log(paste0("Loading tissue-specific Seurat object: ", basename(file_path)))
                                  
                                  # Load the tissue-specific Seurat object
                                  tissue_seurat <- readRDS(file_path)
                                  
                                  # Initial assignments for current worker. Use '<-' for first assignment in this scope.
                                  initial_cell_count_for_return <- ncol(tissue_seurat)
                                  initial_gene_count_for_return <- nrow(tissue_seurat)
                                  
                                  # --- Memory Optimization: Ensure Seurat 'data' assay is a sparse matrix ---
                                  if ("RNA" %in% names(tissue_seurat@assays)) {
                                    current_data_matrix <- tryCatch(
                                      expr = Seurat::GetAssayData(tissue_seurat, layer = "data", assay = "RNA"), 
                                      error = function(e) {
                                        message_to_worker_log(paste0("  Warning: Could not access 'data' layer from 'RNA' assay for '", current_tissue_name, "'. Error: ", e$message))
                                        status_for_return <<- "Skipped" # Use <<- to assign to outer scope for 'for_return' variables
                                        reason_for_return <<- paste0("Error accessing RNA data layer: ", e$message)
                                        continue_processing <<- FALSE
                                        return(NULL) # Indicate failure for the inner tryCatch, but outer tryCatch handles final return
                                      }
                                    )
                                    
                                    if (!is.null(current_data_matrix) && prod(dim(current_data_matrix)) == 0) {
                                      message_to_worker_log(paste0("  Warning: 'data' layer in 'RNA' assay is missing or empty or access error for '", current_tissue_name, "'. Skipping analysis."))
                                      status_for_return <<- "Skipped"
                                      reason_for_return <<- "RNA 'data' layer missing/empty or access error"
                                      continue_processing <<- FALSE
                                    } else if (inherits(current_data_matrix, "sparseMatrix")) {
                                      message_to_worker_log(paste0("  'data' assay (logcounts) for '", current_tissue_name, "' is already sparse. No conversion needed."))
                                    } else if (continue_processing) { # Only convert if still processing
                                      message_to_worker_log(paste0("  Converting 'data' assay (logcounts) for '", current_tissue_name, "' to sparse matrix to save memory."))
                                      tissue_seurat@assays$RNA@data <- Matrix::Matrix(current_data_matrix, sparse = TRUE)
                                    }
                                  } else {
                                    message_to_worker_log(paste0("  Warning: 'RNA' assay not found in Seurat object for '", current_tissue_name, "'. Skipping analysis."))
                                    status_for_return <<- "Skipped"
                                    reason_for_return <<- "RNA assay not found in Seurat object"
                                    continue_processing <<- FALSE
                                  }
                                  
                                  # Skip if loaded object is empty
                                  if (continue_processing && ncol(tissue_seurat) == 0) {
                                    message_to_worker_log(paste0("  Skipping '", current_tissue_name, "': loaded object is empty."))
                                    status_for_return <<- "Skipped"
                                    reason_for_return <<- "Loaded Seurat object is empty"
                                    continue_processing <<- FALSE
                                  }
                                  
                                  # --- Pre-MAST Data Checks and Filtering ---
                                  if (continue_processing) {
                                    message_to_worker_log(paste0("Performing pre-MAST checks for ", current_tissue_name, "."))
                                    if (ncol(tissue_seurat) < min_cells_per_tissue) {
                                      message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(tissue_seurat), " < ", min_cells_per_tissue, ")."))
                                      status_for_return <<- "Skipped"
                                      reason_for_return <<- paste0("Insufficient cells (", ncol(tissue_seurat), " < ", min_cells_per_tissue, ")")
                                      continue_processing <<- FALSE
                                    }
                                  }
                                  
                                  if (continue_processing) {
                                    num_subjects <- length(levels(tissue_seurat@meta.data$donor_id))
                                    num_sex_groups <- length(levels(tissue_seurat@meta.data$sex))
                                    num_distinct_ages <- length(unique(tissue_seurat@meta.data$age))
                                    
                                    if (num_subjects < 2 || num_sex_groups < 2 || num_distinct_ages < 3) {
                                      message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient variation for regression (Subjects: ", num_subjects, ", Sex groups: ", num_sex_groups, ", Distinct ages: ", num_distinct_ages, ")."))
                                      status_for_return <<- "Skipped"
                                      reason_for_return <<- paste0("Insufficient variation for regression (Subjects: ", num_subjects, ", Sex groups: ", num_sex_groups, ", Distinct ages: ", num_distinct_ages, ")")
                                      continue_processing <<- FALSE
                                    }
                                  }
                                  
                                  # --- Core Processing (only if 'continue_processing' is still TRUE) ---
                                  if (continue_processing) {
                                    # Convert Seurat object to SingleCellExperiment (SCE) for MAST compatibility.
                                    message_to_worker_log(paste0("Converting to SingleCellExperiment for ", current_tissue_name, "."))
                                    sce_tissue <- as.SingleCellExperiment(tissue_seurat, assay = "RNA")
                                    rm(tissue_seurat); gc(verbose = FALSE) # Clear Seurat object after conversion
                                    
                                    # Drop unused factor levels to prevent issues in downstream models.
                                    colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                    colData(sce_tissue)$sex <- droplevels(colData(sce_tissue)$sex)
                                    
                                    # Filter out subjects with only one cell for MAST stability.
                                    subject_counts <- table(colData(sce_tissue)$donor_id)
                                    subjects_to_keep <- names(subject_counts[subject_counts > 1])
                                    
                                    if (length(subjects_to_keep) < 2) {
                                      message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient subjects with more than one cell (after filtering)."))
                                      status_for_return <<- "Skipped"
                                      reason_for_return <<- "Insufficient subjects with more than one cell after filtering"
                                      continue_processing <<- FALSE
                                    } else { # Only proceed if subjects are sufficient
                                      sce_tissue <- sce_tissue[, colData(sce_tissue)$donor_id %in% subjects_to_keep]
                                      colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                      
                                      if (ncol(sce_tissue) < min_cells_per_tissue) {
                                        message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(sce_tissue), " < ", min_cells_per_tissue, ") after subject filtering."))
                                        status_for_return <<- "Skipped"
                                        reason_for_return <<- paste0("Insufficient cells (", ncol(sce_tissue), " < ", min_cells_per_tissue, ") after subject filtering")
                                        continue_processing <<- FALSE
                                      } else { # Only proceed if cells are sufficient
                                        # Filter genes: keep only those expressed in a minimum percentage of cells.
                                        expressed_genes <- rowSums(assay(sce_tissue, "logcounts") > 0) / ncol(sce_tissue) > min_expressed_gene_threshold
                                        if (sum(expressed_genes) < min_genes_after_filter) {
                                          message_to_worker_log(paste0("  Skipping '", current_tissue_name, "' due to insufficient highly expressed genes (", sum(expressed_genes), " < ", min_genes_after_filter, ")."))
                                          status_for_return <<- "Skipped"
                                          reason_for_return <<- paste0("Insufficient highly expressed genes (", sum(expressed_genes), " < ", min_genes_after_filter, ")")
                                          continue_processing <<- FALSE
                                        } else { # All filters passed, proceed with MAST preparation
                                          sce_tissue_filtered <- sce_tissue[expressed_genes, ]
                                          rm(sce_tissue); gc(verbose = FALSE) # Clear original SCE object after filtering
                                          
                                          final_cell_count_for_return <<- ncol(sce_tissue_filtered) # Update with 'final' counts
                                          final_gene_count_for_return <<- nrow(sce_tissue_filtered) # Update with 'final' counts
                                          
                                          # Explicitly define primerid (gene ID) and wellKey (cell ID) for MAST.
                                          rowData(sce_tissue_filtered)$primerid <- rownames(sce_tissue_filtered)
                                          colData(sce_tissue_filtered)$wellKey <- colnames(sce_tissue_filtered)
                                          
                                          # --- Perform MAST Analysis and OmicSignature Creation ---
                                          message_to_worker_log(paste0("Running MAST for '", current_tissue_name, "' with ", nrow(sce_tissue_filtered), " genes and ", ncol(sce_tissue_filtered), " cells."))
                                          
                                          sca_mast <- SceToSingleCellAssay(sce_tissue_filtered, class = "SingleCellAssay")
                                          rm(sce_tissue_filtered); gc(verbose = FALSE) # Clear filtered SCE object after conversion to SCA
                                          
                                          # Fit the ZLM model: gene ~ age + sex + donor_id.
                                          if (safe_tissue_name %in% c('ovary', 'prostate_gland', 'testis')) {
                                            zlm_obj <- zlm(~ age + donor_id, sca = sca_mast, method = 'glm', ebayes = TRUE, parallel = TRUE, exprs_value = 'logcounts') 
                                          } else {
                                            zlm_obj <- zlm(~ age + sex + donor_id, sca = sca_mast, method = 'glm', ebayes = TRUE, parallel = TRUE, exprs_value = 'logcounts') 
                                          }
                                          rm(sca_mast); gc(verbose = FALSE) # Clear SCA object after ZLM model creation
                                          
                                          # Get summary results for the 'age' coefficient using a Likelihood Ratio Test (doLRT).
                                          summary_age_results <- summary(zlm_obj, doLRT = "age")
                                          results_table_mast_raw <- summary_age_results$datatable
                                          rm(zlm_obj, summary_age_results); gc(verbose = FALSE) # Clear ZLM object and summary after extracting datatable
                                          
                                          # Filter for the 'age' contrast and calculate FDR.
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
                                          rm(results_table_mast_raw); gc(verbose = FALSE) # Clear raw MAST results after filtering
                                          
                                          # Skip if no differential expression results found for 'age'.
                                          if (is.null(results_table_mast) || nrow(results_table_mast) == 0) {
                                            message_to_worker_log(paste0("  No differential expression results found for 'age' in tissue: ", current_tissue_name))
                                            status_for_return <<- "Skipped"
                                            reason_for_return <<- "No differential expression results found for 'age'"
                                            continue_processing <<- FALSE # Set flag, but let flow continue to final return
                                          } else {
                                            # Prepare results for OmicSignature object (difexp data frame).
                                            results_table_omic <- results_table_mast %>%
                                              dplyr::mutate(
                                                probe_id = PrimerID, feature_name = PrimerID, score = `logFC`,
                                                p_value = `Pvalue`, adj_p = `FDR`
                                              ) %>%
                                              dplyr::select(probe_id, feature_name, score, p_value, adj_p) %>%
                                              dplyr::mutate(
                                                group_label = as.factor(ifelse(score > 0, "Increased_with_Age", "Decreased_with_Age"))
                                              )
                                            rm(results_table_mast); gc(verbose = FALSE) # Clear results_table_mast after conversion to Omic format
                                            
                                            # --- Automated BRENDA ontology lookup for sample_type ---
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
                                            
                                            # Create metadata for the tissue-specific OmicSignature.
                                            metadata_tissue_sig <- OmicSignature::createMetadata(
                                              signature_name = paste0("Aging Signature - ", current_tissue_name),
                                              organism = "Homo sapiens", direction_type = "bi-directional", phenotype = paste0("Aging in ", current_tissue_name),
                                              assay_type = "transcriptomics", covariates = "sex, donor_id", platform = "transcriptomics by single-cell RNA-seq",
                                              sample_type = found_sample_type, adj_p_cutoff = adj_p_cutoff, score_cutoff = score_cutoff,
                                              keywords = c("Aging", current_tissue_name, "Tabula Sapiens", "single-cell", "MAST"),
                                              author = "ChallengeProject2025", PMID = NULL, year = as.numeric(format(Sys.Date(), "%Y")),
                                              description = paste0("Aging signature derived from Tabula Sapiens human single-cell RNA-seq data for the ", current_tissue_name, " tissue. Differential expression calculated with MAST, adjusting for sex and donor_id. Filters: min cells=",min_cells_per_tissue,", min gene expr=",min_expressed_gene_threshold*100,"%, adj.p<=",adj_p_cutoff,", |logFC|>=",score_cutoff,".")
                                            )
                                            
                                            # Filter for significant genes based on defined cutoffs.
                                            sig_genes <- results_table_omic %>%
                                              dplyr::filter(adj_p <= adj_p_cutoff & abs(score) >= score_cutoff) %>%
                                              dplyr::select(probe_id, feature_name, score, group_label)
                                            
                                            # Skip if no significant genes found after filtering.
                                            if (nrow(sig_genes) == 0) {
                                              message_to_worker_log(paste0("  No significant genes found for 'age' in tissue: ", current_tissue_name, " with current cutoffs (adj_p <= ", adj_p_cutoff, ", |logFC| >= ", score_cutoff, "). No OmicSignature object created."))
                                              
                                              status_for_return <<- "Skipped (No Sig Genes)" 
                                              reason_for_return <<- paste0("No significant genes found with current cutoffs (adj_p <= ", adj_p_cutoff, ", |logFC| >= ", score_cutoff, "). OmicSignature object not created.")
                                              significant_genes_count_for_return <<- 0
                                              output_file_path_for_return <<- "N/A"
                                              omic_sig_object_for_return <<- NULL 
                                              rm(results_table_omic, sig_genes); gc(verbose = FALSE) 
                                            } else {
                                              # Create the OmicSignature object for the current tissue.
                                              omic_sig_tissue <- OmicSignature$new(
                                                metadata = metadata_tissue_sig,
                                                signature = sig_genes,
                                                difexp = results_table_omic 
                                              )
                                              
                                              # Save individual OmicSignature object (for easier access) within each worker.
                                              current_output_file <- file.path(omic_signature_output_path, paste0("aging_signature_", safe_tissue_name, "_oSig.rds"))
                                              saveRDS(omic_sig_tissue, file = current_output_file)
                                              message_to_worker_log(paste0("  Successfully created and added aging signature for ", current_tissue_name, ". (", nrow(sig_genes), " significant genes)"))
                                              
                                              status_for_return <<- "Success"
                                              reason_for_return <<- "OmicSignature created successfully"
                                              significant_genes_count_for_return <<- nrow(sig_genes)
                                              output_file_path_for_return <<- current_output_file
                                              omic_sig_object_for_return <<- omic_sig_tissue
                                            }
                                          }
                                        }
                                      }
                                    }
                                  } # End of if (continue_processing) for core logic
                                  
                                  # THIS MUST BE THE FINAL EXPRESSION OF THE tryCatch's expr block
                                  # All outcomes (success, various skips) flow to this consistent call.
                                  build_structured_result_for_foreach(
                                    current_tissue_name = current_tissue_name,
                                    status_val = status_for_return,
                                    reason_val = reason_for_return,
                                    initial_cells = initial_cell_count_for_return,
                                    final_cells_before_mast = final_cell_count_for_return,
                                    initial_genes = initial_gene_count_for_return,
                                    final_genes_before_mast = final_gene_count_for_return,
                                    significant_genes = significant_genes_count_for_return,
                                    omic_signature_file = output_file_path_for_return,
                                    omic_sig_object = omic_sig_object_for_return
                                  )
                                  
                                }, error = function(e) {
                                  # This error handler catches any CRITICAL, unexpected errors during processing.
                                  message_to_worker_log(paste0("  CRITICAL ERROR during processing for tissue '", current_tissue_name, "': ", e$message))
                                  
                                  # Ensure all variables used in build_structured_result_for_foreach are available.
                                  # For initial/final counts, use the values captured before the error occurred.
                                  build_structured_result_for_foreach(
                                    current_tissue_name = current_tissue_name,
                                    status_val = "Error",
                                    reason_val = paste0("CRITICAL ERROR during processing: ", e$message),
                                    initial_cells = initial_cell_count_for_return,
                                    final_cells_before_mast = final_cell_count_for_return,
                                    initial_genes = initial_gene_count_for_return,
                                    final_genes_before_mast = final_gene_count_for_return,
                                    significant_genes = 0, # No significant genes on error
                                    omic_signature_file = "N/A", # No file saved on error
                                    omic_sig_object = NULL # No OmicSignature object on error
                                  )
                                })
                                
                                # Explicitly return the result from the tryCatch block for the foreach combiner
                                return(final_worker_output)
                              } # End foreach loop


# --- Post-Processing of Results ---

# Stop the parallel cluster for the external loop to release resources.
if (exists("cl") && inherits(cl, "cluster")) {
  stopCluster(cl)
  message("\nStopped parallel cluster for external loop.")
}

# Now, 'all_tissue_results' is a list of lists (each inner list is what a worker returned).
# We need to extract the omic_sig and summary parts from each.

# 1. Extract all 'omic_sig' components. This will be a list of named lists (or NULLs).
all_omic_sigs_named_lists <- lapply(all_tissue_results, `[[`, "omic_sig")

# 2. Extract all 'summary' data frames. This will be a list of single-row data frames.
all_summaries_dfs <- lapply(all_tissue_results, `[[`, "summary")

# 3. Combine the list of single-element named lists for OmicSignature objects into a single named list.
# `do.call(c, ...)` flattens a list of lists into a single list.
all_tissue_omicsigs_for_collection <- do.call(c, all_omic_sigs_named_lists)

# 4. Filter out any NULL entries from the combined OmicSignature list (for skipped/errored tissues).
all_tissue_omicsigs_filtered <- Filter(Negate(is.null), all_tissue_omicsigs_for_collection)

# 5. Combine all tissue processing summaries into a single data frame.
# Ensure that if all_summaries_dfs is empty, we still get an empty dataframe with correct columns
if (length(all_summaries_dfs) > 0) {
  tissue_processing_summary_df <- do.call(rbind, all_summaries_dfs)
} else {
  # Create an empty dataframe with the expected columns if no summaries were generated
  tissue_processing_summary_df <- data.frame(
    TissueName = character(0),
    Status = character(0),
    Reason = character(0),
    InitialCells = numeric(0),
    FinalCellsBeforeMAST = numeric(0),
    InitialGenes = numeric(0),
    FinalGenesBeforeMAST = numeric(0),
    SignificantGenes = numeric(0),
    OmicSignatureFile = character(0),
    stringsAsFactors = FALSE
  )
}

# --- Intermediate Results Debugging ---
message("\n--- Intermediate Results Debugging ---")
message("Contents of all_tissue_results structure (list of worker returns):")
str(all_tissue_results, max.level = 2) 
message("\nContents of all_omic_sigs_named_lists structure (list of named OmicSig lists/NULLs):")
str(all_omic_sigs_named_lists, max.level = 2)
message("\nContents of all_summaries_dfs structure (list of summary dataframes):")
str(all_summaries_dfs, max.level = 2)
message("\nContents of all_tissue_omicsigs_for_collection (flattened list, including NULLs):")
str(all_tissue_omicsigs_for_collection, max.level = 2)
message("\nContents of all_tissue_omicsigs_filtered (flattened, named, and non-NULL):")
str(all_tissue_omicsigs_filtered, max.level = 2)
message("\nContents of tissue_processing_summary_df (combined dataframe):")
print(tissue_processing_summary_df) 
message("\n--- End Intermediate Results Debugging ---\n")

# Save the comprehensive tissue processing summary
summary_output_file <- file.path(omic_signature_output_path, "Tabula_Sapiens_Tissue_Processing_Summary.csv")
write.csv(tissue_processing_summary_df, file = summary_output_file, row.names = FALSE)
message(paste0("\nSaved tissue processing summary to '", summary_output_file, "'."))


# --- Save the complete OmicSignatureCollection ---
# This creates a single collection object from all successfully generated tissue signatures.
if (length(all_tissue_omicsigs_filtered) > 0) {
  message("\n--- Creating and Saving OmicSignatureCollection ---") 
  
  aging_signature_collection <- OmicSignatureCollection$new(
    metadata = omicsig_collection_metadata,
    OmicSigList = all_tissue_omicsigs_filtered # Use the correctly named and filtered list
  )
  
  saveRDS(aging_signature_collection, file = file.path(omic_signature_output_path, "Tabula_Sapiens_Aging_OmicSignatureCollection.rds"))
  message(paste0("\nSaved OmicSignatureCollection with ", length(aging_signature_collection$OmicSigList), " tissue signatures to '", omic_signature_output_path, "'."))
} else {
  message("\nNo aging signatures were successfully generated for any tissue and added to the collection. OmicSignatureCollection was not created.")
}

message("\nScript finished.")

