
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

# Number of tissue analyses to run concurrently (outer loop parallelization).
# This divides the total SGE slots into independent parallel tasks.
n_concurrent_tissues <- 1

# Number of CPU cores for MAST zlm to use within each concurrent tissue analysis.
# This ensures each individual MAST run gets a dedicated set of cores.
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
# We register 'n_concurrent_tissues' workers to process different tissue files simultaneously.
if (n_concurrent_tissues > 1) {
  # Using 'FORK' type is generally more memory-efficient on Linux/Unix (like SCC)
  # as it copies the parent R session rather than creating new ones.
  cl <- makeCluster(n_concurrent_tissues, type = "FORK") 
  registerDoParallel(cl)
  message(paste0("  Registered parallel backend for external tissue loop with ", n_concurrent_tissues, " workers."))
} else {
  message("  Running external tissue loop in serial mode (1 worker).")
  registerDoSEQ() # Register a sequential backend if not parallelizing
}


# --- Initialize an OmicSignatureCollection Metadata ---
# This metadata defines the overarching collection of signatures.
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
tissue_files <- list.files(data_input_path, pattern = "TabulaSapiens_.*\\.rds$", full.names = TRUE)
if (length(tissue_files) == 0) {
  stop("No tissue Seurat object files found in ", data_input_path, ".")
}
message(paste0("Found ", length(tissue_files), " tissue files to analyze."))





# --- Loop through individual tissue files and perform analysis (Parallelized with foreach) ---
# Each iteration of this loop runs on a separate worker in parallel.
# .export: Variables needed by each parallel worker from the main R session.
# .packages: Libraries each parallel worker needs to load.
# .combine = 'list': Crucially, this keeps results from each worker as separate list items.
# .init = list(): Initializes the combined result with an empty list.
# .verbose = TRUE: Provides detailed status messages from the foreach loop itself.
all_tissue_results <- foreach(file_path = tissue_files, 
                              .export = c("data_input_path", "omic_signature_output_path",
                                          "min_cells_per_tissue", "min_expressed_gene_threshold", 
                                          "min_genes_after_filter", "adj_p_cutoff", "score_cutoff",
                                          "mast_cores_per_tissue"), 
                              .packages = c("tidyverse", "Seurat", "SingleCellExperiment", "MAST", "OmicSignature", "Biobase", "Matrix"),
                              .combine = 'list',
                              .init = list(),
                              .verbose = TRUE) %dopar% {
                                
                                # Set mc.cores for MAST's internal parallelism for THIS specific worker/tissue task.
                                options(mc.cores = mast_cores_per_tissue) 
                                
                                current_tissue_name <- gsub("_", " ", gsub("TabulaSapiens_|_organ|\\.rds$", "", basename(file_path)))
                                safe_tissue_name <- gsub("[^[:alnum:]_]", "_", current_tissue_name) # Ensure safe name for filename
                                
                                # Define worker-specific log file for immediate progress feedback
                                worker_log_file <- file.path(omic_signature_output_path, paste0("log_worker_", safe_tissue_name, ".txt"))
                                
                                # This variable will collect ALL messages/output for return to the main process
                                captured_output_vec <- character(0) 
                                
                                # Start logging to worker-specific file for immediate feedback
                                cat(paste0(Sys.time(), " --- Starting analysis for tissue: ", current_tissue_name, " ---\n"), file = worker_log_file, append = FALSE)
                                
                                # --- Set up temporary textConnection to capture ALL messages/output for return ---
                                temp_captured_conn <- textConnection("captured_output_vec", "w", local = TRUE)
                                sink(temp_captured_conn, type = "output")
                                sink(temp_captured_conn, type = "message")
                                
                                # Initialize a placeholder for the result. This will be returned at the end.
                                final_result_for_worker <- list(
                                  omicSig = NULL, 
                                  tissueName = current_tissue_name, 
                                  status = "Processing_Failed_Unspecified", # Default status if nothing else sets it
                                  messages = character(0) # Will be filled from captured_output_vec
                                )
                                
                                # Start memory profiling for this worker.
                                # Profile to a unique file for each worker, so we can analyze it if it dies.
                                mem_profile_file <- file.path(omic_signature_output_path, paste0("mem_profile_", safe_tissue_name, ".Rprof"))
                                Rprof(mem_profile_file, memory.profiling = TRUE, interval = 0.1)
                                
                                # Use tryCatch for robust error handling and warning capturing for the main logic
                                tryCatch({
                                  
                                  message(paste0("\n--- Analyzing tissue from file: ", basename(file_path), " (", current_tissue_name, ") ---"))
                                  cat(paste0(Sys.time(), " [PROGRESS] Analyzing tissue: ", current_tissue_name, "\n"), file = worker_log_file, append = TRUE)
                                  
                                  cat(paste0(Sys.time(), " [PROGRESS] Loading tissue-specific Seurat object: ", basename(file_path), "\n"), file = worker_log_file, append = TRUE)
                                  tissue_seurat <- readRDS(file_path)
                                  
                                  # --- Memory Optimization: Ensure Seurat 'data' assay is a sparse matrix ---
                                  if ("RNA" %in% names(tissue_seurat@assays)) {
                                    current_data_matrix <- tryCatch(
                                      expr = Seurat::GetAssayData(tissue_seurat, layer = "data", assay = "RNA"), 
                                      error = function(e) {
                                        message(paste0("  Warning: Could not access 'data' layer from 'RNA' assay for '", current_tissue_name, "'. Error: ", e$message))
                                        final_result_for_worker$status <<- "Skipped_NoData_Access"
                                        stop("ControlledExit") # Use a custom error to exit this block
                                      }
                                    )
                                    
                                    if (!is.null(current_data_matrix) && prod(dim(current_data_matrix)) > 0) {
                                      if (!inherits(current_data_matrix, "sparseMatrix")) {
                                        message(paste0("  Converting 'data' assay (logcounts) for '", current_tissue_name, "' to sparse matrix to save memory."))
                                        tissue_seurat@assays$RNA@data <- Matrix::Matrix(current_data_matrix, sparse = TRUE)
                                      } else {
                                        message(paste0("  'data' assay (logcounts) for '", current_tissue_name, "' is already sparse. No conversion needed."))
                                      }
                                    } else {
                                      message(paste0("  Warning: 'data' layer in 'RNA' assay is missing or empty for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure data is normalized before analysis."))
                                      final_result_for_worker$status <<- "Skipped_NoData"
                                      stop("ControlledExit") 
                                    }
                                  } else {
                                    message(paste0("  Warning: 'RNA' assay not found in Seurat object for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure 'RNA' assay exists."))
                                    final_result_for_worker$status <<- "Skipped_NoRNAAssay"
                                    stop("ControlledExit") 
                                  }
                                  
                                  # Skip if loaded object is empty
                                  if (ncol(tissue_seurat) == 0) {
                                    message(paste0("  Skipping '", current_tissue_name, "': loaded object is empty."))
                                    final_result_for_worker$status <<- "Skipped_EmptyObject"
                                    stop("ControlledExit") 
                                  }
                                  
                                  # --- Pre-MAST Data Checks and Filtering ---
                                  cat(paste0(Sys.time(), " [PROGRESS] Performing pre-MAST checks for ", current_tissue_name, ".\n"), file = worker_log_file, append = TRUE)
                                  if (ncol(tissue_seurat) < min_cells_per_tissue) {
                                    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(tissue_seurat), " < ", min_cells_per_tissue, ")."))
                                    final_result_for_worker$status <<- "Skipped_InsufficientCells"
                                    stop("ControlledExit") 
                                  }
                                  
                                  num_subjects <- length(levels(tissue_seurat@meta.data$donor_id))
                                  num_sex_groups <- length(levels(tissue_seurat@meta.data$sex))
                                  num_distinct_ages <- length(unique(tissue_seurat@meta.data$age))
                                  
                                  if (num_subjects < 2 || num_sex_groups < 2 || num_distinct_ages < 2) {
                                    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient variation for regression (Subjects: ", num_subjects, ", Sex groups: ", num_sex_groups, ", Distinct ages: ", num_distinct_ages, ")."))
                                    final_result_for_worker$status <<- "Skipped_InsufficientVariation"
                                    stop("ControlledExit") 
                                  }
                                  
                                  # Convert Seurat object to SingleCellExperiment (SCE) for MAST compatibility.
                                  cat(paste0(Sys.time(), " [PROGRESS] Converting to SingleCellExperiment for ", current_tissue_name, ".\n"), file = worker_log_file, append = TRUE)
                                  sce_tissue <- as.SingleCellExperiment(tissue_seurat, assay = "RNA")
                                  
                                  # MEMORY OPTIMIZATION: Remove the original Seurat object as it's no longer needed
                                  rm(tissue_seurat)
                                  gc(verbose = FALSE) 
                                  
                                  # Drop unused factor levels to prevent issues in downstream models.
                                  colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                  colData(sce_tissue)$sex <- droplevels(colData(sce_tissue)$sex)
                                  
                                  # Filter out subjects with only one cell for MAST stability.
                                  subject_counts <- table(colData(sce_tissue)$donor_id)
                                  subjects_to_keep <- names(subject_counts[subject_counts > 1])
                                  
                                  if (length(subjects_to_keep) < 2) {
                                    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient subjects with more than one cell (after filtering)."))
                                    final_result_for_worker$status <<- "Skipped_InsufficientSubjectsPostFilter"
                                    stop("ControlledExit") 
                                  }
                                  
                                  sce_tissue <- sce_tissue[, colData(sce_tissue)$donor_id %in% subjects_to_keep]
                                  colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                  
                                  if (ncol(sce_tissue) < min_cells_per_tissue) {
                                    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(sce_tissue), " < ", min_cells_per_tissue, ") after subject filtering."))
                                    final_result_for_worker$status <<- "Skipped_InsufficientCellsPostFilter"
                                    stop("ControlledExit") 
                                  }
                                  
                                  # Filter genes: keep only those expressed in a minimum percentage of cells.
                                  expressed_genes <- rowSums(assay(sce_tissue, "logcounts") > 0) / ncol(sce_tissue) > min_expressed_gene_threshold
                                  if (sum(expressed_genes) < min_genes_after_filter) {
                                    message(paste0("  Skipping '", current_tissue_name, "' due to insufficient highly expressed genes (", sum(expressed_genes), " < ", min_genes_after_filter, ")."))
                                    final_result_for_worker$status <<- "Skipped_InsufficientGenes"
                                    stop("ControlledExit") 
                                  }
                                  sce_tissue_filtered <- sce_tissue[expressed_genes, ]
                                  
                                  # MEMORY OPTIMIZATION: Remove the unfiltered SCE object
                                  rm(sce_tissue)
                                  gc(verbose = FALSE) 
                                  
                                  # Explicitly define primerid (gene ID) and wellKey (cell ID) for MAST.
                                  rowData(sce_tissue_filtered)$primerid <- rownames(sce_tissue_filtered)
                                  colData(sce_tissue_filtered)$wellKey <- colnames(sce_tissue_filtered)
                                  
                                  # --- DEBUGGING CHECKS: Verify data consistency before MAST processing ---
                                  message(paste0("  DEBUG: Dimensions of sce_tissue_filtered: ", paste(dim(sce_tissue_filtered), collapse = "x")))
                                  if (!identical(length(rownames(sce_tissue_filtered)), length(rowData(sce_tissue_filtered)$primerid))) stop("DEBUG ERROR: Rownames/primerid mismatch!")
                                  if (!identical(length(colnames(sce_tissue_filtered)), length(colData(sce_tissue_filtered)$wellKey))) stop("DEBUG ERROR: Colnames/wellKey mismatch!")
                                  if (any(nchar(rownames(sce_tissue_filtered)) == 0)) stop("DEBUG ERROR: Rownames contain empty strings!")
                                  if (any(is.na(rowData(sce_tissue_filtered)$primerid))) stop("DEBUG ERROR: primerid contains NA values!")
                                  # --- END DEBUGGING CHECKS ---
                                  
                                  # --- Perform MAST Analysis and OmicSignature Creation ---
                                  cat(paste0(Sys.time(), " [PROGRESS] Running MAST for '", current_tissue_name, "' with ", nrow(sce_tissue_filtered), " genes and ", ncol(sce_tissue_filtered), " cells.\n"), file = worker_log_file, append = TRUE)
                                  
                                  # Confirm sparsity before conversion to SingleCellAssay (MAST's required format).
                                  if (inherits(assay(sce_tissue_filtered, "logcounts"), "sparseMatrix")) {
                                    message(paste0("  DEBUG: 'logcounts' assay in sce_tissue_filtered is sparse before SceToSingleCellAssay."))
                                  } else {
                                    message(paste0("  DEBUG: 'logcounts' assay in sce_tissue_filtered is dense before SceToSingleCellAssay. This might trigger coercion."))
                                  }
                                  
                                  sca_mast <- SceToSingleCellAssay(sce_tissue_filtered, class = "SingleCellAssay")
                                  
                                  # MEMORY OPTIMIZATION: Remove the filtered SCE object after conversion to SCA
                                  rm(sce_tissue_filtered)
                                  gc(verbose = FALSE) 
                                  
                                  # Confirm sparsity after conversion to SingleCellAssay.
                                  if (inherits(assay(sca_mast, "logcounts"), "sparseMatrix")) {
                                    message(paste0("  DEBUG: 'logcounts' assay in sca_mast is sparse after SceToSingleCellAssay."))
                                  } else {
                                    message(paste0("  DEBUG: 'logcounts' assay in sca_mast is dense after SceToSingleCellAssay. This is where dense conversion for MAST occurs."))
                                  }
                                  
                                  # Fit the ZLM model: gene ~ age + sex + donor_id.
                                  # Explicitly setting exprs_value to 'logcounts' to prevent any ambiguity.
                                  zlm_obj <- zlm(~ age + sex + donor_id, sca = sca_mast, method = 'glm', ebayes = TRUE, parallel = TRUE, exprs_value = 'logcounts') 
                                  
                                  # MEMORY OPTIMIZATION: Remove the SingleCellAssay object now that the model is fit
                                  rm(sca_mast)
                                  gc(verbose = FALSE)
                                  
                                  # Get summary results for the 'age' coefficient using a Likelihood Ratio Test (doLRT).
                                  summary_age_results <- summary(zlm_obj, doLRT = "age")
                                  results_table_mast_raw <- summary_age_results$datatable
                                  
                                  # MEMORY OPTIMIZATION: Remove the zlm_obj now that summary is extracted
                                  rm(zlm_obj)
                                  gc(verbose = FALSE)
                                  
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
                                  
                                  # MEMORY OPTIMIZATION: Remove the raw MAST results table
                                  rm(results_table_mast_raw)
                                  gc(verbose = FALSE)
                                  
                                  # Skip if no differential expression results found for 'age'.
                                  if (is.null(results_table_mast) || nrow(results_table_mast) == 0) {
                                    message(paste0("  No differential expression results found for 'age' in tissue: ", current_tissue_name))
                                    final_result_for_worker$status <<- "Skipped_NoDEResults"
                                    stop("ControlledExit") 
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
                                    
                                    # MEMORY OPTIMIZATION: Remove results_table_mast
                                    rm(results_table_mast)
                                    gc(verbose = FALSE)
                                    
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
                                      message(paste0("  Warning: No suitable BRENDA ontology term found for '", current_tissue_name, "'. Using '", found_sample_type, "'."))
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
                                      message(paste0("  No significant genes found for 'age' in tissue: ", current_tissue_name, " with current cutoffs (adj_p <= ", adj_p_cutoff, ", |logFC| >= ", score_cutoff, ")."))
                                      # MEMORY OPTIMIZATION: Remove results_table_omic before returning NULL
                                      rm(results_table_omic)
                                      gc(verbose = FALSE)
                                      final_result_for_worker$status <<- "Skipped_NoSignificantGenes"
                                      stop("ControlledExit") 
                                    } else {
                                      # Create the OmicSignature object for the current tissue.
                                      omic_sig_object_local <- OmicSignature$new( # Use a local name
                                        metadata = metadata_tissue_sig,
                                        signature = sig_genes,
                                        difexp = results_table_omic # Store the full differential expression results
                                      )
                                      
                                      # MEMORY OPTIMIZATION: results_table_omic and sig_genes are now inside omic_sig_object
                                      rm(results_table_omic, sig_genes)
                                      gc(verbose = FALSE)
                                      
                                      # Save individual OmicSignature object (for easier access) within each worker.
                                      saveRDS(omic_sig_object_local, file = file.path(omic_signature_output_path, paste0("aging_signature_", safe_tissue_name, "_oSig.rds")))
                                      cat(paste0(Sys.time(), " [SAVED] Successfully created and added aging signature for ", current_tissue_name, ". (", nrow(omic_sig_object_local$signature), " significant genes)\n"), file = worker_log_file, append = TRUE)
                                      
                                      # Final success assignment to the placeholder
                                      final_result_for_worker$omicSig <<- omic_sig_object_local # Assign to the outer placeholder
                                      final_result_for_worker$status <<- "Success"
                                    }
                                  } 
                                }, error = function(e) {
                                  if (grepl("ControlledExit", e$message)) {
                                    # This is our controlled early exit, status is already set.
                                    # No action needed, final_result_for_worker is already set.
                                  } else {
                                    # This is a genuine, unhandled error within the tryCatch block
                                    error_message <- paste0(Sys.time(), "  ERROR: An unhandled error occurred for tissue '", current_tissue_name, "': ", e$message, "\n")
                                    message(error_message) # This message goes into `captured_output_vec`
                                    cat(error_message, file = worker_log_file, append = TRUE) # Write to worker log file for immediate visibility
                                    final_result_for_worker$status <<- "Error"
                                  }
                                }, warning = function(w) {
                                  # Just log warnings, do not stop or change status here.
                                  warning_message <- paste0(Sys.time(), "  WARNING: for tissue '", current_tissue_name, "': ", w$message, "\n")
                                  message(warning_message) # This message goes into `captured_output_vec`
                                  cat(warning_message, file = worker_log_file, append = TRUE) # Write to worker log file for immediate visibility
                                }) # End of tryCatch
                                
                                # Stop memory profiling.
                                Rprof(NULL)
                                
                                # IMPORTANT: Close the sinks for this worker's `textConnection`
                                sink(type = "message")
                                sink(type = "output")
                                close(temp_captured_conn)
                                
                                # Finalize the result to be returned using the captured output
                                final_result_for_worker$messages <- paste(captured_output_vec, collapse = "\n")
                                
                                cat(paste0(Sys.time(), " --- Worker finished for ", current_tissue_name, " (", final_result_for_worker$status, ") ---\n"), file = worker_log_file, append = TRUE)
                                
                                # Final MEMORY OPTIMIZATION for the worker process
                                rm(list=ls(all.names=TRUE)) 
                                gc(verbose = FALSE) 
                                
                                return(final_result_for_worker) # Always returns a valid, structured list
                              } # End foreach loop


# --- Post-Processing of Results ---

# Separate the OmicSignature objects from the captured messages.
all_tissue_omicsigs <- list()
all_captured_messages <- list()
tissue_processing_summary <- data.frame(Tissue = character(), Status = character(), stringsAsFactors = FALSE)

# Iterate through each worker's distinct result item
for (worker_result in all_tissue_results) {
  # Collect messages from this worker
  if (!is.null(worker_result$messages)) {
    all_captured_messages <- c(all_captured_messages, worker_result$messages)
  }
  
  # If an OmicSignature object was successfully generated, add it to the main collection list
  if (!is.null(worker_result$omicSig) && !is.null(worker_result$tissueName)) {
    all_tissue_omicsigs[[worker_result$tissueName]] <- worker_result$omicSig
  }
  
  # Add summary to data frame
  if (!is.null(worker_result$tissueName) && !is.null(worker_result$status)) {
    tissue_processing_summary <- rbind(tissue_processing_summary, 
                                       data.frame(Tissue = worker_result$tissueName, 
                                                  Status = worker_result$status, 
                                                  stringsAsFactors = FALSE))
  }
}

# Print all captured messages from workers to the main log file.
if (length(all_captured_messages) > 0) {
  message("\n--- Captured messages from parallel workers ---")
  cat(paste(all_captured_messages, collapse = "\n"), "\n") 
  message("--- End captured messages ---\n")
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
# This creates a single collection object from all successfully generated tissue signatures.
if (length(all_tissue_omicsigs) > 0) {
  message("\n--- Creating and Saving OmicSignatureCollection ---") 
  # 'all_tissue_omicsigs' is already a named list where names are tissue names.
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

