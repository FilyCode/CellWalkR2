
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

# Re-enable concurrent workers to test parallel behavior
n_concurrent_tissues <- 2

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
all_tissue_results <- foreach(file_path = tissue_files,
                              .export = c("omic_signature_output_path", "min_cells_per_tissue", "min_expressed_gene_threshold", 
                                          "min_genes_after_filter", "adj_p_cutoff", "score_cutoff",
                                          "mast_cores_per_tissue", "build_structured_result_for_foreach"), # Export the helper
                              .packages = c("tidyverse", "Seurat", "SingleCellExperiment", "MAST", "OmicSignature", "Biobase", "Matrix"),
                              .combine = 'c',
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
                                
                                # Initialize ALL variables LOCAL to this worker. These will be returned by process_tissue_inner.
                                omic_sig_object_out <- NULL
                                status_out <- "Processing" 
                                reason_out <- "Initialized"
                                initial_cell_count_out <- 0
                                final_cell_count_out <- 0
                                initial_gene_count_out <- 0
                                final_gene_count_out <- 0
                                significant_genes_count_out <- 0
                                output_file_path_out <- "N/A"
                                
                                # Inner function to encapsulate processing logic and return all outcomes
                                process_tissue_inner <- function(current_tissue_name, safe_tissue_name, file_path) {
                                  # Declare variables that will hold the results for this worker, default values
                                  omic_sig_object <- NULL
                                  status <- "Skipped" 
                                  reason <- "Early exit"
                                  initial_cell_count <- 0
                                  final_cell_count <- 0
                                  initial_gene_count <- 0
                                  final_gene_count <- 0
                                  significant_genes_count <- 0
                                  output_file_path <- "N/A"
                                  
                                  tryCatch({
                                    message_to_worker_log(paste0("Analyzing tissue from file: ", basename(file_path), " (", current_tissue_name, ")"))
                                    message_to_worker_log(paste0("Loading tissue-specific Seurat object: ", basename(file_path)))
                                    
                                    tissue_seurat <- readRDS(file_path)
                                    
                                    initial_cell_count <- ncol(tissue_seurat)
                                    initial_gene_count <- nrow(tissue_seurat)
                                    
                                    if (!("RNA" %in% names(tissue_seurat@assays))) {
                                      reason <- "RNA assay not found in Seurat object"
                                      stop(reason) # Use stop to jump to error handler and set status
                                    }
                                    
                                    current_data_matrix <- tryCatch(
                                      expr = Seurat::GetAssayData(tissue_seurat, layer = "data", assay = "RNA"), 
                                      error = function(e) {
                                        stop(paste0("Error accessing RNA data layer: ", e$message))
                                      }
                                    )
                                    
                                    if (is.null(current_data_matrix) || prod(dim(current_data_matrix)) == 0) {
                                      stop("RNA 'data' layer missing/empty or access error")
                                    }
                                    
                                    if (!inherits(current_data_matrix, "sparseMatrix")) {
                                      message_to_worker_log(paste0("  Converting 'data' assay (logcounts) for '", current_tissue_name, "' to sparse matrix to save memory."))
                                      tissue_seurat@assays$RNA@data <- Matrix::Matrix(current_data_matrix, sparse = TRUE)
                                    } else {
                                      message_to_worker_log(paste0("  'data' assay (logcounts) for '", current_tissue_name, "' is already sparse. No conversion needed."))
                                    }
                                    
                                    if (ncol(tissue_seurat) == 0) {
                                      stop("Loaded Seurat object is empty")
                                    }
                                    
                                    message_to_worker_log(paste0("Performing pre-MAST checks for ", current_tissue_name, "."))
                                    if (ncol(tissue_seurat) < min_cells_per_tissue) {
                                      stop(paste0("Insufficient cells (", ncol(tissue_seurat), " < ", min_cells_per_tissue, ")"))
                                    }
                                    
                                    num_subjects <- length(levels(tissue_seurat@meta.data$donor_id))
                                    num_sex_groups <- length(levels(tissue_seurat@meta.data$sex))
                                    num_distinct_ages <- length(unique(tissue_seurat@meta.data$age))
                                    
                                    if (num_subjects < 2 || num_sex_groups < 2 || num_distinct_ages < 3) {
                                      stop(paste0("Insufficient variation for regression (Subjects: ", num_subjects, ", Sex groups: ", num_sex_groups, ", Distinct ages: ", num_distinct_ages, ")"))
                                    }
                                    
                                    sce_tissue <- as.SingleCellExperiment(tissue_seurat, assay = "RNA")
                                    rm(tissue_seurat); gc(verbose = FALSE)
                                    
                                    colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                    colData(sce_tissue)$sex <- droplevels(colData(sce_tissue)$sex)
                                    
                                    subject_counts <- table(colData(sce_tissue)$donor_id)
                                    subjects_to_keep <- names(subject_counts[subject_counts > 1])
                                    
                                    if (length(subjects_to_keep) < 2) {
                                      stop("Insufficient subjects with more than one cell after filtering")
                                    }
                                    
                                    sce_tissue <- sce_tissue[, colData(sce_tissue)$donor_id %in% subjects_to_keep]
                                    colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
                                    
                                    if (ncol(sce_tissue) < min_cells_per_tissue) {
                                      stop(paste0("Insufficient cells (", ncol(sce_tissue), " < ", min_cells_per_tissue, ") after subject filtering"))
                                    }
                                    
                                    expressed_genes <- rowSums(assay(sce_tissue, "logcounts") > 0) / ncol(sce_tissue) > min_expressed_gene_threshold
                                    if (sum(expressed_genes) < min_genes_after_filter) {
                                      stop(paste0("Insufficient highly expressed genes (", sum(expressed_genes), " < ", min_genes_after_filter, ")"))
                                    }
                                    sce_tissue_filtered <- sce_tissue[expressed_genes, ]
                                    rm(sce_tissue); gc(verbose = FALSE)
                                    
                                    final_cell_count <- ncol(sce_tissue_filtered)
                                    final_gene_count <- nrow(sce_tissue_filtered)
                                    
                                    rowData(sce_tissue_filtered)$primerid <- rownames(sce_tissue_filtered)
                                    colData(sce_tissue_filtered)$wellKey <- colnames(sce_tissue_filtered)
                                    
                                    message_to_worker_log(paste0("Running MAST for '", current_tissue_name, "' with ", nrow(sce_tissue_filtered), " genes and ", ncol(sce_tissue_filtered), " cells."))
                                    
                                    sca_mast <- SceToSingleCellAssay(sce_tissue_filtered, class = "SingleCellAssay")
                                    rm(sce_tissue_filtered); gc(verbose = FALSE)
                                    
                                    if (safe_tissue_name %in% c('ovary', 'prostate_gland', 'testis')) {
                                      zlm_obj <- zlm(~ age + donor_id, sca = sca_mast, method = 'glm', ebayes = TRUE, parallel = TRUE, exprs_value = 'logcounts') 
                                    } else {
                                      zlm_obj <- zlm(~ age + sex + donor_id, sca = sca_mast, method = 'glm', ebayes = TRUE, parallel = TRUE, exprs_value = 'logcounts') 
                                    }
                                    rm(sca_mast); gc(verbose = FALSE)
                                    
                                    summary_age_results <- summary(zlm_obj, doLRT = "age")
                                    results_table_mast_raw <- summary_age_results$datatable
                                    rm(zlm_obj, summary_age_results); gc(verbose = FALSE)
                                    
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
                                    rm(results_table_mast_raw); gc(verbose = FALSE)
                                    
                                    if (is.null(results_table_mast) || nrow(results_table_mast) == 0) {
                                      stop("No differential expression results found for 'age'")
                                    }
                                    
                                    results_table_omic <- results_table_mast %>%
                                      dplyr::mutate(
                                        probe_id = PrimerID, feature_name = PrimerID, score = `logFC`,
                                        p_value = `Pvalue`, adj_p = `FDR`
                                      ) %>%
                                      dplyr::select(probe_id, feature_name, score, p_value, adj_p) %>%
                                      dplyr::mutate(
                                        group_label = as.factor(ifelse(score > 0, "Increased_with_Age", "Decreased_with_Age"))
                                      )
                                    rm(results_table_mast); gc(verbose = FALSE)
                                    
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
                                      stop(paste0("No significant genes found with current cutoffs (adj_p <= ", adj_p_cutoff, ", |logFC| >= ", score_cutoff, "). No OmicSignature object created."))
                                    } else {
                                      omic_sig_tissue <- OmicSignature$new(
                                        metadata = metadata_tissue_sig,
                                        signature = sig_genes,
                                        difexp = results_table_omic 
                                      )
                                      output_file_path <- file.path(omic_signature_output_path, paste0("aging_signature_", safe_tissue_name, "_oSig.rds"))
                                      saveRDS(omic_sig_tissue, file = output_file_path)
                                      message_to_worker_log(paste0("  Successfully created and added aging signature for ", current_tissue_name, ". (", nrow(sig_genes), " significant genes)"))
                                      
                                      status <- "Success"
                                      reason <- "OmicSignature created successfully"
                                      significant_genes_count <- nrow(sig_genes)
                                      omic_sig_object <- omic_sig_tissue
                                    }
                                    
                                  }, error = function(e) {
                                    # This error handler captures any specific "stop" message or critical error
                                    message_to_worker_log(paste0("  Error during processing for tissue '", current_tissue_name, "': ", e$message))
                                    status <- "Skipped" # Default to Skipped for known issues
                                    reason <- e$message # The reason is the error message itself
                                    # Ensure all 'final' counts are 0/N/A for skipped/error state
                                    final_cell_count <- 0
                                    final_gene_count <- 0
                                    significant_genes_count <- 0
                                    output_file_path <- "N/A"
                                    omic_sig_object <- NULL
                                  })
                                  
                                  # Always return a complete list of results, regardless of success or failure
                                  return(list(
                                    omic_sig_object = omic_sig_object,
                                    status = status,
                                    reason = reason,
                                    initial_cell_count = initial_cell_count,
                                    final_cell_count = final_cell_count,
                                    initial_gene_count = initial_gene_count,
                                    final_gene_count = final_gene_count,
                                    significant_genes_count = significant_genes_count,
                                    output_file_path = output_file_path
                                  ))
                                } # End process_tissue_inner function
                                
                                # Call the inner processing function and capture all its returned values
                                worker_output_list <- process_tissue_inner(current_tissue_name, safe_tissue_name, file_path)
                                
                                # Now, format these into the structured list for the foreach combiner
                                final_structured_result <- build_structured_result_for_foreach(
                                  current_tissue_name = current_tissue_name,
                                  status_val = worker_output_list$status,
                                  reason_val = worker_output_list$reason,
                                  initial_cells = worker_output_list$initial_cell_count,
                                  final_cells_before_mast = worker_output_list$final_cell_count,
                                  initial_genes = worker_output_list$initial_gene_count,
                                  final_genes_before_mast = worker_output_list$final_gene_count,
                                  significant_genes = worker_output_list$significant_genes_count,
                                  omic_signature_file = worker_output_list$output_file_path,
                                  omic_sig_object = worker_output_list$omic_sig_object
                                )
                                
                                # Explicitly return this structured result for the foreach combiner
                                return(final_structured_result)
                              } # End foreach loop


# --- Post-Processing of Results ---

# Stop the parallel cluster for the external loop to release resources.
if (exists("cl") && inherits(cl, "cluster")) {
  stopCluster(cl)
  message("\nStopped parallel cluster for external loop.")
}


# 1. Extract all elements named "omic_sig"
all_omic_sigs_list_of_named_lists <- all_tissue_results[names(all_tissue_results) == "omic_sig"]

# 2. Flatten this list of named lists into a single named list of OmicSignature objects (or NULLs)
all_tissue_omicsigs_for_collection <- do.call(c, all_omic_sigs_list_of_named_lists)

# 3. Filter out any NULL entries from the combined OmicSignature list
all_tissue_omicsigs_filtered <- Filter(Negate(is.null), all_tissue_omicsigs_for_collection)

# 4. Extract all elements named "summary"
all_summaries_dfs_list <- all_tissue_results[names(all_tissue_results) == "summary"]

# 5. Combine all tissue processing summaries into a single data frame
if (length(all_summaries_dfs_list) > 0) {
  tissue_processing_summary_df <- do.call(rbind, all_summaries_dfs_list)
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

