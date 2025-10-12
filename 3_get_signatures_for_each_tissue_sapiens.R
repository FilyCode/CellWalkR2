
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



# --- PARALLELISM CONFIGURATION: Setting up nested parallelization ---
# 1. Outer loop: Multiple tissues processed concurrently.
# 2. Inner loop: MAST (zlm) uses multiple cores for each tissue.

sge_total_slots <- as.numeric(Sys.getenv("NSLOTS", unset = 1)) # Get total CPU slots from SGE job

# Number of tissue analyses to run concurrently (outer loop parallelization).
# This divides the total SGE slots into independent parallel tasks.
n_concurrent_tissues <- 4 

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
                              .combine = 'list', # <<< CHANGED: Use 'list' to keep worker results separate
                              .init = list(),
                              .verbose = TRUE) %dopar% { 
                                
      # Set mc.cores for MAST's internal parallelism for THIS specific worker/tissue task.
      # This must be done inside the foreach loop for each worker's session to be effective.
      options(mc.cores = mast_cores_per_tissue) 
      
      current_tissue_name <- gsub("_", " ", gsub("TabulaSapiens_|_organ|\\.rds$", "", basename(file_path)))
      
      # Store messages specific to this tissue, to be returned and printed by the main process
      # This captures all messages, warnings, and standard output.
      output_messages <- c() 
      
      # Define a custom message handler to capture messages
      capture_message <- function(m) {
        output_messages <<- c(output_messages, conditionMessage(m))
        invokeRestart("muffleMessage") # Prevents default message printing
      }
      
      # Initialize placeholders for the worker's return values
      omic_sig_object <- NULL 
      
      # Use tryCatch for robust error handling: if an error occurs for one tissue,
      # it won't stop the entire parallel job. Messages inside are handled by withCallingHandlers.
      tryCatch({
        withCallingHandlers(
          { # All code for a single tissue analysis goes here
            message(paste0("\n--- Analyzing tissue from file: ", basename(file_path), " (", current_tissue_name, ") ---"))
            
            # Load the tissue-specific Seurat object
            tissue_seurat <- readRDS(file_path)
            
            # --- Memory Optimization: Ensure Seurat 'data' assay is a sparse matrix ---
            # This prevents high memory usage from dense matrix conversions during processing.
            if ("RNA" %in% names(tissue_seurat@assays)) {
              # Safely get the 'data' layer (log-normalized counts) from the Seurat object
              current_data_matrix <- tryCatch(
                expr = Seurat::GetAssayData(tissue_seurat, layer = "data", assay = "RNA"), 
                error = function(e) {
                  message(paste0("  Warning: Could not access 'data' layer from 'RNA' assay for '", current_tissue_name, "'. Error: ", e$message))
                  return(NULL) # Return NULL to indicate failure to retrieve data
                }
              )
              
              # If data is found and not empty, check if it's dense and convert to sparse if needed.
              if (!is.null(current_data_matrix) && prod(dim(current_data_matrix)) > 0) {
                if (!inherits(current_data_matrix, "sparseMatrix")) {
                  message(paste0("  Converting 'data' assay (logcounts) for '", current_tissue_name, "' to sparse matrix to save memory."))
                  tissue_seurat@assays$RNA@data <- Matrix::Matrix(current_data_matrix, sparse = TRUE)
                } else {
                  message(paste0("  'data' assay (logcounts) for '", current_tissue_name, "' is already sparse. No conversion needed."))
                }
              } else {
                message(paste0("  Warning: 'data' layer in 'RNA' assay is missing or empty for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure data is normalized before analysis."))
                # Skip further analysis for this tissue if data is truly missing/empty
                return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
              }
            } else {
              message(paste0("  Warning: 'RNA' assay not found in Seurat object for '", current_tissue_name, "'. Skipping sparse conversion check. Please ensure 'RNA' assay exists."))
              # Skip further analysis for this tissue
              return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
            }
            
            # Skip if loaded object is empty
            if (ncol(tissue_seurat) == 0) {
              message(paste0("  Skipping '", current_tissue_name, "': loaded object is empty."))
              return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
            }
            
            # --- Pre-MAST Data Checks and Filtering ---
            # Ensure sufficient cells, subjects, and genes for robust MAST analysis.
            if (ncol(tissue_seurat) < min_cells_per_tissue) {
              message(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(tissue_seurat), " < ", min_cells_per_tissue, ")."))
              return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
            }
            
            num_subjects <- length(levels(tissue_seurat@meta.data$donor_id))
            num_sex_groups <- length(levels(tissue_seurat@meta.data$sex))
            num_distinct_ages <- length(unique(tissue_seurat@meta.data$age))
            
            if (num_subjects < 2 || num_sex_groups < 2 || num_distinct_ages < 2) {
              message(paste0("  Skipping '", current_tissue_name, "' due to insufficient variation for regression (Subjects: ", num_subjects, ", Sex groups: ", num_sex_groups, ", Distinct ages: ", num_distinct_ages, ")."))
              return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
            }
            
            # Convert Seurat object to SingleCellExperiment (SCE) for MAST compatibility.
            sce_tissue <- as.SingleCellExperiment(tissue_seurat, assay = "RNA")
            
            # Drop unused factor levels to prevent issues in downstream models.
            colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
            colData(sce_tissue)$sex <- droplevels(colData(sce_tissue)$sex)
            
            # Filter out subjects with only one cell for MAST stability.
            subject_counts <- table(colData(sce_tissue)$donor_id)
            subjects_to_keep <- names(subject_counts[subject_counts > 1])
            
            if (length(subjects_to_keep) < 2) {
              message(paste0("  Skipping '", current_tissue_name, "' due to insufficient subjects with more than one cell (after filtering)."))
              return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
            }
            
            sce_tissue <- sce_tissue[, colData(sce_tissue)$donor_id %in% subjects_to_keep]
            colData(sce_tissue)$donor_id <- droplevels(colData(sce_tissue)$donor_id)
            
            if (ncol(sce_tissue) < min_cells_per_tissue) {
              message(paste0("  Skipping '", current_tissue_name, "' due to insufficient cells (", ncol(sce_tissue), " < ", min_cells_per_tissue, ") after subject filtering."))
              return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
            }
            
            # Filter genes: keep only those expressed in a minimum percentage of cells.
            expressed_genes <- rowSums(assay(sce_tissue, "logcounts") > 0) / ncol(sce_tissue) > min_expressed_gene_threshold
            if (sum(expressed_genes) < min_genes_after_filter) {
              message(paste0("  Skipping '", current_tissue_name, "' due to insufficient highly expressed genes (", sum(expressed_genes), " < ", min_genes_after_filter, ")."))
              return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
            }
            sce_tissue_filtered <- sce_tissue[expressed_genes, ]
            
            # Explicitly define primerid (gene ID) and wellKey (cell ID) for MAST.
            rowData(sce_tissue_filtered)$primerid <- rownames(sce_tissue_filtered)
            colData(sce_tissue_filtered)$wellKey <- colnames(sce_tissue_filtered)
            
            # --- DEBUGGING CHECKS: Verify data consistency before MAST processing ---
            # These checks help catch dimension or identifier mismatches that can lead to errors.
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
            
            # Confirm sparsity before conversion to SingleCellAssay (MAST's required format).
            if (inherits(assay(sce_tissue_filtered, "logcounts"), "sparseMatrix")) {
              message(paste0("  DEBUG: 'logcounts' assay in sce_tissue_filtered is sparse before SceToSingleCellAssay."))
            } else {
              message(paste0("  DEBUG: 'logcounts' assay in sce_tissue_filtered is dense before SceToSingleCellAssay. This might trigger coercion."))
            }
            
            sca_mast <- SceToSingleCellAssay(sce_tissue_filtered, class = "SingleCellAssay")
            
            # Confirm sparsity after conversion to SingleCellAssay.
            if (inherits(assay(sca_mast, "logcounts"), "sparseMatrix")) {
              message(paste0("  DEBUG: 'logcounts' assay in sca_mast is sparse after SceToSingleCellAssay."))
            } else {
              message(paste0("  DEBUG: 'logcounts' assay in sca_mast is dense after SceToSingleCellAssay. This is where dense conversion for MAST occurs."))
            }
            
            # Fit the ZLM model: gene ~ age + sex + donor_id.
            # 'parallel = TRUE' tells MAST to use the cores set by options(mc.cores) for this worker.
            zlm_obj <- zlm(~ age + sex + donor_id, sca = sca_mast, method = 'glm', ebayes = TRUE, parallel = TRUE) 
            
            # Get summary results for the 'age' coefficient using a Likelihood Ratio Test (doLRT).
            summary_age_results <- summary(zlm_obj, doLRT = "age")
            results_table_mast_raw <- summary_age_results$datatable
            
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
            
            # Skip if no differential expression results found for 'age'.
            if (is.null(results_table_mast) || nrow(results_table_mast) == 0) {
              message(paste0("  No differential expression results found for 'age' in tissue: ", current_tissue_name))
              return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
            } else {
              # Prepare results for OmicSignature object (difexp data frame).
              results_table_omic <- results_table_mast %>%
                dplyr::mutate(
                  probe_id = PrimerID, feature_name = PrimerID, score = `logFC`,
                  p_value = `Pvalue`, adj_p = `FDR`
                ) %>%
                dplyr::select(probe_id, feature_name, score, p_value, adj_p) %>%
                dplyr::mutate(
                  # Define group_label for bi-directional signature (required by OmicSignature).
                  group_label = as.factor(ifelse(score > 0, "Increased_with_Age", "Decreased_with_Age"))
                )
              
              # --- Automated BRENDA ontology lookup for sample_type ---
              # Tries to find a suitable biological sample type term for the tissue.
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
                return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages))
              } else {
                # Create the OmicSignature object for the current tissue.
                omic_sig_object <- OmicSignature$new( # Assign to omic_sig_object, not directly return
                  metadata = metadata_tissue_sig,
                  signature = sig_genes,
                  difexp = results_table_omic # Store the full differential expression results
                )
                
                # Save individual OmicSignature object (for easier access) within each worker.
                safe_tissue_name <- gsub("[^[:alnum:]_]", "_", current_tissue_name) 
                saveRDS(omic_sig_object, file = file.path(omic_signature_output_path, paste0("aging_signature_", safe_tissue_name, "_oSig.rds")))
                message(paste0("  Successfully created and added aging signature for ", current_tissue_name, ". (", nrow(sig_genes), " significant genes)"))
                
                # Final return for a successful worker task
                return(list(omicSig = omic_sig_object, tissueName = current_tissue_name, messages = output_messages))
              }
            } # End if/else for results_table_mast
          }, # End withCallingHandlers code block
          message = capture_message, warning = capture_message # Capture messages and warnings
        ) # End withCallingHandlers
      }, error = function(e) {
        # Catches any error during tissue analysis and logs it without stopping the script.
        # Ensures error messages are also captured and returned.
        output_messages <<- c(output_messages, paste0("  ERROR: ", e$message))
        message(paste0("  Error during MAST or OmicSignature creation for tissue '", current_tissue_name, "': ", e$message))
        return(list(omicSig = NULL, tissueName = current_tissue_name, messages = output_messages)) for errors
      }) 
    } # End foreach loop




# --- Post-Processing of Results ---

# Separate the OmicSignature objects from the captured messages.
all_tissue_omicsigs <- list()
all_captured_messages <- list()

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
}

# Print all captured messages from workers to the main log file.
if (length(all_captured_messages) > 0) {
  message("\n--- Captured messages from parallel workers ---")
  for (msg in all_captured_messages) {
    message(msg)
  }
  message("--- End captured messages ---\n")
}


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

