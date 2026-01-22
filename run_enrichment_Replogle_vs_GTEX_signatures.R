suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(OmicSignature)
  library(ggplot2)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(cowplot)
  library(BiocParallel)
  library(purrr)
  library(data.table)
})

# --- Configuration ---
# File paths
work_dir <- "/restricted/projectnb/agedisease/projects/challenge2025/"
perturbation_collection_file <- paste0(work_dir, "results/perturbational_omic_sigs/replogle_2022/Replogle_Perturb_Combined_OmicSignatureCollection.rds")
gtex_collection_file <- paste0(work_dir, "data/GTEX/GTEX_aging_omic_col_stat_v114.rds")
output_dir <- paste0(work_dir, "results/fgsea_combined_analysis/Replogle_GTEX")
gene_map_file <- paste0(work_dir, "data/Homo_sapiens_GRCh38_114_genemap.rds")
essential_genes_file <- paste0(work_dir, "data/Gene_Dependency_Profile_Summary.csv")

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Filtering parameters for gene sets
logFC_filter_val <- 0.25
adj_pval_filter_val <- 0.05
geneset_top_n <- 500
top_n_genesets_to_plot <- 20 # For existing 4-panel plots
top_n_combined_plots <- 20 # For new combined NES dot plots

# fgsea parameters
fgsea_min_size <- 15 
fgsea_max_size <- Inf

# Parallelization strategy
total_num_cores <- 28


# Load essential genes from DepMap
message("--- Loading Essential Genes from DepMap ---")
full_essential_genes_df <- read.csv(essential_genes_file)
nr_all_depmap_genes = nrow(unique(full_essential_genes_df))
essential_genes_df <- read.csv(essential_genes_file) %>%
  dplyr::filter(Dataset == "DependencyEnum.Chronos_Combined", Common.Essential == 'True')
essential_gene_list <- unique(essential_genes_df$Gene)
nr_ess_depmap_genes = length(essential_gene_list)
message(paste0("Loaded ", length(essential_gene_list), " common essential genes from DepMap for highlighting."))


message("--- Starting GSEA Combined Analysis ---")
message(paste0("Output directory: ", output_dir))
message(paste0("Gene set filtering: |logFC| > ", logFC_filter_val, ", adj.pval < ", adj_pval_filter_val))
message(paste0("fgsea parameters: minSize = ", fgsea_min_size, ", maxSize = ", fgsea_max_size, 
               ", total_num_cores = ", total_num_cores, 
               " (each fgsea task will run on 1 core, with up to ", total_num_cores, " tasks concurrently)"))
message(paste0("Top N results for plotting (traditional): ", top_n_genesets_to_plot))
message(paste0("Top N results for plotting (combined NES): ", top_n_combined_plots))

# --- Helper Functions ---

#' Simplifies gene names, taking the shortest alphabetical entry from a ' /// ' delimited string.
#'
#' @param x A string representing one or more gene names, possibly ' /// ' delimited.
#' @return A simplified, uppercase gene name.
simplify_entry <- function(x) {
  if (is.null(x) || is.na(x) || x == "" || x == "<NA>") { # Added "<NA>" check
    return(NA_character_)
  }
  # Split on ' /// '
  parts <- strsplit(x, " /// ")[[1]]
  # Filter out empty strings if any
  parts <- parts[parts != ""]
  if (length(parts) == 0) {
    return(NA_character_)
  }
  # Get lengths
  len <- nchar(parts)
  shortest <- parts[len == min(len)]
  # If tie, sort alphabetically
  chosen <- sort(shortest)[1]
  toupper(chosen)
}


#' Extracts ranked gene lists from an OmicSignatureCollection.
#' Ranks genes by 'score' without p-value/logFC filtering.
#' Prioritizes 'difexp' table, falls back to 'signature' if 'difexp' is missing.
#' Maps probe_id to gene_name using a provided *simplified* gene map.
#'
#' @param omic_collection An OmicSignatureCollection object.
#' @param collection_name A string identifying the source collection (for warnings/messages).
#' @param simplified_gene_map A named character vector for mapping probe_id to gene_name (already simplified).
#' @return A named list of numeric vectors, where names are gene_name and values are scores.
extract_ranked_lists <- function(omic_collection, collection_name, simplified_gene_map) { 
  ranked_lists <- list()
  if (is.null(omic_collection) || is.null(omic_collection$OmicSigList)) {
    stop(paste("OmicSignatureCollection is NULL or empty for", collection_name))
  }
  
  message(paste0("Preparing ranked lists from '", collection_name, "' collection..."))
  for (sig_name in names(omic_collection$OmicSigList)) {
    sig_obj <- omic_collection$OmicSigList[[sig_name]]
    
    df_for_ranks <- NULL
    if (!is.null(sig_obj$difexp) && all(c("probe_id", "score") %in% names(sig_obj$difexp))) {
      df_for_ranks <- sig_obj$difexp
    } else if (!is.null(sig_obj$signature) && all(c("probe_id", "score") %in% names(sig_obj$signature))) {
      warning(paste("Using 'signature' for ranked list instead of 'difexp' for", collection_name, sig_name, "as 'difexp' is missing or incomplete."))
      df_for_ranks <- sig_obj$signature
    } else {
      warning(paste("Skipping", collection_name, sig_name, ": 'difexp' or 'signature' missing or does not contain required columns (probe_id, score)."))
      next
    }
    
    # Use the pre-simplified gene_map for direct lookup if we have ENSG* names and not gene symbols, otherwise just use probe_id
    # preserve original order
    df_for_ranks$orig_row <- seq_len(nrow(df_for_ranks))
    
    # detect ENSG probe_ids
    is_ensg <- grepl("^ENSG", df_for_ranks$probe_id)
    
    # vectorized lookup for ENSG rows (fast named-vector indexing)
    gene_name <- rep(NA_character_, nrow(df_for_ranks))
    gene_name[is_ensg] <- simplified_gene_map[df_for_ranks$probe_id[is_ensg]]
    
    # fallback: use probe_id when mapping failed or for non-ENSG rows
    na_idx <- is.na(gene_name)
    gene_name[na_idx] <- df_for_ranks$probe_id[na_idx]
    
    # assign gene_name
    df_for_ranks$gene_name <- gene_name
    
    # deduplicate only among successfully mapped ENSG rows (keep first occurrence)
    mapped_idx <- is_ensg & df_for_ranks$gene_name != df_for_ranks$probe_id
    keep <- rep(TRUE, nrow(df_for_ranks))
    dup_pos <- which(mapped_idx)[duplicated(df_for_ranks$gene_name[mapped_idx])]
    keep[dup_pos] <- FALSE
    
    # filter and restore original order; drop helper column
    df_for_ranks <- df_for_ranks[keep, , drop = FALSE]
    df_for_ranks <- df_for_ranks[order(df_for_ranks$orig_row), ]
    df_for_ranks$orig_row <- NULL
    
    current_ranks <- df_for_ranks %>%
      dplyr::select(gene_name, score) %>% 
      drop_na(score, gene_name) %>% # Remove NAs in score and gene_name 
      distinct(gene_name, .keep_all = TRUE) %>% # Keep one entry per gene_name if duplicates exist 
      arrange(desc(score)) %>%
      tibble::deframe() # Convert to named numeric vector
    
    if (length(current_ranks) > 0) {
      ranked_lists[[sig_name]] <- current_ranks
    } else {
      warning(paste("Skipping", collection_name, sig_name, ": ranked list is empty after processing."))
    }
  }
  
  if (length(ranked_lists) == 0) {
    stop(paste("No valid ranked lists could be prepared from", collection_name, "collection."))
  }
  message(paste0("Prepared ", length(ranked_lists), " ranked lists from '", collection_name, "' collection."))
  return(ranked_lists)
}



#' Extracts gene sets (upregulated and downregulated) from an OmicSignatureCollection.
#' Filters genes based on logFC and adjusted p-value thresholds.
#' Prioritizes 'signature' table, falls back to 'difexp' if 'signature' is missing.
#' Maps probe_id to gene_name using a provided *simplified* gene map.
#'
#' @param omic_collection An OmicSignatureCollection object.
#' @param logFC_thresh Numeric, absolute logFC threshold for filtering.
#' @param pval_thresh Numeric, adjusted p-value threshold for filtering.
#' @param geneset_top_n Numeric or NULL. If a number, takes top N genes by absolute logFC after other filters.
#' @param collection_name A string identifying the source collection (for warnings/messages).
#' @param simplified_gene_map A named character vector for mapping probe_id to gene_name (already simplified).
#' @return A list containing two named lists: 'up_gene_sets' and 'dn_gene_sets'.
extract_gene_sets_up_dn <- function(omic_collection, logFC_thresh, pval_thresh, geneset_top_n, collection_name, simplified_gene_map) { 
  up_gene_sets <- list()
  dn_gene_sets <- list()
  
  if (is.null(omic_collection) || is.null(omic_collection$OmicSigList)) {
    stop(paste("OmicSignatureCollection is NULL or empty for", collection_name))
  }
  
  message(paste0("Preparing UP and DN gene sets from '", collection_name, "' collection..."))
  
  # Define expected column names. logFC is now conditionally used.
  required_id_col <- "probe_id"
  required_score_col <- "score"
  possible_logfc_cols <- c("logFC") # Standard logFC column name
  possible_group_label_cols <- c("group_label") # Standard group_label column name
  possible_pval_cols <- c("adj.pval", "adj_p", "adj_pval", "adj_p_val", "pval_adj") 
  
  # Helper function to get a valid dataframe and standardize p-value column
  # A dataframe is considered valid if it has probe_id, score, and any of the possible p-value columns.
  get_valid_df <- function(candidate_df, source_str) {
    if (!is.null(candidate_df)) {
      current_cols <- names(candidate_df)
      has_basic_cols <- all(c(required_id_col, required_score_col) %in% current_cols)
      has_pval_col <- any(possible_pval_cols %in% current_cols)
      
      if (has_basic_cols && has_pval_col) {
        # Standardize p-value column name to 'adj.pval'
        pval_col_name <- intersect(possible_pval_cols, current_cols)[1]
        if (pval_col_name != "adj.pval") {
          candidate_df <- candidate_df %>% dplyr::rename(adj.pval = !!sym(pval_col_name))
        }
        return(candidate_df)
      }
    }
    return(NULL)
  }
  
  for (sig_name in names(omic_collection$OmicSigList)) {
    sig_obj <- omic_collection$OmicSigList[[sig_name]]
    
    df_to_process <- NULL
    used_source_name <- "none" 
    
    
    # Attempt to get dataframe, preferring 'difexp' over 'signature'
    df_to_process <- get_valid_df(sig_obj$difexp, "difexp")
    if (!is.null(df_to_process)) {
      used_source_name <- "difexp"
    } else {
      df_to_process <- get_valid_df(sig_obj$signature, "signature")
      if (!is.null(df_to_process)) {
        used_source_name <- "signature"
        # Only warn about using 'signature' if 'difexp' was actually present but invalid
        if (!is.null(sig_obj$difexp)) { 
          warning(paste("Using 'signature' for gene set filtering for", collection_name, sig_name, 
                        "because 'difexp' was present but lacked required columns (probe_id, score, or adjusted p-value)."))
        }
      }
    }
    
    # If no suitable dataframe was found after checking both, skip this signature
    if (is.null(df_to_process)) {
      warning(paste("Skipping", collection_name, sig_name, 
                    ": Neither 'difexp' nor 'signature' contains the minimum required columns (probe_id, score, and any of", paste(possible_pval_cols, collapse="/"), ")."))
      next
    }
    
    # Check for logFC availability
    has_logFC <- any(possible_logfc_cols %in% names(df_to_process))
    if (!has_logFC) {
      warning(paste("Warning: 'logFC' column missing for", collection_name, sig_name, 
                    "from", used_source_name, ". Will use 'score' for determining UP/DN gene set direction.",
                    "logFC_thresh will not be applied for this signature."))
    }
    
    # Determine if logFC is available and if group_label is a valid fallback for directionality
    has_logFC <- any(possible_logfc_cols %in% names(df_to_process)) && 
      !all(is.na(df_to_process[[possible_logfc_cols[1]]])) # logFC column exists and is not all NA
    
    use_group_label_for_direction <- FALSE
    group_label_col_name <- NULL
    
    if (!has_logFC) { # If logFC is not available or all NA
      if (any(possible_group_label_cols %in% names(df_to_process))) {
        group_label_col_name <- intersect(possible_group_label_cols, names(df_to_process))[1]
        # Check for presence of *any* non-NA group label
        if (!all(is.na(df_to_process[[group_label_col_name]]))) { 
          use_group_label_for_direction <- TRUE
          warning(paste("Warning: 'logFC' column missing or all NA for", collection_name, sig_name, 
                        "from", used_source_name, ". Using '", group_label_col_name, "' for determining UP/DN gene set direction.",
                        "logFC_thresh will not be applied for this signature."))
        } else {
          warning(paste("Skipping", collection_name, sig_name, ": 'logFC' column missing or all NA, and '", group_label_col_name, "' is either missing required values ('Older', 'Younger') or is all NA. Directionality for gene sets cannot be determined."))
          next
        }
      } else {
        warning(paste("Skipping", collection_name, sig_name, ": Missing both 'logFC' and 'group_label' columns. Directionality for gene sets cannot be determined."))
        next
      }
    }
    
    # Filter for NAs in critical columns and ensure unique gene IDs
    df_filtered <- df_to_process %>% 
      dplyr::filter(!is.na(!!sym(required_score_col)), !is.na(adj.pval)) %>% 
      dplyr::distinct(!!sym(required_id_col), .keep_all = TRUE) 
    
    # Filter NAs for directionality columns if they are being used
    if (has_logFC) {
      df_filtered <- df_filtered %>% dplyr::filter(!is.na(!!sym(possible_logfc_cols[1])))
    } else if (use_group_label_for_direction) {
      df_filtered <- df_filtered %>% dplyr::filter(!is.na(!!sym(group_label_col_name)))
    }
    
    # Apply p-value filtering
    significant_genes <- df_filtered %>%
      dplyr::filter(adj.pval < pval_thresh)
    
    # Apply logFC threshold ONLY IF logFC is present and being used for directionality
    if (has_logFC) { 
      significant_genes <- significant_genes %>% dplyr::filter(abs(!!sym(possible_logfc_cols[1])) > logFC_thresh)
    }
    
    if (nrow(significant_genes) == 0) {
      message(paste("No significant genes found for", collection_name, sig_name, "after filtering. Skipping gene set creation."))
      next
    }
    
    # Apply top_n filter if specified
    if (!is.null(geneset_top_n) && is.numeric(geneset_top_n) && geneset_top_n > 0) {
      significant_genes <- significant_genes %>% 
        dplyr::arrange(desc(abs(!!sym(required_score_col)))) %>% # Rank by absolute score.
        dplyr::slice_head(n = geneset_top_n) # Take top N genes
    } 
    
    # Use the pre-simplified gene_map for direct lookup if we have ENSG* names and not gene symbols, otherwise just use probe_id
    # preserve original order
    significant_genes$orig_row <- seq_len(nrow(significant_genes))
    
    # detect ENSG probe_ids
    is_ensg <- grepl("^ENSG", significant_genes$probe_id)
    
    # vectorized lookup for ENSG rows (fast named-vector indexing)
    gene_name <- rep(NA_character_, nrow(significant_genes))
    gene_name[is_ensg] <- simplified_gene_map[significant_genes$probe_id[is_ensg]]
    
    # fallback: use probe_id when mapping failed or for non-ENSG rows
    na_idx <- is.na(gene_name)
    gene_name[na_idx] <- significant_genes$probe_id[na_idx]
    
    # assign gene_name
    significant_genes$gene_name <- gene_name
    
    # deduplicate only among successfully mapped ENSG rows (keep first occurrence)
    mapped_idx <- is_ensg & significant_genes$gene_name != significant_genes$probe_id
    keep <- rep(TRUE, nrow(significant_genes))
    dup_pos <- which(mapped_idx)[duplicated(significant_genes$gene_name[mapped_idx])]
    keep[dup_pos] <- FALSE
    
    # filter and restore original order; drop helper column
    significant_genes <- significant_genes[keep, , drop = FALSE]
    significant_genes <- significant_genes[order(significant_genes$orig_row), ]
    significant_genes$orig_row <- NULL
    
    if (nrow(significant_genes) == 0) {
      message(paste("No significant genes with valid gene names found for", collection_name, sig_name, "after mapping. Skipping gene set creation."))
      next
    }
    
    # Filter for upregulated genes and downregulated genes (dynamic directionality)
    up_genes <- character(0)
    dn_genes <- character(0)
    
    if (has_logFC) { 
      up_genes <- significant_genes %>% dplyr::filter(!!sym(possible_logfc_cols[1]) > 0) %>% dplyr::pull(gene_name) %>% unique() 
      dn_genes <- significant_genes %>% dplyr::filter(!!sym(possible_logfc_cols[1]) < 0) %>% dplyr::pull(gene_name) %>% unique()
    } else if (use_group_label_for_direction) {
      # Handle different possible group_label values (Aging vs. Perturbation)
      if (any(grepl("Increased_with_Age|Increased_by_Perturbation", unique(significant_genes[[group_label_col_name]])))) {
        up_genes <- significant_genes %>% dplyr::filter(grepl("Increased_with_Age|Increased_by_Perturbation", !!sym(group_label_col_name))) %>% dplyr::pull(gene_name) %>% unique()
        dn_genes <- significant_genes %>% dplyr::filter(grepl("Decreased_with_Age|Decreased_by_Perturbation", !!sym(group_label_col_name))) %>% dplyr::pull(gene_name) %>% unique()
      } else if (any(c("Older", "Younger") %in% unique(significant_genes[[group_label_col_name]]))) {
        up_genes <- significant_genes %>% dplyr::filter(!!sym(group_label_col_name) == "Older") %>% dplyr::pull(gene_name) %>% unique()
        dn_genes <- significant_genes %>% dplyr::filter(!!sym(group_label_col_name) == "Younger") %>% dplyr::pull(gene_name) %>% unique()
      } else {
        warning(paste("Skipping", collection_name, sig_name, ": Unrecognized group_label values for directionality. Expected 'Increased/Decreased' or 'Older/Younger'."))
        next
      }
    }
    
    # Add to gene sets list
    if (length(up_genes) > 0) {
      up_gene_sets[[sig_name]] <- up_genes
    }
    if (length(dn_genes) > 0) {
      dn_gene_sets[[sig_name]] <- dn_genes
    }
  }
  
  if (length(up_gene_sets) == 0 && length(dn_gene_sets) == 0) {
    stop(paste("No valid gene sets could be prepared from", collection_name, "collection after filtering."))
  }
  message(paste0("Prepared ", length(up_gene_sets), " UP gene sets and ", length(dn_gene_sets), " DN gene sets from '", collection_name, "' collection."))
  
  return(list(up_gene_sets = up_gene_sets, dn_gene_sets = dn_gene_sets))
}



#' Performs fgsea analysis for a set of ranked lists against up/down gene sets and combines the results.
#'
#' @param ranked_lists A named list of ranked gene score vectors.
#' @param gene_sets_up A named list of upregulated gene sets.
#' @param gene_sets_dn A named list of downregulated gene sets.
#' @param ranked_source_name A string, e.g., "Aging" or "Perturbation".
#' @param geneset_source_name A string, e.g., "Aging" or "Perturbation".
#' @param num_cores_per_fgsea_task Numeric, number of CPU cores for *each* fgsea call. 
#' @param num_concurrent_fgsea_tasks Numeric, number of fgsea calls to run concurrently. 
#' @return A combined data frame of fgsea results.
perform_fgsea_and_combine <- function(ranked_lists, gene_sets_up, gene_sets_dn, 
                                      ranked_source_name, geneset_source_name,
                                      total_num_cores, fgsea_min_size, fgsea_max_size) {
  
  message(paste0("Initiating fgsea analysis for '", ranked_source_name, "' ranked lists vs. '", geneset_source_name, "' gene sets."))
  
  # Create a list of all individual fgsea tasks to be executed
  fgsea_tasks <- list()
  
  # Add tasks for UP gene sets
  if (length(gene_sets_up) > 0) {
    for (rl_name in names(ranked_lists)) {
      fgsea_tasks[[length(fgsea_tasks) + 1]] <- list(
        pathways = gene_sets_up,
        stats = ranked_lists[[rl_name]],
        ranked_list_name = rl_name,
        geneset_source_name = geneset_source_name,
        geneset_direction = "UP",
        ranked_source_name = ranked_source_name
      )
    }
  }
  
  # Add tasks for DN gene sets
  if (length(gene_sets_dn) > 0) {
    for (rl_name in names(ranked_lists)) {
      fgsea_tasks[[length(fgsea_tasks) + 1]] <- list(
        pathways = gene_sets_dn,
        stats = ranked_lists[[rl_name]],
        ranked_list_name = rl_name,
        geneset_source_name = geneset_source_name,
        geneset_direction = "DN",
        ranked_source_name = ranked_source_name
      )
    }
  }
  
  if (length(fgsea_tasks) == 0) {
    warning(paste("No fgsea tasks could be prepared for ", ranked_source_name, " vs ", geneset_source_name, " analysis type."))
    return(NULL)
  }
  
  message(paste0("  Prepared ", length(fgsea_tasks), " individual fgsea tasks. Running with ", total_num_cores, " concurrent workers."))
  
  # Store the current default BPPARAM to restore it later to avoid side effects
  old_bpparam_registered <- BiocParallel::bpparam()
  on.exit(BiocParallel::register(old_bpparam_registered, default = TRUE), add = TRUE)
  
  # Register MulticoreParam for bplapply itself to use multiple cores for multiple tasks.
  # Change progressbar=FALSE to progressbar=TRUE here to enable the overall progress bar.
  BiocParallel::register(BiocParallel::MulticoreParam(workers = total_num_cores, progressbar = TRUE), default = TRUE)
  
  # Execute all fgsea tasks in parallel
  all_fgsea_results <- BiocParallel::bplapply(fgsea_tasks, function(task) {
    # Each fgsea call runs serially (nproc = 1) and explicitly disables its own progress bar
    res <- fgsea(pathways = task$pathways,
                 stats    = task$stats,
                 minSize  = fgsea_min_size,
                 maxSize  = fgsea_max_size,
                 nproc    = 1, 
                 BPPARAM  = BiocParallel::SerialParam(progressbar = FALSE)) # Explicitly disable fgsea's internal progressbar
    
    # Ensure 'res' is a data.table to use set()
    # This check is defensive; fgsea generally returns data.table
    if (!inherits(res, "data.table")) {
      res <- data.table(res)
    }
    
    # Add metadata to the results using data.table::set
    # This is a more explicit and robust way to add columns to a data.table
    data.table::set(res, j = "ranked_list_name", value = task$ranked_list_name)
    data.table::set(res, j = "geneset_source_name", value = task$geneset_source_name)
    data.table::set(res, j = "geneset_direction", value = task$geneset_direction)
    data.table::set(res, j = "ranked_source_name", value = task$ranked_source_name) # Ensure this column is explicitly added
    
    return(res)
  }, BPPARAM = BiocParallel::bpparam()) # Use the global bpparam for this bplapply
  
  # Filter out NULL results (if any failed) and combine
  combined_results_df <- dplyr::bind_rows(all_fgsea_results[!sapply(all_fgsea_results, is.null)])
  
  if (is.null(combined_results_df) || nrow(combined_results_df) == 0) {
    warning(paste("No fgsea results generated for ", ranked_source_name, " vs ", geneset_source_name, " analysis type."))
    return(NULL)
  }
  message(paste0("fgsea analysis complete for '", ranked_source_name, "' vs. '", geneset_source_name, "'. Total results: ", nrow(combined_results_df), " rows."))
  return(combined_results_df)
}


#' Generates and saves a 2x2 grid of ggplot dot plots for fgsea results.
#' The grid shows combinations of UP/DOWN gene sets with Positive/Negative NES.
#'
#' @param fgsea_df A data frame of combined fgsea results.
#' @param analysis_title_prefix A string for plot titles, e.g., "Age-Centered Analysis".
#' @param output_dir Path to save the plots.
#' @param top_n Numeric, number of top results to plot for each panel.
plot_fgsea_results <- function(fgsea_df, analysis_title_prefix, output_dir, top_n = 20) {
  
  if (is.null(fgsea_df) || nrow(fgsea_df) == 0) {
    message(paste("No fgsea results for", analysis_title_prefix, "to display or plot."))
    return(invisible(NULL))
  }
  
  message(paste0("Generating 4-panel dot plots for: ", analysis_title_prefix))
  
  # Ensure output directory exists (already done in main, but good for standalone)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # First, filter for all significant results
  all_significant_results <- fgsea_df %>%
    filter(padj < 0.05)
  
  if (nrow(all_significant_results) == 0) {
    message(paste("No significant interactions (padj < 0.05) found for", analysis_title_prefix, "to plot."))
    return(invisible(NULL))
  } else {
    message(paste0("Found ", nrow(all_significant_results), " significant interactions (padj < 0.05) for ", analysis_title_prefix, "."))
    
    # Check for required columns - NES is used for plotting here
    if (!("geneset_direction" %in% colnames(all_significant_results) &&
          "NES" %in% colnames(all_significant_results) &&
          "ranked_list_name" %in% colnames(all_significant_results) &&
          "pathway" %in% colnames(all_significant_results))) {
      stop("Error: fgsea_df must contain 'geneset_direction', 'NES', 'ranked_list_name', and 'pathway' columns for 4-panel plot.")
    }
    
    # Determine what type of gene sets are being used (e.g., "Perturbation Gene Set", "Aging Gene Set")
    gene_set_label <- unique(all_significant_results$geneset_source_name)[1]
    if (is.null(gene_set_label) || is.na(gene_set_label)) {
      gene_set_label <- "Gene Set" # Fallback
    } else {
      gene_set_label <- paste0(gene_set_label, " Gene Set")
    }
    
    # Helper function to create a single GSEA plot panel
    create_gsea_plot_panel <- function(data, panel_title, nes_pos = TRUE) {
      if (nrow(data) == 0) {
        return(ggplot() + geom_text(aes(x=0.5, y=0.5, label=paste0("No significant data for\n", panel_title)), size=4, color="grey50") + theme_void())
      }
      
      # Define color gradient based on NES direction
      if (nes_pos) {
        color_scale <- scale_color_gradient(low = "yellow", high = "red", name = "NES")
      } else { # Negative NES
        color_scale <- scale_color_gradient(low = "darkblue", high = "lightblue", name = "NES")
      }
      
      p <- ggplot(data, aes(x = reorder(pathway, NES), y = ranked_list_name)) +
        geom_point(aes(size = size, color = NES)) +
        scale_size_continuous(name = "Gene Set Size") +
        color_scale +
        coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(face = "bold", hjust = 0.5, size = 10), # Adjust title size for grid
              axis.title = element_text(size = 8),
              axis.text = element_text(size = 7),
              legend.text = element_text(size = 7),
              legend.title = element_text(size = 8),
              legend.position = "bottom", # Place legend at bottom for more plot space
              plot.margin = margin(5, 5, 5, 5, "pt")) + # Add a small margin around each plot
        labs(x = gene_set_label, # Dynamic x-axis label (e.g., "Perturbation Gene Set")
             y = paste0(unique(data$ranked_source_name), " Ranked List"), # Dynamic y-axis label (e.g., "Aging Ranked List")
             title = panel_title)
      return(p)
    }
    
    plot_list <- list() # This list will hold the 4 ggplot objects
    
    # --- 1. Gene Sets UP, Positive NES enrichment ---
    data_up_pos <- all_significant_results %>%
      dplyr::filter(geneset_direction == 'UP', NES > 0) %>%
      dplyr::arrange(dplyr::desc(NES)) %>%
      dplyr::slice(1:min(dplyr::n(), top_n))
    
    plot_list$up_pos <- create_gsea_plot_panel(data_up_pos, 
                                               paste0("Top ", top_n, " UP ", gene_set_label, ", Positive NES"), 
                                               nes_pos = TRUE)
    message(paste0("Prepared plot for Top ", top_n, " UP ", gene_set_label, ", Positive NES (", nrow(data_up_pos), " results)."))
    
    # --- 2. Gene Sets UP, Negative NES enrichment ---
    data_up_neg <- all_significant_results %>%
      dplyr::filter(geneset_direction == 'UP', NES < 0) %>%
      dplyr::arrange(NES) %>% # Arrange by NES ascending for most negative first
      dplyr::slice(1:min(dplyr::n(), top_n))
    
    plot_list$up_neg <- create_gsea_plot_panel(data_up_neg, 
                                               paste0("Top ", top_n, " UP ", gene_set_label, ", Negative NES"), 
                                               nes_pos = FALSE)
    message(paste0("Prepared plot for Top ", top_n, " UP ", gene_set_label, ", Negative NES (", nrow(data_up_neg), " results)."))
    
    # --- 3. Gene Sets DOWN, Positive NES enrichment ---
    data_down_pos <- all_significant_results %>%
      dplyr::filter(geneset_direction == 'DN', NES > 0) %>%
      dplyr::arrange(dplyr::desc(NES)) %>% # Arrange by NES descending for most positive first
      dplyr::slice(1:min(dplyr::n(), top_n))
    
    plot_list$down_pos <- create_gsea_plot_panel(data_down_pos, 
                                                 paste0("Top ", top_n, " DOWN ", gene_set_label, ", Positive NES"), 
                                                 nes_pos = TRUE)
    message(paste0("Prepared plot for Top ", top_n, " DOWN ", gene_set_label, ", Positive NES (", nrow(data_down_pos), " results)."))
    
    # --- 4. Gene Sets DOWN, Negative NES enrichment ---
    data_down_neg <- all_significant_results %>%
      dplyr::filter(geneset_direction == 'DN', NES < 0) %>%
      dplyr::arrange(NES) %>% # Arrange by NES ascending for most negative first
      dplyr::slice(1:min(dplyr::n(), top_n))
    
    plot_list$down_neg <- create_gsea_plot_panel(data_down_neg, 
                                                 paste0("Top ", top_n, " DOWN ", gene_set_label, ", Negative NES"), 
                                                 nes_pos = FALSE)
    message(paste0("Prepared plot for Top ", top_n, " DOWN ", gene_set_label, ", Negative NES (", nrow(data_down_neg), " results)."))
    
    # --- Combine and save the plots into a 2x2 grid ---
    if (length(plot_list) > 0) {
      combined_plot <- plot_grid(plot_list$up_pos, plot_list$up_neg, 
                                 plot_list$down_pos, plot_list$down_neg, 
                                 ncol = 2, nrow = 2, align = "hv", 
                                 labels = c("A", "B", "C", "D"), label_size = 10)
      
      # Add an overall main title for the entire grid
      final_title_text <- paste0("GSEA ", analysis_title_prefix, " (Top ", top_n, " Pathways)")
      final_title <- ggdraw() + 
        draw_label(final_title_text, 
                   fontface = 'bold', size = 16, x = 0.02, hjust = 0) + # Adjust x and hjust for left alignment
        theme(plot.margin = margin(0, 0, 0, 7, "pt"))
      
      # Combine the main title with the grid of plots
      combined_plot_with_title <- plot_grid(final_title, combined_plot, ncol = 1, rel_heights = c(0.05, 1))
      
      # Save as png
      plot_filename <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_top_", top_n, "_gsea_4_panel_enrichment.png"))
      ggsave(plot_filename, combined_plot_with_title, width = 16, height = 12) # Increased width/height for 4 plots
      
      # Save as SVG (vector file)
      plot_filename <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_top_", top_n, "_gsea_4_panel_enrichment.svg"))
      ggsave(plot_filename, combined_plot_with_title, width = 16, height = 12) # Increased width/height for 4 plots
      
      message("  4-panel GSEA plot generated and saved to: ", plot_filename)
    } else {
      message("No plots generated due to lack of significant results in any category for ", analysis_title_prefix, ".")
    }
  }
}


#' Generates and saves a clustered heatmap from fgsea results.
#'
#' @param fgsea_df A data frame of fgsea results for a single direction (UP or DN).
#' @param analysis_title_prefix A string for plot titles, e.g., "Age-Centered Analysis".
#' @param geneset_source_name The source of the gene sets (e.g., "Perturbation", "Aging").
#' @param geneset_direction The direction of the gene sets (e.g., "UP", "DN").
#' @param output_dir Path to save the plots.
plot_clustered_heatmap <- function(fgsea_df, analysis_title_prefix, geneset_source_name, geneset_direction, output_dir) {
  
  if (is.null(fgsea_df) || nrow(fgsea_df) == 0) {
    message(paste("No fgsea results for heatmap in", analysis_title_prefix, "with", geneset_source_name, geneset_direction, "gene sets."))
    return(invisible(NULL))
  }
  
  message(paste0("Generating clustered heatmap for: ", analysis_title_prefix, " with ", geneset_source_name, " ", geneset_direction, " gene sets."))
  
  # Filter for significant results (padj < 0.05)
  significant_results <- fgsea_df %>%
    dplyr::filter(padj < 0.05)
  
  if (nrow(significant_results) == 0) {
    message(paste("  No significant interactions (padj < 0.05) found for heatmap in", analysis_title_prefix, "with", geneset_source_name, geneset_direction, "gene sets."))
    return(invisible(NULL))
  }
  
  # Ensure there's enough data for a meaningful heatmap (at least 2 pathways and 2 ranked lists)
  num_pathways <- length(unique(significant_results$pathway))
  num_ranked_lists <- length(unique(significant_results$ranked_list_name))
  
  if (num_pathways < 2 || num_ranked_lists < 2) {
    message(paste("  Not enough unique pathways (", num_pathways, ") or ranked lists (", num_ranked_lists, ") for a meaningful heatmap in", analysis_title_prefix, "with", geneset_source_name, geneset_direction, "gene sets. Skipping heatmap."))
    return(invisible(NULL))
  }
  
  # Reshape data for heatmap: pathways as rows, ranked lists as columns, ES as values # CHANGE NES to ES
  heatmap_data <- significant_results %>%
    dplyr::select(pathway, ranked_list_name, NES) %>% 
    tidyr::pivot_wider(names_from = ranked_list_name, values_from = NES, values_fill = 0) # Fill non-significant/missing with 0
  
  mat_for_filtering <- as.matrix(heatmap_data %>% dplyr::select(-pathway))
  rownames(mat_for_filtering) <- heatmap_data$pathway
  
  # Filter for pathways with at least 5 non-zero NES values to avoid overly sparse heatmaps
  row_non_zero_nes_counts <- rowSums(mat_for_filtering != 0, na.rm = TRUE) # Count non-zero NES
  heatmap_data <- heatmap_data[row_non_zero_nes_counts >= 5, ]
  
  if (nrow(heatmap_data) == 0) {
    message(paste("  No pathways left after filtering for at least 5 positive ES enrichments for heatmap in", analysis_title_prefix, "with", geneset_source_name, geneset_direction, "gene sets. Skipping heatmap."))
    return(invisible(NULL))
  }
  
  # Convert to matrix, setting row names
  mat <- as.matrix(heatmap_data %>% dplyr::select(-pathway))
  rownames(mat) <- heatmap_data$pathway
  
  # Determine a symmetric color range based on max absolute NES
  max_abs_nes <- max(abs(mat), na.rm = TRUE) # CHANGE NES to ES
  # Define the color function using circlize::colorRamp2 for a diverging palette
  col_fun <- colorRamp2(c(-max_abs_nes, -max_abs_nes/2, 0, max_abs_nes/2, max_abs_nes), 
                        c("darkblue", "lightblue", "white", "pink2", "darkred")) # Matches her code's color scheme
  
  # Construct the plot title and filename
  hm_title <- paste0("Clustered Heatmap: ", analysis_title_prefix, "\n(", geneset_source_name, " ", geneset_direction, " Gene Sets vs. Ranked Lists, NES)",
                     "\nRows: ", nrow(mat), ", Cols: ", ncol(mat)) # Added row/column counts
  file_name_base <- paste0(gsub(" ", "_", analysis_title_prefix), "_", geneset_source_name, "_", geneset_direction, "_heatmap_NES") # Added NES to filename
  heatmap_output_path <- file.path(output_dir, paste0(file_name_base, ".png"))
  heatmap_output_path_svg <- file.path(output_dir, paste0(file_name_base, ".svg"))
  
  # Generate heatmap using ComplexHeatmap
  hm <- Heatmap(
    mat,
    name = "NES", # Legend name
    col = col_fun,
    na_col = "grey90", # Color for NA values
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 90
  )
  
  # Save as PNG
  png(heatmap_output_path, width = 2200, height = 1800, res = 300)
  draw(hm, column_title = hm_title)
  dev.off()
  
  # Save as SVG (vector file)
  svg(heatmap_output_path_svg, width = 7.33, height = 6) # width and height in inches
  draw(hm, column_title = hm_title)
  dev.off()
  
  message("  Clustered heatmap generated and saved to: ", heatmap_output_path)
}


#' Calculates combined NES and ES scores (NES_UP - NES_DN or ES_UP - ES_DN)
#' and a combined p-value (using Fisher's method) for fgsea results.
#'
#' @param fgsea_df A data frame of combined fgsea results (containing 'pathway', 'ranked_list_name', 'geneset_direction', 'ES', 'padj').
#' @param analysis_name A string for logging/message purposes.
#' @return A data frame with combined NES, combined ES, combined p-values, and original metadata.
calculate_combined_scores <- function(fgsea_df, analysis_name) {
  if (is.null(fgsea_df) || nrow(fgsea_df) == 0) {
    message(paste("No fgsea results to calculate combined scores for", analysis_name))
    return(NULL)
  }
  
  message(paste0("Calculating combined NES/ES and p-values for '", analysis_name, "'...")) 
  
  required_cols <- c("pathway", "ranked_list_name", "geneset_direction", "NES", "ES", "padj", "size", "geneset_source_name", "ranked_source_name")
  if (!all(required_cols %in% colnames(fgsea_df))) {
    stop(paste("Missing required columns for combined score calculation:", setdiff(required_cols, colnames(fgsea_df))))
  }
  
  combined_df <- fgsea_df %>%
    dplyr::select(pathway, ranked_list_name, geneset_direction, NES, ES, padj, size, geneset_source_name, ranked_source_name) %>%
    dplyr::mutate(geneset_direction = factor(geneset_direction, levels = c("UP", "DN"))) %>% 
    tidyr::pivot_wider(
      id_cols = c(pathway, ranked_list_name, geneset_source_name, ranked_source_name), 
      names_from = geneset_direction,
      values_from = c(NES, ES, padj, size),
      names_glue = "{.value}_{.name}",
      values_fill = NA_real_ 
    ) %>%
    # Use rename_with to conditionally rename if the double prefix exists (e.g., NES_NES_UP -> NES_UP)
    dplyr::rename_with(.fn = ~gsub("NES_NES_", "NES_", .x), .cols = dplyr::starts_with("NES_NES_")) %>%
    dplyr::rename_with(.fn = ~gsub("ES_ES_", "ES_", .x), .cols = dplyr::starts_with("ES_ES_")) %>%
    dplyr::rename_with(.fn = ~gsub("padj_padj_", "padj_", .x), .cols = dplyr::starts_with("padj_padj_")) %>%
    dplyr::rename_with(.fn = ~gsub("size_size_", "size_", .x), .cols = dplyr::starts_with("size_size_")) %>%
    dplyr::mutate(
      # Fill NAs for UP/DN NES, ES, and Size with 0 for calculation to allow subtraction
      NES_UP = if_else(is.na(NES_UP), 0, NES_UP),
      NES_DN = if_else(is.na(NES_DN), 0, NES_DN),
      ES_UP = if_else(is.na(ES_UP), 0, ES_UP),
      ES_DN = if_else(is.na(ES_DN), 0, ES_DN),
      size_UP = if_else(is.na(size_UP), 0, size_UP), # Default missing sizes to 0
      size_DN = if_else(is.na(size_DN), 0, size_DN), # Default missing sizes to 0
      
      # Classification label
      classification = case_when(
        (padj_UP < adj_pval_filter_val & padj_DN < adj_pval_filter_val & NES_UP > 0 & NES_DN < 0) ~ "Gero-Advancer",
        (padj_UP < adj_pval_filter_val & (is.na(padj_DN) | padj_DN >= adj_pval_filter_val) & NES_UP > 0 & NES_DN < 0) ~ "Partial UP Gero-Advancer",
        (padj_DN < adj_pval_filter_val & (is.na(padj_UP) | padj_UP >= adj_pval_filter_val) & NES_UP > 0 & NES_DN < 0) ~ "Partial DN Gero-Advancer",
        (padj_UP < adj_pval_filter_val & padj_DN < adj_pval_filter_val & NES_UP < 0 & NES_DN > 0) ~ "Gero-Protector",
        (padj_UP < adj_pval_filter_val & (is.na(padj_DN) | padj_DN >= adj_pval_filter_val) & NES_UP < 0 & NES_DN > 0) ~ "Partial UP Gero-Protector",
        (padj_DN < adj_pval_filter_val & (is.na(padj_UP) | padj_UP >= adj_pval_filter_val) & NES_UP < 0 & NES_DN > 0) ~ "Partial DN Gero-Protector",
        TRUE ~ "Not Significant"
      ),
      
      ## Compute combined scores
      combined_NES = NES_UP - NES_DN,
      # combined_ES = ES_UP - ES_DN,
      # Take the size, prioritizing UP if available, otherwise DN. If both 0/NA, size is 0.
      size = if_else(size_UP > 0, size_UP, size_DN),
      size = if_else(size == 0 & size_UP == 0 & size_DN == 0, 0, size), # Ensure it's 0 if no valid size
      
      # Combined P-val
      combined_padj = purrr::pmap_dbl(list(padj_UP, padj_DN, NES_UP, NES_DN, classification), function(p_up, p_dn, nes_up, nes_dn, class) {
        p_values_to_combine <- c()
        
        
        if (class == "Gero-Advancer" || class == "Gero-Protector") {
          p_values_to_combine <- c(p_up, p_dn) 
        } 
        else if (class == "Partial UP Gero-Advancer" || class == "Partial UP Gero-Protector") {
          p_values_to_combine <- p_up
        }
        else if (class == "Partial DN Gero-Advancer" || class == "Partial DN Gero-Protector") {
          p_values_to_combine <- p_dn
        }
  
        ## Combined p value estimation
        # Set to 1 if no p-values meet criteria or are available
        if (length(p_values_to_combine) == 0) {
          return(1) 
        } else if(length(p_values_to_combine) == 1) {
          return(p_values_to_combine[1])
        } else {
          # Both p-values available
          # Ensure p-values are not exactly zero for log calculation
          p_values_to_combine <- pmax(p_values_to_combine, .Machine$double.xmin)
        
          # Compute Fisher
          chisq <- -2 * sum(log(p_values_to_combine))
          df_fisher <- 2 * length(p_values_to_combine)
        
          return(stats::pchisq(chisq, df = df_fisher, lower.tail = FALSE))
        }
        })
      ) %>%
    # Select final desired columns
    dplyr::select(
      pathway, ranked_list_name, geneset_source_name, ranked_source_name,
      NES_UP, padj_UP, NES_DN, padj_DN,
      #ES_UP, ES_DN,combined_ES,
      combined_NES, combined_padj, size, classification
    ) %>%
    # Filter out rows where combined_padj could not be calculated (e.g., if both padj_UP and padj_DN were NA and didn't meet criteria)
    dplyr::filter(!is.na(combined_padj)) %>%
    # Annotate pathways based on combined p-value
    dplyr::mutate(
      combined_padj_status = case_when(
        combined_padj == 0 ~ "p=0", # Handle exact zero p-values
        combined_padj < 0.05 ~ "Significant",
        TRUE ~ "Not Significant"
      ),
      # Factor the status for consistent ordering in plots
      combined_padj_status = factor(combined_padj_status, levels = c("Significant", "p=0", "Not Significant"))
    )
  
  message(paste0("Combined scores calculated and annotated for '", analysis_name, "'. Total combined results: ", nrow(combined_df), " rows."))
  return(combined_df)
}


#' Generates and saves a dot plot for combined fgsea results (ES_UP - ES_DN).
#'
#' @param combined_fgsea_df A data frame of combined fgsea results from calculate_combined_scores.
#' @param analysis_title_prefix A string for plot titles, e.g., "Age-Centered Analysis".
#' @param output_dir Path to save the plots.
#' @param top_n Numeric, number of top results to plot.
#' @param essential_gene_list Character vector of essential genes for highlighting.
plot_combined_fgsea_dotplot <- function(combined_fgsea_df, analysis_title_prefix, output_dir, top_n = 20, essential_gene_list = NULL) {
  
  if (is.null(combined_fgsea_df) || nrow(combined_fgsea_df) == 0) {
    message(paste("No combined fgsea results for", analysis_title_prefix, "to display or plot."))
    return(invisible(NULL))
  }
  
  message(paste0("Generating combined dot plot for: ", analysis_title_prefix))
  
  all_significant_combined_results <- combined_fgsea_df %>%
    dplyr::filter(combined_padj < 0.05)
  
  if (nrow(all_significant_combined_results) == 0) {
    message(paste("No significant combined interactions (combined_padj < 0.05) found for", analysis_title_prefix, "to plot."))
    return(invisible(NULL))
  } else {
    message(paste0("Found ", nrow(all_significant_combined_results), " significant combined interactions (combined_padj < 0.05) for ", analysis_title_prefix, "."))
  }
  
  gene_set_label <- unique(all_significant_combined_results$geneset_source_name)[1]
  if (is.null(gene_set_label) || is.na(gene_set_label)) {
    gene_set_label <- "Gene Set" 
  } else {
    gene_set_label <- paste0(gene_set_label, " Gene Set")
  }
  
  ranked_list_label <- unique(all_significant_combined_results$ranked_source_name)[1]
  if (is.null(ranked_list_label) || is.na(ranked_list_label)) {
    ranked_list_label <- "Ranked List" 
  } else {
    ranked_list_label <- paste0(ranked_list_label, " Ranked List")
  }
  
  # Extract perturbation gene symbol and label essential genes
  plot_data_base <- all_significant_combined_results %>%
    dplyr::mutate(
      perturb_gene_symbol = gsub("^(.*?)\\s+Knockdown Signature - .*", "\\1", pathway) # Extract gene symbol from pathway name
    )
  
  if (!is.null(essential_gene_list) && length(essential_gene_list) > 0) {
    plot_data_base <- plot_data_base %>%
      dplyr::mutate(
        is_essential = perturb_gene_symbol %in% essential_gene_list
      )
  } else {
    plot_data_base$is_essential <- FALSE # Default to FALSE if no list provided
  }
  
  # --- Plot 1: Top Positive combined_NES (Gero-advancer) ---
  data_pos <- plot_data_base %>%
    dplyr::filter(combined_NES > 0) %>%
    dplyr::arrange(dplyr::desc(combined_NES)) %>%
    dplyr::slice(1:min(dplyr::n(), top_n))
  
  plot_pos <- if (nrow(data_pos) > 0) {
    ggplot(data_pos, aes(x = reorder(pathway, combined_NES), y = ranked_list_name)) +
      geom_point(aes(size = size, color = combined_NES, shape = is_essential)) + # 'is_essential' is now available
      scale_size_continuous(name = "Gene Set Size") +
      scale_color_gradient(low = "yellow", high = "red", name = "Combined NES\n(NES_UP - NES_DN)") +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), name = "Essential Gene", labels = c("No", "Yes")) +
      coord_flip() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
            axis.title = element_text(size = 8), axis.text = element_text(size = 7),
            legend.text = element_text(size = 7), legend.title = element_text(size = 8),
            legend.position = "bottom", plot.margin = margin(5, 5, 5, 5, "pt")) +
      labs(x = gene_set_label, y = ranked_list_label, title = paste0("Top ", top_n, " Gero-Advancer Pathways (Positive Combined NES)"))
  } else {
    ggplot() + geom_text(aes(x=0.5, y=0.5, label="No significant Gero-Advancer pathways"), size=4, color="grey50") + theme_void()
  }
  message(paste0("Prepared plot for Top ", top_n, " Gero-Advancer pathways (", nrow(data_pos), " results)."))
  
  # --- Plot 2: Top Negative combined_NES (Gero-protector) ---
  data_neg <- plot_data_base %>%
    dplyr::filter(combined_NES < 0) %>%
    dplyr::arrange(combined_NES) %>% # Arrange ascending for most negative first
    dplyr::slice(1:min(dplyr::n(), top_n))
  
  plot_neg <- if (nrow(data_neg) > 0) {
    ggplot(data_neg, aes(x = reorder(pathway, combined_NES), y = ranked_list_name)) +
      geom_point(aes(size = size, color = combined_NES, shape = is_essential)) + # 'is_essential' is now available
      scale_size_continuous(name = "Gene Set Size") +
      scale_color_gradient(low = "darkblue", high = "lightblue", name = "Combined NES\n(NES_UP - NES_DN)") + # Gero-protector colors
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), name = "Essential Gene", labels = c("No", "Yes")) +
      coord_flip() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
            axis.title = element_text(size = 8), axis.text = element_text(size = 7),
            legend.text = element_text(size = 7), legend.title = element_text(size = 8),
            legend.position = "bottom", plot.margin = margin(5, 5, 5, 5, "pt")) +
      labs(x = gene_set_label, y = ranked_list_label, title = paste0("Top ", top_n, " Gero-Protector Pathways (Negative Combined NES)"))
  } else {
    ggplot() + geom_text(aes(x=0.5, y=0.5, label="No significant Gero-Protector pathways"), size=4, color="grey50") + theme_void()
  }
  message(paste0("Prepared plot for Top ", top_n, " Gero-Protector pathways (", nrow(data_neg), " results)."))
  
  # Combine and save
  if (nrow(data_pos) > 0 || nrow(data_neg) > 0) {
    combined_plot <- plot_grid(plot_pos, plot_neg, ncol = 2, align = "hv", 
                               labels = c("A", "B"), label_size = 10)
    
    final_title_text <- paste0("GSEA Combined NES: ", analysis_title_prefix, " (Top ", top_n, " Pathways)") # Changed Combined ES to Combined NES
    final_title <- ggdraw() + 
      draw_label(final_title_text, fontface = 'bold', size = 16, x = 0.02, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7, "pt"))
    
    combined_plot_with_title <- plot_grid(final_title, combined_plot, ncol = 1, rel_heights = c(0.05, 1))
    
    plot_filename_png <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_top_", top_n, "_gsea_combined_NES_dotplot.png")) 
    ggsave(plot_filename_png, combined_plot_with_title, width = 16, height = 10)
    
    plot_filename_svg <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_top_", top_n, "_gsea_combined_NES_dotplot.svg")) 
    ggsave(plot_filename_svg, combined_plot_with_title, width = 16, height = 10)
    
    message("  Combined NES dot plot generated and saved to: ", plot_filename_png) 
  } else {
    message("No combined dot plots generated for ", analysis_title_prefix, " due to lack of significant results.")
  }
}

#' Generates and saves a clustered heatmap from combined fgsea results (combined_NES).
#' Highlights essential perturbation genes as a column annotation.
#'
#' @param combined_fgsea_df A data frame of combined fgsea results (from calculate_combined_scores).
#' @param analysis_title_prefix A string for plot titles, e.g., "Age-Centered Analysis".
#' @param output_dir Path to save the plots.
#' @param essential_gene_list Character vector of essential genes for highlighting.
#' @param filter Logical, if TRUE, filters out essential genes from columns.
#' @param ranked_source_name A string, e.g., "Aging" or "Perturbation", indicating the source of the ranked lists.
#' @return A list containing the heatmap's row dendrogram (`row_dend`) and column dendrogram (`col_dend`) if generated, else a list of NULLs.
plot_combined_heatmap <- function(combined_fgsea_df, analysis_title_prefix, output_dir, essential_gene_list = NULL, filter = FALSE, ranked_source_name) {
  
  if (is.null(combined_fgsea_df) || nrow(combined_fgsea_df) == 0) {
    message(paste("No combined fgsea results for heatmap in", analysis_title_prefix, "."))
    return(list(row_dend = NULL, col_dend = NULL))
  }
  
  message(paste0("Generating clustered heatmap for combined NES: ", analysis_title_prefix)) 
  
  # Filter for significant results (combined_padj < 0.05)
  significant_combined_results <- combined_fgsea_df %>%
    dplyr::filter(combined_padj < 0.05)
  
  if (nrow(significant_combined_results) == 0) {
    message(paste("  No significant combined interactions (combined_padj < 0.05) found for heatmap in", analysis_title_prefix, "."))
    return(list(row_dend = NULL, col_dend = NULL))
  }
  
  num_pathways <- length(unique(significant_combined_results$pathway))
  num_ranked_lists <- length(unique(significant_combined_results$ranked_list_name))
  
  if (num_pathways < 2 || num_ranked_lists < 2) {
    message(paste("  Not enough unique pathways (", num_pathways, ") or ranked lists (", num_ranked_lists, ") for a meaningful combined heatmap in", analysis_title_prefix, ". Skipping heatmap."))
    return(list(row_dend = NULL, col_dend = NULL))
  }
  
  # Reshape data for heatmap: combined_NES as values
  heatmap_data <- significant_combined_results %>%
    dplyr::select(pathway, ranked_list_name, combined_NES) %>%
    tidyr::pivot_wider(names_from = ranked_list_name, values_from = combined_NES, values_fill = 0) 
  
  # Filter heatmap_data rows where combined_NES != 0 for at least 5 columns (similar to original logic)
  mat_for_filtering <- as.matrix(heatmap_data %>% dplyr::select(-pathway))
  rownames(mat_for_filtering) <- heatmap_data$pathway
  row_non_zero_nes_counts <- rowSums(mat_for_filtering != 0, na.rm = TRUE)
  heatmap_data <- heatmap_data[row_non_zero_nes_counts >= 5, ] 
  
  if (nrow(heatmap_data) == 0) {
    message(paste("  No pathways left after filtering for at least 5 non-zero combined NES enrichments for heatmap in", analysis_title_prefix, ". Skipping heatmap."))
    return(list(row_dend = NULL, col_dend = NULL))
  }
  
  mat <- as.matrix(heatmap_data %>% dplyr::select(-pathway))
  rownames(mat) <- heatmap_data$pathway
  
  # Determine which axis represents perturbation genes for annotation
  # In Age-Centered analysis, perturbation signatures are pathways (rows).
  # In Perturbation-Centered analysis, perturbation signatures are ranked_list_name (columns).
  annotate_on_rows <- (ranked_source_name == "Aging")
  
  # Determine the gene names for initial essentiality check based on the annotation axis
  if (annotate_on_rows) {
    gene_names_for_annotation_initial <- rownames(mat)
  } else { # Annotate on columns
    gene_names_for_annotation_initial <- colnames(mat)
  }
  
  # Ensure essential_gene_list is a character vector for robust %in% operation
  if (!is.null(essential_gene_list)) {
    essential_gene_list <- as.character(essential_gene_list)
  }
  
  initial_perturb_gene_symbols <- gsub("^(.*?)\\s+Knockdown Signature - .*", "\\1", gene_names_for_annotation_initial)
  initial_is_essential_vector <- (initial_perturb_gene_symbols %in% essential_gene_list)
  
  # Apply filtering to 'mat' if filter is TRUE and essential_gene_list is provided
  if (filter == TRUE && !is.null(essential_gene_list) && length(essential_gene_list) > 0) {
    if (annotate_on_rows) {
      mat <- mat[!initial_is_essential_vector, , drop = FALSE]
    } else { # Filter columns
      mat <- mat[, !initial_is_essential_vector, drop = FALSE]
    }
    
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      message(paste("  No non-essential perturbation genes remaining after filtering for heatmap in", analysis_title_prefix, ". Skipping heatmap."))
      return(list(row_dend = NULL, col_dend = NULL))
    }
  }
  
  # After (potential) filtering, determine the final gene names for annotation
  if (annotate_on_rows) {
    final_gene_names_for_annotation <- rownames(mat)
  } else {
    final_gene_names_for_annotation <- colnames(mat)
  }
  
  final_perturb_gene_symbols <- gsub("^(.*?)\\s+Knockdown Signature - .*", "\\1", final_gene_names_for_annotation)
  final_is_essential_vector <- (final_perturb_gene_symbols %in% essential_gene_list)
  
  # Determine symmetric color range based on max absolute combined_NES
  max_abs_nes <- max(abs(mat), na.rm = TRUE) 
  col_fun <- colorRamp2(c(-max_abs_nes, -max_abs_nes/2, 0, max_abs_nes/2, max_abs_nes), 
                        c("darkblue", "lightblue", "white", "pink2", "darkred"))
  
  annotation_obj <- NULL
  # Create annotation_obj ONLY if essential_gene_list was provided AND is not empty
  if (!is.null(essential_gene_list) && length(essential_gene_list) > 0) {
    # Explicit check for annotation length consistency
    actual_length <- if (annotate_on_rows) nrow(mat) else ncol(mat)
    if (length(final_is_essential_vector) != actual_length) {
      stop(paste0("Annotation length mismatch: Annotation vector has ", length(final_is_essential_vector), 
                  " elements, but matrix has ", actual_length, " ", ifelse(annotate_on_rows, "rows", "columns"), ". This should not happen."))
    }
    
    is_ess_chr <- ifelse(final_is_essential_vector, "TRUE", "FALSE")
    ann_colors <- c("FALSE" = "grey90", "TRUE" = "darkgreen")
    
    if (annotate_on_rows) {
      annotation_obj <- rowAnnotation(
        is_essential = is_ess_chr,
        col = list(is_essential = ann_colors),
        annotation_name_side = "bottom", # Place row annotation name at the bottom
        width = unit(5, "mm") # Adjust width for better appearance
      )
    } else {
      annotation_obj <- HeatmapAnnotation(
        is_essential = is_ess_chr,
        col = list(is_essential = ann_colors),
        annotation_name_side = "left", # Place column annotation name at the left
        height = unit(5, "mm") # Adjust height for better appearance
      )
    }
  }
  
  # Title to reflect combined NES and include row/column counts
  hm_title <- paste0("Clustered Heatmap: ", analysis_title_prefix, "\n(Combined NES: NES_UP - NES_DN)",
                     "\nRows: ", nrow(mat), ", Cols: ", ncol(mat))
  file_name_base <- paste0(gsub(" ", "_", analysis_title_prefix), "_combined_NES_heatmap")
  heatmap_output_path_png <- file.path(output_dir, paste0(file_name_base, ".png"))
  heatmap_output_path_svg <- file.path(output_dir, paste0(file_name_base, ".svg"))
  
  # Heatmap call with conditional annotation placement
  if (annotate_on_rows) {
    hm <- Heatmap(
      mat,
      name = "Combined NES", 
      col = col_fun,
      na_col = "grey90",
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = FALSE,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_names_rot = 90,
      left_annotation = if (filter == TRUE) NULL else annotation_obj # Add annotation if not filtering
    )
  } else {
    hm <- Heatmap(
      mat,
      name = "Combined NES", 
      col = col_fun,
      na_col = "grey90",
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = FALSE,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_names_rot = 90,
      top_annotation = if (filter == TRUE) NULL else annotation_obj # Add annotation if not filtering
    )
  }
  
  # Draw the heatmap object to extract dendrograms
  ht_list <- draw(hm, column_title = hm_title)
  row_dend_res <- row_dend(ht_list)
  col_dend_res <- column_dend(ht_list)
  
  png(heatmap_output_path_png, width = 2200, height = 1800, res = 300)
  draw(ht_list, column_title = hm_title) # Redraw for saving purposes
  dev.off()
  
  svg(heatmap_output_path_svg, width = 7.33, height = 6)
  draw(ht_list, column_title = hm_title) # Redraw for saving purposes
  dev.off()
  
  message("  Combined NES clustered heatmap generated and saved to: ", heatmap_output_path_png)
  
  return(list(row_dend = row_dend_res, col_dend = col_dend_res)) # Return dendrograms
}


#' Generates and saves a clustered heatmap from combined fgsea results using a signed -log10(p-value) for coloring.
#' The sign is determined by combined_NES. Highlights essential perturbation genes as a column annotation.
#'
#' @param combined_fgsea_df A data frame of combined fgsea results (from calculate_combined_scores).
#' @param analysis_title_prefix A string for plot titles, e.g., "Age-Centered Analysis".
#' @param output_dir Path to save the plots.
#' @param essential_gene_list Character vector of essential genes for highlighting.
#' @param filter Logical, if TRUE, filters out essential genes from columns.
#' @param ranked_source_name A string, e.g., "Aging" or "Perturbation", indicating the source of the ranked lists.
#' @param row_dend Optional, dendrogram for row clustering (from NES heatmap).
#' @param col_dend Optional, dendrogram for column clustering (from NES heatmap).
plot_combined_heatmap_signed_pvalue_NES <- function(combined_fgsea_df, analysis_title_prefix, output_dir, essential_gene_list = NULL, filter = FALSE, ranked_source_name, row_dend = NULL, col_dend = NULL) {
  
  if (is.null(combined_fgsea_df) || nrow(combined_fgsea_df) == 0) {
    message(paste("No combined fgsea results for signed -log10(p-value) heatmap (NES-based) in", analysis_title_prefix, "."))
    return(invisible(NULL))
  }
  
  message(paste0("Generating clustered heatmap for signed -log10(p-value) (NES-based): ", analysis_title_prefix))
  
  significant_combined_results <- combined_fgsea_df %>%
    dplyr::filter(combined_padj < 0.05)
  
  if (nrow(significant_combined_results) == 0) {
    message(paste("  No significant combined interactions (combined_padj < 0.05) found for signed -log10(p-value) heatmap (NES-based) in", analysis_title_prefix, "."))
    return(invisible(NULL))
  }
  
  num_pathways <- length(unique(significant_combined_results$pathway))
  num_ranked_lists <- length(unique(significant_combined_results$ranked_list_name))
  
  if (num_pathways < 2 || num_ranked_lists < 2) {
    message(paste("  Not enough unique pathways (", num_pathways, ") or ranked lists (", num_ranked_lists, ") for a meaningful signed -log10(p-value) heatmap (NES-based) in", analysis_title_prefix, ". Skipping heatmap."))
    return(invisible(NULL))
  }
  
  # Calculate signed_log10_p using combined_NES and categorize p-values
  heatmap_data <- significant_combined_results %>%
    dplyr::mutate(
      signed_log10_p = sign(combined_NES) * (-log10(combined_padj)),
      # Discretize the combined_padj values into categories for coloring
      p_category = case_when(
        combined_padj < 0.001 & combined_NES > 0 ~ "p < 0.001 (UP)",
        combined_padj < 0.01 & combined_NES > 0 ~ "p < 0.01 (UP)",
        combined_padj >= 0.01 & combined_padj < 0.05 & combined_NES > 0 ~ "p < 0.05 (UP)",
        combined_padj < 0.001 & combined_NES < 0 ~ "p < 0.001 (DN)",
        combined_padj < 0.01 & combined_NES < 0 ~ "p < 0.01 (DN)",
        combined_padj >= 0.01 & combined_padj < 0.05 & combined_NES < 0 ~ "p < 0.05 (DN)",
        TRUE ~ "Not significant" # For combined_padj >= 0.05 or NES == 0
      ),
      # Ensure explicit factor levels for consistent legend order (DN -> Not Sig -> UP)
      p_category = factor(p_category, levels = c(
        "p < 0.001 (DN)",
        "p < 0.01 (DN)", 
        "p < 0.05 (DN)", 
        "Not significant", 
        "p < 0.05 (UP)", 
        "p < 0.01 (UP)",
        "p < 0.001 (UP)"
      )) 
    ) %>%
    dplyr::select(pathway, ranked_list_name, signed_log10_p, p_category) %>% # Keep p_category for coloring
    tidyr::pivot_wider(names_from = ranked_list_name, values_from = c(signed_log10_p, p_category), names_glue = "{.value}_{ranked_list_name}", values_fill = NA)
  
  # Extract matrices for signed_log10_p (for values/filtering) and p_category (for coloring)
  mat_signed_p_val <- as.matrix(heatmap_data %>% dplyr::select(starts_with("signed_log10_p_")))
  colnames(mat_signed_p_val) <- gsub("signed_log10_p_", "", colnames(mat_signed_p_val))
  rownames(mat_signed_p_val) <- heatmap_data$pathway
  
  mat_p_category <- as.matrix(heatmap_data %>% dplyr::select(starts_with("p_category_")))
  colnames(mat_p_category) <- gsub("p_category_", "", colnames(mat_p_category))
  rownames(mat_p_category) <- heatmap_data$pathway
  
  # Filter based on mat_signed_p_val (if values are 0, they are not considered enriched)
  mat_for_filtering <- mat_signed_p_val
  row_non_zero_scores_counts <- rowSums(mat_for_filtering != 0, na.rm = TRUE)
  
  # Apply the row filtering to both matrices
  mat_signed_p_val <- mat_signed_p_val[row_non_zero_scores_counts >= 5, , drop = FALSE]
  mat_p_category <- mat_p_category[row_non_zero_scores_counts >= 5, , drop = FALSE]
  
  if (nrow(mat_signed_p_val) == 0) {
    message(paste("  No pathways left after filtering for at least 5 non-zero signed -log10(p-value) enrichments (NES-based) for heatmap in", analysis_title_prefix, ". Skipping heatmap."))
    return(invisible(NULL))
  }
  
  # Determine which axis represents perturbation genes for annotation
  annotate_on_rows <- (ranked_source_name == "Aging")
  
  # Determine the gene names for initial essentiality check based on the annotation axis
  if (annotate_on_rows) {
    gene_names_for_annotation_initial <- rownames(mat_signed_p_val)
  } else { # Annotate on columns
    gene_names_for_annotation_initial <- colnames(mat_signed_p_val)
  }
  
  # Ensure essential_gene_list is a character vector for robust %in% operation
  if (!is.null(essential_gene_list)) {
    essential_gene_list <- as.character(essential_gene_list)
  }
  
  initial_perturb_gene_symbols <- gsub("^(.*?)\\s+Knockdown Signature - .*", "\\1", gene_names_for_annotation_initial)
  initial_is_essential_vector <- (initial_perturb_gene_symbols %in% essential_gene_list)
  
  # Apply filtering to matrices if filter is TRUE and essential_gene_list is provided
  if (filter == TRUE && !is.null(essential_gene_list) && length(essential_gene_list) > 0) {
    if (annotate_on_rows) {
      mat_signed_p_val <- mat_signed_p_val[!initial_is_essential_vector, , drop = FALSE]
      mat_p_category <- mat_p_category[!initial_is_essential_vector, , drop = FALSE]
    } else { # Filter columns
      mat_signed_p_val <- mat_signed_p_val[, !initial_is_essential_vector, drop = FALSE]
      mat_p_category <- mat_p_category[, !initial_is_essential_vector, drop = FALSE]
    }
    
    if (nrow(mat_signed_p_val) == 0 || ncol(mat_signed_p_val) == 0) {
      message(paste("  No non-essential perturbation genes remaining after filtering for heatmap in", analysis_title_prefix, ". Skipping heatmap."))
      return(invisible(NULL))
    }
  }
  
  # After (potential) filtering, determine the final gene names for annotation
  if (annotate_on_rows) {
    final_gene_names_for_annotation <- rownames(mat_signed_p_val)
  } else {
    final_gene_names_for_annotation <- colnames(mat_signed_p_val)
  }
  
  final_perturb_gene_symbols <- gsub("^(.*?)\\s+Knockdown Signature - .*", "\\1", final_gene_names_for_annotation)
  final_is_essential_vector <- (final_perturb_gene_symbols %in% essential_gene_list)
  
  # Define discrete color palette for p_category, matching volcano plot colors
  p_category_colors <- c(
    "p < 0.001 (DN)" = "darkblue",       # Highest significance DN
    "p < 0.01 (DN)" = "blue",       # Highest significance DN
    "p < 0.05 (DN)" = "lightblue",      # Medium significance DN
    "Not significant" = "grey90",       # No significance
    "p < 0.05 (UP)" = "lightcoral",     # Medium significance UP
    "p < 0.01 (UP)" = "red",         # Highest significance UP
    "p < 0.001 (UP)" = "darkred"         # Highest significance UP
  )
  
  annotation_obj <- NULL
  # Create annotation_obj ONLY if essential_gene_list was provided AND is not empty
  if (!is.null(essential_gene_list) && length(essential_gene_list) > 0) {
    # Check length consistency after filtering
    actual_length <- if (annotate_on_rows) nrow(mat_signed_p_val) else ncol(mat_signed_p_val)
    if (length(final_is_essential_vector) != actual_length) {
      stop(paste0("Annotation length mismatch: Annotation vector has ", length(final_is_essential_vector), 
                  " elements, but matrix has ", actual_length, " ", ifelse(annotate_on_rows, "rows", "columns"), ". This should not happen."))
    }
    
    is_ess_chr <- ifelse(final_is_essential_vector, "TRUE", "FALSE")
    ann_colors <- c("FALSE" = "grey90", "TRUE" = "darkgreen")
    
    if (annotate_on_rows) {
      annotation_obj <- rowAnnotation(
        is_essential = is_ess_chr,
        col = list(is_essential = ann_colors),
        annotation_name_side = "bottom", # Place row annotation name at the bottom
        width = unit(5, "mm") # Adjust width for better appearance
      )
    } else {
      annotation_obj <- HeatmapAnnotation(
        is_essential = is_ess_chr,
        col = list(is_essential = ann_colors),
        annotation_name_side = "left", # Place column annotation name at the left
        height = unit(5, "mm") # Adjust height for better appearance
      )
    }
  }
  
  hm_title <- paste0("Clustered Heatmap: ", analysis_title_prefix, "\n(Signed -log10(Combined P-value) Categories, by NES)", 
                     "\nRows: ", nrow(mat_signed_p_val), ", Cols: ", ncol(mat_signed_p_val)) 
  file_name_base <- paste0(gsub(" ", "_", analysis_title_prefix), "_signed_log10_pvalue_NES_categories_heatmap") 
  heatmap_output_path_png <- file.path(output_dir, paste0(file_name_base, ".png"))
  heatmap_output_path_svg <- file.path(output_dir, paste0(file_name_base, ".svg"))
  
  # Heatmap call with conditional annotation placement and clustering
  if (annotate_on_rows) {
    hm <- Heatmap(
      mat_p_category, # Use the categorical matrix for coloring
      name = "Combined P-value\nCategory", # Legend name reflects categories
      col = p_category_colors, # Use discrete color mapping
      na_col = "grey90",
      cluster_rows = row_dend, # Use provided row dendrogram
      cluster_columns = col_dend, # Use provided column dendrogram
      show_row_names = FALSE,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_names_rot = 90,
      left_annotation = if (filter == TRUE) NULL else annotation_obj # Add annotation if not filtering
    )
  } else {
    hm <- Heatmap(
      mat_p_category, # Use the categorical matrix for coloring
      name = "Combined P-value\nCategory", # Legend name reflects categories
      col = p_category_colors, # Use discrete color mapping
      na_col = "grey90",
      cluster_rows = row_dend, # Use provided row dendrogram
      cluster_columns = col_dend, # Use provided column dendrogram
      show_row_names = FALSE,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_names_rot = 90,
      top_annotation = if (filter == TRUE) NULL else annotation_obj # Add annotation if not filtering
    )
  }
  
  png(heatmap_output_path_png, width = 2200, height = 1800, res = 300)
  draw(hm, column_title = hm_title)
  dev.off()
  
  svg(heatmap_output_path_svg, width = 7.33, height = 6)
  draw(hm, column_title = hm_title)
  dev.off()
  
  message("  Signed -log10(p-value) (NES-based) clustered heatmap generated and saved to: ", heatmap_output_path_png)
}


#' Generates an exploratory plot of adjusted p-value vs. ES.
#'
#' @param fgsea_df A data frame of fgsea results (can be raw or combined, specify ES_col and pval_col).
#' @param analysis_title_prefix A string for plot titles.
#' @param output_dir Path to save the plots.
#' @param es_col The name of the ES column (e.g., "NES", "ES", "combined_NES").
#' @param pval_col The name of the p-value column (e.g., "padj", "combined_padj").
#' @param max_plot_points Optional: maximum number of points to plot for performance, samples if more.
plot_pval_vs_es <- function(fgsea_df, analysis_title_prefix, output_dir, 
                            es_col = "NES", pval_col = "padj", max_plot_points = 50000) {
  
  if (is.null(fgsea_df) || nrow(fgsea_df) == 0) {
    message(paste("No fgsea results for p-value vs ES plot in", analysis_title_prefix, "."))
    return(invisible(NULL))
  }
  
  message(paste0("Generating p-value vs ", es_col, " plot for: ", analysis_title_prefix))
  
  if (!all(c(es_col, pval_col) %in% colnames(fgsea_df))) {
    warning(paste("Missing required columns ('", es_col, "', '", pval_col, "') for p-value vs ES plot in", analysis_title_prefix, ". Skipping."))
    return(invisible(NULL))
  }
  
  plot_data <- fgsea_df %>%
    dplyr::filter(!is.na(!!sym(es_col)), !is.na(!!sym(pval_col))) %>%
    dplyr::mutate(log10_pval = -log10(!!sym(pval_col)))
  
  if (nrow(plot_data) == 0) {
    message(paste("No valid data points for p-value vs ES plot in", analysis_title_prefix, "after NA filtering. Skipping."))
    return(invisible(NULL))
  }
  
  if (nrow(plot_data) > max_plot_points) {
    message(paste0("  Sampling ", max_plot_points, " points for p-value vs ES plot to improve performance."))
    plot_data <- plot_data %>% dplyr::sample_n(max_plot_points)
  }
  
  p <- ggplot(plot_data, aes(x = !!sym(es_col), y = log10_pval)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(
      title = paste0(analysis_title_prefix, ": -log10(", pval_col, ") vs ", es_col),
      x = es_col,
      y = paste0("-log10(", pval_col, ")")
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  plot_filename_png <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_", pval_col, "_vs_", es_col, ".png"))
  ggsave(plot_filename_png, p, width = 8, height = 6)
  
  plot_filename_svg <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_", pval_col, "_vs_", es_col, ".svg"))
  ggsave(plot_filename_svg, p, width = 8, height = 6)
  
  message("  P-value vs ES plot generated and saved to: ", plot_filename_png)
}

#' Generates an exploratory plot of gene set size vs. ES.
#'
#' @param fgsea_df A data frame of fgsea results (can be raw or combined, specify ES_col).
#' @param analysis_title_prefix A string for plot titles.
#' @param output_dir Path to save the plots.
#' @param es_col The name of the ES column (e.g., "NES", "ES", "combined_NES").
#' @param max_plot_points Optional: maximum number of points to plot for performance, samples if more.
plot_geneset_size_vs_es <- function(fgsea_df, analysis_title_prefix, output_dir, 
                                    es_col = "NES", max_plot_points = 50000) {
  
  if (is.null(fgsea_df) || nrow(fgsea_df) == 0) {
    message(paste("No fgsea results for gene set size vs ES plot in", analysis_title_prefix, "."))
    return(invisible(NULL))
  }
  
  message(paste0("Generating gene set size vs ", es_col, " plot for: ", analysis_title_prefix))
  
  if (!all(c("size", es_col) %in% colnames(fgsea_df))) {
    warning(paste("Missing required columns ('size', '", es_col, "') for gene set size vs ES plot in", analysis_title_prefix, ". Skipping."))
    return(invisible(NULL))
  }
  
  plot_data <- fgsea_df %>%
    dplyr::filter(!is.na(size), !is.na(!!sym(es_col)))
  
  if (nrow(plot_data) == 0) {
    message(paste("No valid data points for gene set size vs ES plot in", analysis_title_prefix, "after NA filtering. Skipping."))
    return(invisible(NULL))
  }
  
  if (nrow(plot_data) > max_plot_points) {
    message(paste0("  Sampling ", max_plot_points, " points for gene set size vs ES plot to improve performance."))
    plot_data <- plot_data %>% dplyr::sample_n(max_plot_points)
  }
  
  p <- ggplot(plot_data, aes(x = size, y = !!sym(es_col))) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "loess", color = "blue", se = FALSE) + 
    labs(
      title = paste0(analysis_title_prefix, ": Gene Set Size vs ", es_col),
      x = "Gene Set Size",
      y = es_col
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  plot_filename_png <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_size_vs_", es_col, ".png"))
  ggsave(plot_filename_png, p, width = 8, height = 6)
  
  plot_filename_svg <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_size_vs_", es_col, ".svg"))
  ggsave(plot_filename_svg, p, width = 8, height = 6)
  
  message("  Gene set size vs ES plot generated and saved to: ", plot_filename_png)
}


#' Generates and saves a volcano plot of combined NES against adjusted p-value,
#' with points colored by specific categories based on individual and combined significance.
#'
#' @param combined_fgsea_df A data frame of combined fgsea results from calculate_combined_scores.
#' @param analysis_title_prefix A string for plot titles, e.g., "Age-Centered Analysis".
#' @param output_dir Path to save the plots.
#' @param indiv_padj_thresh Numeric, threshold for individual padj significance (e.g., 0.05).
plot_volcano_with_coloring <- function(combined_fgsea_df, analysis_title_prefix, output_dir, indiv_padj_thresh = 0.05) {
  
  if (is.null(combined_fgsea_df) || nrow(combined_fgsea_df) == 0) {
    message(paste("No combined fgsea results for volcano plot with coloring in", analysis_title_prefix, "."))
    return(invisible(NULL))
  }
  
  message(paste0("Generating volcano plot with coloring for: ", analysis_title_prefix))
  
  # Ensure necessary columns exist for plotting
  required_cols <- c("combined_NES", "combined_padj", "padj_UP", "padj_DN", "NES_UP", "NES_DN", "combined_padj_status")
  if (!all(required_cols %in% colnames(combined_fgsea_df))) {
    warning(paste("Missing required columns for volcano plot with coloring:", setdiff(required_cols, colnames(combined_fgsea_df)), ". Skipping plot."))
    return(invisible(NULL))
  }
  
  # Prepare plot_data based on user's provided logic
  plot_data <- combined_fgsea_df %>%
    dplyr::filter(!is.na(combined_NES), !is.na(combined_padj), combined_padj != 0) %>% # Filter out NA/zero p-values
    dplyr::mutate(
      log10_combined_padj = -log10(combined_padj),
      
      # Determine significance for individual UP/DN padj (handling NAs)
      padj_UP_sig = !is.na(padj_UP) & (padj_UP < indiv_padj_thresh),
      padj_DN_sig = !is.na(padj_DN) & (padj_DN < indiv_padj_thresh),
      
      # Define coloring_category with more interpretable labels
      coloring_category = case_when(
        # Both individual padj are significant, and ES have same sign (or one is zero)
        combined_padj_status == "Significant" & padj_UP_sig & padj_DN_sig & (sign(NES_UP) == sign(NES_DN) | NES_UP == 0 | NES_DN == 0) ~ "Highly Significant & Consistent", # Dark Red
        
        # Both individual padj are significant, and ES have opposite signs
        combined_padj_status == "Significant" & padj_UP_sig & padj_DN_sig & sign(NES_UP) != sign(NES_DN) ~ "Highly Significant & Bidirectional", # Dark Green
        
        # Combined significant, but one individual padj is NOT significant (or NA) and ES have same sign (or one is zero)
        combined_padj_status == "Significant" & ((!padj_UP_sig & !is.na(padj_UP)) | (!padj_DN_sig & !is.na(padj_DN))) & (sign(NES_UP) == sign(NES_DN) | NES_UP == 0 | NES_DN == 0) ~ "Significant (Weak Individual P)", # Light Red
        
        # Any other significant case not covered above (e.g. one-sided padj_UP_sig but NES_DN=0)
        combined_padj_status == "Significant" ~ "Significant (Other Cases)", # Purple
        
        TRUE ~ "Not Significant" # Default for non-significant combined_padj (Grey)
      ),
      # Ensure factor order for legend consistency
      coloring_category = factor(coloring_category, levels = c(
        "Highly Significant & Consistent",
        "Highly Significant & Bidirectional",
        "Significant (Weak Individual P)",
        "Significant (Other Cases)",
        "Not Significant"
      ))
    )
  
  # Define plot_colors
  plot_colors <- c(
    "Highly Significant & Consistent" = "darkred",
    "Highly Significant & Bidirectional" = "darkgreen",
    "Significant (Weak Individual P)" = "lightcoral",
    "Significant (Other Cases)" = "purple",
    "Not Significant" = "grey80"
  )
  
  if (nrow(plot_data) == 0) {
    message(paste("No valid data points for volcano plot with coloring in", analysis_title_prefix, "after filtering. Skipping."))
    return(invisible(NULL))
  }
  
  # Create the Volcano Plot
  volcano_plot <- ggplot(plot_data, aes(x = combined_NES, y = log10_combined_padj, color = coloring_category)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = plot_colors, name = "Combined Status") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # Combined padj threshold
    labs(
      title = paste0("Volcano Plot: ", analysis_title_prefix, " Perturbation-Aging Signature Analysis"),
      x = "Combined Normalized Enrichment Score (NES)",
      y = "-log10(Combined Adjusted p-value)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plot_filename_png <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_volcano_plot_colored.png"))
  ggsave(plot_filename_png, volcano_plot, width = 10, height = 8)
  
  plot_filename_svg <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_volcano_plot_colored.svg"))
  ggsave(plot_filename_svg, volcano_plot, width = 10, height = 8)
  
  message("  Volcano plot with coloring generated and saved to: ", plot_filename_png)
}


#' Generates scatter plots comparing combined NES between age-centered and perturbation-centered analyses.
#' Each point represents a unique (Perturbation_Gene, Cell_Line, Tissue) triplet.
#' The plots include a y=x reference line, and highlight essential perturbation genes.
#'
#' @param age_combined_df A data frame of combined fgsea results from age-centered analysis.
#' @param perturb_combined_df A data frame of combined fgsea results from perturbation-centered analysis.
#' @param output_dir Path to save the plots.
#' @param essential_gene_list Character vector of essential genes for highlighting.
#' @param top_n_to_label Numeric, number of top/bottom points to label by default.
plot_cross_analysis_scatterplots <- function(age_combined_df, perturb_combined_df, output_dir, essential_gene_list = NULL, top_n_to_label = 20, filter = TRUE) {
  message("Generating cross-analysis scatter plots (Combined NES)...") # Removed ES from message
  
  if (is.null(age_combined_df) || nrow(age_combined_df) == 0) {
    message("Age-centered combined results are empty. Skipping cross-analysis scatter plots.")
    return(invisible(NULL))
  }
  if (is.null(perturb_combined_df) || nrow(perturb_combined_df) == 0) {
    message("Perturbation-centered combined results are empty. Skipping cross-analysis scatter plots.")
    return(invisible(NULL))
  }
  
  # Filter for significant results using the new 'combined_padj_status'
  age_filt <- age_combined_df %>%
    {
      if (filter) {
        dplyr::filter(., combined_padj_status == 'Significant')
      } else {
        . # If 'filter' is FALSE, pass the data through unchanged
      }
    } %>%
    dplyr::rename(
      combined_ES_age = combined_ES,
      combined_NES_age = combined_NES,
      combined_padj_age = combined_padj,
      size_age = size
    )
  
  perturb_filt <- perturb_combined_df %>%
    {
      if (filter) {
        dplyr::filter(., combined_padj_status == 'Significant')
      } else {
        . # If 'filter' is FALSE, pass the data through unchanged
      }
    } %>%
    dplyr::rename(
      combined_ES_perturb = combined_ES,
      combined_NES_perturb = combined_NES,
      combined_padj_perturb = combined_padj,
      size_perturb = size
    )
  
  if (nrow(age_filt) == 0 || nrow(perturb_filt) == 0) {
    message("No significant data after filtering for cross-analysis scatter plots. Skipping.")
    return(invisible(NULL))
  }
  
  # Helper function to extract base name and direction
  extract_components_vectorized <- function(signature_strings) {
    direction <- rep("UNKNOWN", length(signature_strings))
    base_name <- signature_strings
    
    # Identify entries ending with _UP or _DN
    up_idx <- grepl("_UP$", signature_strings)
    dn_idx <- grepl("_DN$", signature_strings)
    
    # Assign direction and remove suffix for _UP entries
    direction[up_idx] <- "UP"
    base_name[up_idx] <- gsub("_UP$", "", signature_strings[up_idx])
    
    # Assign direction and remove suffix for _DN entries
    direction[dn_idx] <- "DN"
    base_name[dn_idx] <- gsub("_DN$", "", signature_strings[dn_idx])
    
    list(base = base_name, direction = direction)
  }
  
  # Create symmetric identifiers for joining the two dataframes.
  # The goal is to join on (Perturbation_Signature_Base, Aging_Tissue_Base, Direction).
  
  age_filt_for_join <- age_filt %>%
    dplyr::mutate(
      # 'pathway' column contains the perturbation signature with _UP/_DN suffix
      # Directly extract base and direction
      perturb_base_signature = extract_components_vectorized(pathway)$base,
      direction = extract_components_vectorized(pathway)$direction,
      
      # 'ranked_list_name' column contains the clean aging tissue name
      aging_tissue = ranked_list_name,
      
      # Create a consistent identifier string for joining
      identifier = paste0(perturb_base_signature, "__", aging_tissue, "__", direction)
    ) %>%
    dplyr::select(identifier, combined_ES_age, combined_NES_age, combined_padj_age, size_age)
  
  perturb_filt_for_join <- perturb_filt %>%
    dplyr::mutate(
      # 'ranked_list_name' column contains the clean perturbation signature name
      perturb_base_signature = ranked_list_name,
      
      # 'pathway' column contains the aging tissue name with _UP/_DN suffix
      # Directly extract base and direction
      aging_tissue = extract_components_vectorized(pathway)$base,
      direction = extract_components_vectorized(pathway)$direction,
      
      # Create a consistent identifier string for joining
      identifier = paste0(perturb_base_signature, "__", aging_tissue, "__", direction)
    ) %>%
    dplyr::select(identifier, combined_ES_perturb, combined_NES_perturb, combined_padj_perturb, size_perturb)
  
  merged_data <- dplyr::inner_join(
    age_filt_for_join,
    perturb_filt_for_join,
    by = "identifier"
  )
  
  if (nrow(merged_data) == 0) {
    message("No common significant (perturbation_signature, aging_tissue) pairs found for cross-analysis scatter plots. Skipping.")
    return(invisible(NULL))
  }
  
  # Add metadata for plotting (perturb_gene_cl, tissue, direction, perturb_gene_symbol, is_essential)
  merged_data <- merged_data %>%
    # Extract perturb_gene_cl, tissue, and direction back from the identifier
    tidyr::separate(identifier, into = c("perturb_gene_cl", "tissue", "direction"), sep = "__", remove = FALSE) %>%
    dplyr::mutate(
      # Label for text annotation, now including direction
      label = paste0(perturb_gene_cl, " (", direction, "), Tissue: ", tissue),
      # Extract gene symbol (e.g., "UBL5" from "UBL5 Knockdown Signature - K562")
      perturb_gene_symbol = gsub("^(.*?)\\s+Knockdown Signature - .*", "\\1", perturb_gene_cl)
    )
  
  # Label essential genes
  if (!is.null(essential_gene_list) && length(essential_gene_list) > 0) {
    merged_data <- merged_data %>%
      dplyr::mutate(
        is_essential = perturb_gene_symbol %in% essential_gene_list
      )
  } else {
    merged_data$is_essential <- FALSE # Default to FALSE if no list provided
  }
  
  # Identify points to label (top/bottom correlated points)
  # Using product of NES for correlation, similar to your original scatterplot logic
  top_cor <- merged_data %>%
    dplyr::arrange(dplyr::desc(combined_NES_age * combined_NES_perturb)) %>%
    dplyr::slice(1:min(dplyr::n(), top_n_to_label))
  
  bottom_cor <- merged_data %>%
    dplyr::arrange(combined_NES_age * combined_NES_perturb) %>%
    dplyr::slice(1:min(dplyr::n(), top_n_to_label))
  
  labels_to_show <- unique(rbind(top_cor, bottom_cor)) %>% dplyr::pull(label)
  
  plot_data_final <- merged_data %>%
    dplyr::mutate(is_labeled = label %in% labels_to_show)
  
  if (filter) {
    title = "Comparison of Combined NES: Perturbation-Centered vs. Age-Centered Analysis"
  } else {
    title = "Comparison of Combined NES: Perturbation-Centered vs. Age-Centered Analysis, non-filtered"
  }
  
  # --- Plot 1: Combined NES Comparison ---
  plot_nes_comparison <- ggplot(plot_data_final, aes(x = combined_NES_perturb, y = combined_NES_age)) +
    geom_point(aes(color = combined_NES_perturb, size = size_perturb, shape = is_essential), alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") + # Add y=x line
    geom_text(data = dplyr::filter(plot_data_final, is_labeled), aes(label = label),
              size = 2.5, vjust = -0.8, hjust = 0.5, check_overlap = TRUE) +
    scale_color_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0, name = "Perturb. NES") + # Symmetric colors
    scale_size_continuous(range = c(1, 5), name = "Perturb. Geneset Size") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), name = "Essential Perturbation", labels = c("No", "Yes")) +
    labs(
      title = title,
      subtitle = "Each point: (Perturbation Gene + Cell Line, Tissue) pair; Shape: Essential Gene",
      x = "Perturbation-Centered Combined NES (Perturbation's effect on Aging)",
      y = "Age-Centered Combined NES (Aging's effect on Perturbation)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  if (filter) {
    plot_filename_nes_png <- file.path(output_dir, "cross_analysis_combined_NES_scatterplot.png")
    plot_filename_nes_svg <- file.path(output_dir, "cross_analysis_combined_NES_scatterplot.svg")
  } else {
    plot_filename_nes_png <- file.path(output_dir, "cross_analysis_combined_NES_scatterplot_non-filtered.png")
    plot_filename_nes_svg <- file.path(output_dir, "cross_analysis_combined_NES_scatterplot_non-filtered.svg")
  }
  ggsave(plot_filename_nes_png, plot_nes_comparison, width = 12, height = 10)
  ggsave(plot_filename_nes_svg, plot_nes_comparison, width = 12, height = 10)
  message("  Cross-analysis Combined NES scatter plot generated and saved to: ", plot_filename_nes_png)
}



#' Performs Fisher's exact test to check for enrichment of essential genes in
#' geroadvancer or geroprotector categories.
#'
#' @param combined_fgsea_df A data frame of combined fgsea results from calculate_combined_scores.
#' @param essential_gene_list Character vector of essential genes.
#' @param analysis_name A string for logging/message purposes.
#' @param nr_all_depmap_genes Numeric, total unique genes present in the full DepMap summary file.
#' @param nr_ess_depmap_genes Numeric, count of common essential genes from DepMap.
#' @return A list of lists, containing results for geroadvancers and geroprotectors.
perform_essential_gene_enrichment_test <- function(combined_fgsea_df, essential_gene_list, analysis_name, nr_all_depmap_genes, nr_ess_depmap_genes) {
  message(paste0("\n--- Performing Essential Gene Enrichment Test for '", analysis_name, "' ---"))
  message(paste0("  DepMap Reference: Total unique genes in DepMap summary: ", nr_all_depmap_genes, "; Common essential genes: ", nr_ess_depmap_genes, "."))
  
  if (is.null(combined_fgsea_df) || nrow(combined_fgsea_df) == 0) {
    message(paste("No combined fgsea results for essential gene enrichment in", analysis_name, ". Skipping."))
    return(NULL)
  }
  
  if (is.null(essential_gene_list) || length(essential_gene_list) == 0) {
    message("Essential gene list is empty. Skipping essential gene enrichment test.")
    return(NULL)
  }
  
  # Extract all unique perturbation gene symbols from the pathways present in the combined_fgsea_df.
  # This serves as the universe of perturbation genes considered in this specific analysis run.
  all_perturb_gene_symbols_in_analysis <- unique(gsub("^(.*?)\\s+Knockdown Signature - .*", "\\1", combined_fgsea_df$pathway))
  
  if (length(all_perturb_gene_symbols_in_analysis) == 0) {
    message("No perturbation gene symbols could be extracted from pathways. Skipping essential gene enrichment test.")
    return(NULL)
  }
  
  total_genes_in_universe <- length(all_perturb_gene_symbols_in_analysis)
  essential_in_universe <- sum(all_perturb_gene_symbols_in_analysis %in% essential_gene_list)
  non_essential_in_universe <- total_genes_in_universe - essential_in_universe
  
  message(paste0("  Analysis Universe: Total unique perturbation genes in current combined results: ", total_genes_in_universe, "."))
  message(paste0("  Analysis Universe: Essential genes in this universe: ", essential_in_universe, " (", round(essential_in_universe/total_genes_in_universe*100, 2), "%)", "; Non-essential: ", non_essential_in_universe, "."))
  
  if (essential_in_universe == 0 || non_essential_in_universe == 0) {
    message(paste0("  Warning: Only essential (", essential_in_universe, ") or only non-essential (", non_essential_in_universe, ") genes in the analysis universe. Fisher's exact test might be trivial or fail if categories are also zero. Proceeding anyway."))
  }
  
  results <- list()
  
  # --- Gero-advancers (Positive combined_NES, significant combined_padj) ---
  sig_geroadvancers <- combined_fgsea_df %>%
    dplyr::filter(combined_padj < 0.05, combined_NES > 0)
  
  if (nrow(sig_geroadvancers) > 0) {
    geroadvancer_genes <- unique(gsub("^(.*?)\\s+Knockdown Signature - .*", "\\1", sig_geroadvancers$pathway))
    
    in_category_essential <- sum(geroadvancer_genes %in% essential_gene_list)
    in_category_non_essential <- length(geroadvancer_genes) - in_category_essential
    
    # These counts are for genes in the *total analysis universe* but *not* in the current category
    not_in_category_essential <- essential_in_universe - in_category_essential
    not_in_category_non_essential <- non_essential_in_universe - in_category_non_essential
    
    # Ensure no negative counts in case of small sample sizes or edge cases
    not_in_category_essential <- max(0, not_in_category_essential)
    not_in_category_non_essential <- max(0, not_in_category_non_essential)
    
    contingency_table_advancer <- matrix(c(in_category_essential, not_in_category_essential,
                                           in_category_non_essential, not_in_category_non_essential),
                                         nrow = 2, byrow = TRUE,
                                         dimnames = list(Essentiality = c("Essential", "Non-Essential"),
                                                         Category = c("Gero-advancer", "Not Gero-advancer")))
    
    message(paste0("\n  Gero-advancers (Positive Combined NES) in ", analysis_name, ":"))
    message("    Contingency Table:")
    print(contingency_table_advancer)
    
    # Check if the table is valid for Fisher's Exact Test (no zero margin sums)
    if (min(rowSums(contingency_table_advancer)) > 0 && min(colSums(contingency_table_advancer)) > 0) {
      fisher_res_advancer <- fisher.test(contingency_table_advancer)
      message(paste0("    Fisher's Exact Test p-value: ", format(fisher_res_advancer$p.value, digits = 3)))
      message(paste0("    Odds Ratio (Essential in Gero-advancer vs. Not Gero-advancer): ", format(fisher_res_advancer$estimate, digits = 3)))
      results$geroadvancers <- list(table = contingency_table_advancer, test = fisher_res_advancer)
    } else {
      message("    Fisher's Exact Test skipped due to insufficient data (zero margin sums in contingency table).")
      results$geroadvancers <- list(table = contingency_table_advancer, test = "Skipped (Insufficient data)")
    }
  } else {
    message(paste0("\n  No significant Gero-advancer genes found for '", analysis_name, "' (combined_padj < 0.05, combined_NES > 0)."))
    results$geroadvancers <- NULL
  }
  
  # --- Gero-protectors (Negative combined_NES, significant combined_padj) ---
  sig_geroprotectors <- combined_fgsea_df %>%
    dplyr::filter(combined_padj < 0.05, combined_NES < 0)
  
  if (nrow(sig_geroprotectors) > 0) {
    geroprotector_genes <- unique(gsub("^(.*?)\\s+Knockdown Signature - .*", "\\1", sig_geroprotectors$pathway))
    
    in_category_essential <- sum(geroprotector_genes %in% essential_gene_list)
    in_category_non_essential <- length(geroprotector_genes) - in_category_essential
    
    not_in_category_essential <- essential_in_universe - in_category_essential
    not_in_category_non_essential <- non_essential_in_universe - in_category_non_essential
    
    # Ensure no negative counts
    not_in_category_essential <- max(0, not_in_category_essential)
    not_in_category_non_essential <- max(0, not_in_category_non_essential)
    
    contingency_table_protector <- matrix(c(in_category_essential, not_in_category_essential,
                                            in_category_non_essential, not_in_category_non_essential),
                                          nrow = 2, byrow = TRUE,
                                          dimnames = list(Essentiality = c("Essential", "Non-Essential"),
                                                          Category = c("Gero-protector", "Not Gero-protector")))
    
    message(paste0("\n  Gero-protectors (Negative Combined NES) in ", analysis_name, ":"))
    message("    Contingency Table:")
    print(contingency_table_protector)
    
    # Check if the table is valid for Fisher's Exact Test (no zero margin sums)
    if (min(rowSums(contingency_table_protector)) > 0 && min(colSums(contingency_table_protector)) > 0) {
      fisher_res_protector <- fisher.test(contingency_table_protector)
      message(paste0("    Fisher's Exact Test p-value: ", format(fisher_res_protector$p.value, digits = 3)))
      message(paste0("    Odds Ratio (Essential in Gero-protector vs. Not Gero-protector): ", format(fisher_res_protector$estimate, digits = 3)))
      results$geroprotectors <- list(table = contingency_table_protector, test = fisher_res_protector)
    } else {
      message("    Fisher's Exact Test skipped due to insufficient data (zero margin sums in contingency table).")
      results$geroprotectors <- list(table = contingency_table_protector, test = "Skipped (Insufficient data)")
    }
  } else {
    message(paste0("\n  No significant Gero-protector genes found for '", analysis_name, "' (combined_padj < 0.05, combined_NES < 0)."))
    results$geroprotectors <- NULL
  }
  
  message(paste0("\n--- Essential Gene Enrichment Test for '", analysis_name, "' Complete ---"))
  return(results)
}






# --- Main Script Execution ---

# 1. Load OmicSignature Collections
perturb_collection <- readRDS(perturbation_collection_file)
message(paste0("Loaded ", length(perturb_collection$OmicSigList), " perturbation signatures."))

gtex_collection <- readRDS(gtex_collection_file)
message(paste0("Loaded ", length(gtex_collection$OmicSigList), " GTEX aging signatures."))

gene_map <- readRDS(gene_map_file)
message(paste0("Loaded gene map from: ", gene_map_file))

# Pre-process gene_map once
message("Pre-processing gene_map for faster lookups (applying simplify_entry once)...")
gene_map_simplified <- sapply(gene_map, simplify_entry, USE.NAMES = FALSE)
names(gene_map_simplified) <- names(gene_map) 
message("Gene map pre-processing complete.")



# 2. Age-Centered Analysis
message("\n--- Running Age-Centered Analysis ---")
message("Goal: Compare perturbation gene sets against aging ranked lists.")

# Prepare ranked lists from GTEX aging data
ranked_lists_age <- extract_ranked_lists(gtex_collection, "Aging", gene_map_simplified) 

# Prepare UP/DN gene sets from perturbation data
perturb_gene_sets_up_dn <- extract_gene_sets_up_dn(perturb_collection, logFC_filter_val, adj_pval_filter_val, geneset_top_n, "Perturbation", gene_map_simplified) 

# Add messages for scale diagnostics
message(paste0("  Age-Centered: Number of Aging ranked lists: ", length(ranked_lists_age)))
message(paste0("  Age-Centered: Number of Perturbation UP gene sets: ", length(perturb_gene_sets_up_dn$up_gene_sets)))
message(paste0("  Age-Centered: Number of Perturbation DN gene sets: ", length(perturb_gene_sets_up_dn$dn_gene_sets)))

# Perform fgsea and combine results
fgsea_res_age_centered <- perform_fgsea_and_combine(
  ranked_lists = ranked_lists_age,
  gene_sets_up = perturb_gene_sets_up_dn$up_gene_sets,
  gene_sets_dn = perturb_gene_sets_up_dn$dn_gene_sets,
  ranked_source_name = "Aging",
  geneset_source_name = "Perturbation",
  total_num_cores = total_num_cores,
  fgsea_min_size = fgsea_min_size,
  fgsea_max_size = fgsea_max_size
)

# Save the full results table
if (!is.null(fgsea_res_age_centered)) {
  data.table::fwrite(fgsea_res_age_centered, file.path(output_dir, "fgsea_results_age_centered.csv"))
  message("Full Age-Centered fgsea results saved to: ", file.path(output_dir, "fgsea_results_age_centered.csv"))
}

# Visualize results
plot_fgsea_results(fgsea_res_age_centered, "Age-Centered Analysis", output_dir, top_n_genesets_to_plot)
plot_fgsea_results(fgsea_res_age_centered, "Age-Centered Analysis", output_dir, 100) # Plot with bigger number of top genes

fgsea_res_age_centered_up_gs <- fgsea_res_age_centered %>% filter(geneset_direction == "UP")
plot_clustered_heatmap(fgsea_res_age_centered_up_gs, "Age-Centered Analysis", "Perturbation", "UP", output_dir)

fgsea_res_age_centered_dn_gs <- fgsea_res_age_centered %>% filter(geneset_direction == "DN")
plot_clustered_heatmap(fgsea_res_age_centered_dn_gs, "Age-Centered Analysis", "Perturbation", "DN", output_dir)


if (!is.null(fgsea_res_age_centered)) {
  # Calculate combined scores
  fgsea_res_age_centered_combined <- calculate_combined_scores(fgsea_res_age_centered, "Age-Centered Analysis")
  
  if (!is.null(fgsea_res_age_centered_combined) && nrow(fgsea_res_age_centered_combined) > 0) {
    message("\n--- Generating NEW Age-Centered Visualizations ---")
    plot_combined_fgsea_dotplot(fgsea_res_age_centered_combined, "Age-Centered Analysis", output_dir, top_n_combined_plots, essential_gene_list)
    plot_combined_fgsea_dotplot(fgsea_res_age_centered_combined, "Age-Centered Analysis", output_dir, 100, essential_gene_list)
    
    # Store heatmap clustering results to pass to the p-value heatmap
    nes_heatmap_clusters_age <- plot_combined_heatmap(
      combined_fgsea_df = fgsea_res_age_centered_combined,
      analysis_title_prefix = "Age-Centered Analysis",
      output_dir = output_dir,
      essential_gene_list = essential_gene_list,
      ranked_source_name = "Aging" # Specify that ranked lists are from Aging, so perturbation genes are in rows
    )

    # Pass clustering results to the signed p-value heatmap
    plot_combined_heatmap_signed_pvalue_NES(
      combined_fgsea_df = fgsea_res_age_centered_combined,
      analysis_title_prefix = "Age-Centered Analysis",
      output_dir = output_dir,
      essential_gene_list = essential_gene_list,
      ranked_source_name = "Aging", # Specify that ranked lists are from Aging, so perturbation genes are in rows
      row_dend = nes_heatmap_clusters_age$row_dend,
      col_dend = nes_heatmap_clusters_age$col_dend
    )
    
    # New volcano plot with coloring
    plot_volcano_with_coloring(fgsea_res_age_centered_combined, "Age-Centered Analysis", output_dir)
    
    # Exploratory plots for raw fgsea results
    plot_pval_vs_es(fgsea_res_age_centered, "Age-Centered Analysis (Raw FGSEA)", output_dir, es_col = "NES", pval_col = "padj")
    plot_geneset_size_vs_es(fgsea_res_age_centered, "Age-Centered Analysis (Raw FGSEA)", output_dir, es_col = "NES")
    # Exploratory plots for combined scores
    plot_pval_vs_es(fgsea_res_age_centered_combined, "Age-Centered Analysis (Combined NES)", output_dir, es_col = "combined_NES", pval_col = "combined_padj")
    plot_geneset_size_vs_es(fgsea_res_age_centered_combined, "Age-Centered Analysis (Combined NES)", output_dir, es_col = "combined_NES")
    
    # Essential gene enrichment test for Age-Centered Analysis
    essential_enrichment_age <- perform_essential_gene_enrichment_test(
      combined_fgsea_df = fgsea_res_age_centered_combined,
      essential_gene_list = essential_gene_list,
      analysis_name = "Age-Centered Analysis",
      nr_all_depmap_genes = nr_all_depmap_genes,
      nr_ess_depmap_genes = nr_ess_depmap_genes
    )
    
    # Heatmaps filtered for non-essential genes in Age-Centered Analysis
    message("\n--- Generating Age-Centered Heatmaps (Only Non-Essential Perturbation Genes) ---")
    nes_heatmap_clusters_age_filtered <- plot_combined_heatmap(
      combined_fgsea_df = fgsea_res_age_centered_combined,
      analysis_title_prefix = "Age-Centered Analysis only non-essential genes",
      output_dir = output_dir,
      essential_gene_list = essential_gene_list,
      filter = TRUE,
      ranked_source_name = "Aging" # Perturbation genes are in rows for Age-Centered
    )
    plot_combined_heatmap_signed_pvalue_NES(
      combined_fgsea_df = fgsea_res_age_centered_combined,
      analysis_title_prefix = "Age-Centered Analysis only non-essential genes",
      output_dir = output_dir,
      essential_gene_list = essential_gene_list,
      filter = TRUE,
      ranked_source_name = "Aging", # Perturbation genes are in rows for Age-Centered
      row_dend = nes_heatmap_clusters_age_filtered$row_dend,
      col_dend = nes_heatmap_clusters_age_filtered$col_dend
    )
    
  } else {
    message("Skipping new Age-Centered visualizations due to no combined results.")
  }
}



# 3. Perturbation-Centered Analysis
message("\n--- Running Perturbation-Centered Analysis ---")
message("Goal: Compare aging gene sets against perturbation ranked lists.")

# Prepare ranked lists from perturbation data
ranked_lists_perturb <- extract_ranked_lists(perturb_collection, "Perturbation", gene_map_simplified) 

# Prepare UP/DN gene sets from GTEX aging data
age_gene_sets_up_dn <- extract_gene_sets_up_dn(gtex_collection, logFC_filter_val, adj_pval_filter_val, geneset_top_n, "Aging", gene_map_simplified) 

# Add messages for scale diagnostics
message(paste0("  Perturbation-Centered: Number of Perturbation ranked lists: ", length(ranked_lists_perturb)))
message(paste0("  Perturbation-Centered: Number of Aging UP gene sets: ", length(age_gene_sets_up_dn$up_gene_sets)))
message(paste0("  Perturbation-Centered: Number of Aging DN gene sets: ", length(age_gene_sets_up_dn$dn_gene_sets)))

# Perform fgsea and combine results
fgsea_res_perturb_centered <- perform_fgsea_and_combine(
  ranked_lists = ranked_lists_perturb,
  gene_sets_up = age_gene_sets_up_dn$up_gene_sets,
  gene_sets_dn = age_gene_sets_up_dn$dn_gene_sets,
  ranked_source_name = "Perturbation",
  geneset_source_name = "Aging",
  total_num_cores = total_num_cores,
  fgsea_min_size = fgsea_min_size,
  fgsea_max_size = fgsea_max_size
)

# Save the full results table
if (!is.null(fgsea_res_perturb_centered)) {
  data.table::fwrite(fgsea_res_perturb_centered, file.path(output_dir, "fgsea_results_perturbation_centered.csv"))
  message("Full Perturbation-Centered fgsea results saved to: ", file.path(output_dir, "fgsea_results_perturbation_centered.csv"))
}

# Visualize results
plot_fgsea_results(fgsea_res_perturb_centered, "Perturbation-Centered Analysis", output_dir, top_n_genesets_to_plot)
plot_fgsea_results(fgsea_res_perturb_centered, "Perturbation-Centered Analysis", output_dir, 100)

fgsea_res_perturb_centered_up_gs <- fgsea_res_perturb_centered %>% filter(geneset_direction == "UP")
plot_clustered_heatmap(fgsea_res_perturb_centered_up_gs, "Perturbation-Centered Analysis", "Aging", "UP", output_dir)

fgsea_res_perturb_centered_dn_gs <- fgsea_res_perturb_centered %>% filter(geneset_direction == "DN")
plot_clustered_heatmap(fgsea_res_perturb_centered_dn_gs, "Perturbation-Centered Analysis", "Aging", "DN", output_dir)


if (!is.null(fgsea_res_perturb_centered)) {
  # Calculate combined scores
  fgsea_res_perturb_centered_combined <- calculate_combined_scores(fgsea_res_perturb_centered, "Perturbation-Centered Analysis")
  
  if (!is.null(fgsea_res_perturb_centered_combined) && nrow(fgsea_res_perturb_centered_combined) > 0) {
    message("\n--- Generating NEW Perturbation-Centered Visualizations ---")
    plot_combined_fgsea_dotplot(fgsea_res_perturb_centered_combined, "Perturbation-Centered Analysis", output_dir, top_n_combined_plots, essential_gene_list)
    plot_combined_fgsea_dotplot(fgsea_res_perturb_centered_combined, "Perturbation-Centered Analysis", output_dir, 100, essential_gene_list)
    
    # Store heatmap clustering results to pass to the p-value heatmap
    nes_heatmap_clusters_perturb <- plot_combined_heatmap(
      combined_fgsea_df = fgsea_res_perturb_centered_combined,
      analysis_title_prefix = "Perturbation-Centered Analysis",
      output_dir = output_dir,
      essential_gene_list = essential_gene_list,
      ranked_source_name = "Perturbation" # Specify that ranked lists are from Perturbation, so perturbation genes are in columns
    )

    # Pass clustering results to the signed p-value heatmap
    plot_combined_heatmap_signed_pvalue_NES(
      combined_fgsea_df = fgsea_res_perturb_centered_combined,
      analysis_title_prefix = "Perturbation-Centered Analysis",
      output_dir = output_dir,
      essential_gene_list = essential_gene_list,
      ranked_source_name = "Perturbation", # Specify that ranked lists are from Perturbation, so perturbation genes are in columns
      row_dend = nes_heatmap_clusters_perturb$row_dend,
      col_dend = nes_heatmap_clusters_perturb$col_dend
    )
    
    # New volcano plot with coloring
    plot_volcano_with_coloring(fgsea_res_perturb_centered_combined, "Perturbation-Centered Analysis", output_dir)
    
    # Exploratory plots for raw fgsea results
    plot_pval_vs_es(fgsea_res_perturb_centered, "Perturbation-Centered Analysis (Raw FGSEA)", output_dir, es_col = "NES", pval_col = "padj")
    plot_geneset_size_vs_es(fgsea_res_perturb_centered, "Perturbation-Centered Analysis (Raw FGSEA)", output_dir, es_col = "NES")
    # Exploratory plots for combined scores
    plot_pval_vs_es(fgsea_res_perturb_centered_combined, "Perturbation-Centered Analysis (Combined NES)", output_dir, es_col = "combined_NES", pval_col = "combined_padj")
    plot_geneset_size_vs_es(fgsea_res_perturb_centered_combined, "Perturbation-Centered Analysis (Combined NES)", output_dir, es_col = "combined_NES")
    
    # Essential gene enrichment test for Perturbation-Centered Analysis
    essential_enrichment_perturb <- perform_essential_gene_enrichment_test(
      combined_fgsea_df = fgsea_res_perturb_centered_combined,
      essential_gene_list = essential_gene_list,
      analysis_name = "Perturbation-Centered Analysis",
      nr_all_depmap_genes = nr_all_depmap_genes,
      nr_ess_depmap_genes = nr_ess_depmap_genes
    )
    
    # filtered for non-essential genes - also pass clustering
    message("\n--- Generating Perturbation-Centered Heatmaps (Only Non-Essential Perturbation Genes) ---")
    nes_heatmap_clusters_perturb_filtered <- plot_combined_heatmap(
      combined_fgsea_df = fgsea_res_perturb_centered_combined,
      analysis_title_prefix = "Perturbation-Centered Analysis only non-essential genes",
      output_dir = output_dir,
      essential_gene_list = essential_gene_list,
      filter = TRUE,
      ranked_source_name = "Perturbation" # Specify that ranked lists are from Perturbation, so perturbation genes are in columns
    )
    plot_combined_heatmap_signed_pvalue_NES(
      combined_fgsea_df = fgsea_res_perturb_centered_combined,
      analysis_title_prefix = "Perturbation-Centered Analysis only non-essential genes",
      output_dir = output_dir,
      essential_gene_list = essential_gene_list,
      filter = TRUE,
      ranked_source_name = "Perturbation", # Specify that ranked lists are from Perturbation, so perturbation genes are in columns
      row_dend = nes_heatmap_clusters_perturb_filtered$row_dend,
      col_dend = nes_heatmap_clusters_perturb_filtered$col_dend
    )
    
    
  } else {
    message("Skipping new Perturbation-Centered visualizations due to no combined results.")
  }
}



# Cross-Analysis Scatter Plot (Pass essential gene list, now plots NES only)
if (!is.null(fgsea_res_age_centered_combined) && !is.null(fgsea_res_perturb_centered_combined)) {
  plot_cross_analysis_scatterplots(fgsea_res_age_centered_combined, fgsea_res_perturb_centered_combined, output_dir, essential_gene_list = essential_gene_list)
  plot_cross_analysis_scatterplots(fgsea_res_age_centered_combined, fgsea_res_perturb_centered_combined, output_dir, essential_gene_list = essential_gene_list, filter = FALSE)
} else {
  message("\nSkipping Cross-Analysis scatter plots due to missing combined results from one or both analyses.")
}


message("\n--- GSEA Combined Analysis Complete ---")

