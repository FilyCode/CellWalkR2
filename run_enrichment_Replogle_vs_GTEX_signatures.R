library(tidyverse)
library(fgsea)
library(OmicSignature)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# --- Configuration ---
# File paths
work_dir <- "/restricted/projectnb/agedisease/projects/challenge2025/"
perturbation_collection_file <- paste0(work_dir, "results/perturbational_omic_sigs/replogle_2022/Replogle_Perturb_Combined_OmicSignatureCollection.rds")
gtex_collection_file <- paste0(work_dir, "data/GTEX/GTEX_aging_omic_col_stat_v114.rds")
output_dir <- paste0(work_dir, "results/fgsea_combined_analysis/Replogle_GTEX")

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Filtering parameters for gene sets
logFC_filter_val <- 0.25
adj_pval_filter_val <- 0.05
geneset_top_n <- 500
top_n_genesets_to_plot <- 20

# fgsea parameters
fgsea_min_size <- 1
fgsea_max_size <- 500

message("--- Starting GSEA Combined Analysis ---")
message(paste0("Output directory: ", output_dir))
message(paste0("Gene set filtering: |logFC| > ", logFC_filter_val, ", adj.pval < ", adj_pval_filter_val))
message(paste0("fgsea parameters: minSize = ", fgsea_min_size, ", maxSize = ", fgsea_max_size))
message(paste0("Top N results for plotting: ", top_n_genesets_to_plot))

# --- Helper Functions ---

#' Extracts ranked gene lists from an OmicSignatureCollection.
#' Ranks genes by 'score' without p-value/logFC filtering.
#' Prioritizes 'difexp' table, falls back to 'signature' if 'difexp' is missing.
#'
#' @param omic_collection An OmicSignatureCollection object.
#' @param collection_name A string identifying the source collection (for warnings/messages).
#' @return A named list of numeric vectors, where names are probe_id and values are scores.
extract_ranked_lists <- function(omic_collection, collection_name) {
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
    
    current_ranks <- df_for_ranks %>%
      dplyr::select(probe_id, score) %>%
      drop_na(score) %>% # Remove NAs in score
      distinct(probe_id, .keep_all = TRUE) %>% # Keep one entry per probe_id if duplicates exist
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
#'
#' @param omic_collection An OmicSignatureCollection object.
#' @param logFC_thresh Numeric, absolute logFC threshold for filtering.
#' @param pval_thresh Numeric, adjusted p-value threshold for filtering.
#' @param geneset_top_n Numeric or NULL. If a number, takes top N genes by absolute logFC after other filters.
#' @param collection_name A string identifying the source collection (for warnings/messages).
#' @return A list containing two named lists: 'up_gene_sets' and 'dn_gene_sets'.
extract_gene_sets_up_dn <- function(omic_collection, logFC_thresh, pval_thresh, geneset_top_n, collection_name) {
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
  
  for (sig_name in names(omic_collection$OmicSigList)) {
    sig_obj <- omic_collection$OmicSigList[[sig_name]]
    
    # For debugging, uncomment to see column names for each signature
    # message(paste0("  Processing signature: ", sig_name))
    # if (!is.null(sig_obj$signature)) { print(paste("    signature cols:", paste(names(sig_obj$signature), collapse=", "))) }
    # if (!is.null(sig_obj$difexp)) { print(paste("    difexp cols:", paste(names(sig_obj$difexp), collapse=", "))) }
    
    
    df_to_process <- NULL
    used_source_name <- "none" 
    
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
    
    # Filter for NAs in critical columns and ensure unique gene IDs
    # Filter NAs for mandatory columns (score, adj.pval)
    df_filtered <- df_to_process %>% 
      dplyr::filter(!is.na(!!sym(required_score_col)), !is.na(adj.pval)) %>% 
      dplyr::distinct(!!sym(required_id_col), .keep_all = TRUE) 
    
    # If logFC is present, filter its NAs too
    if (has_logFC) {
      df_filtered <- df_filtered %>% dplyr::filter(!is.na(!!sym(possible_logfc_cols[1])))
    }
    
    
    # Determine if logFC is available and if group_label is a valid fallback for directionality
    has_logFC <- any(possible_logfc_cols %in% names(df_to_process)) && 
      !all(is.na(df_to_process[[possible_logfc_cols[1]]])) # logFC column exists and is not all NA
    
    use_group_label_for_direction <- FALSE
    group_label_col_name <- NULL
    
    if (!has_logFC) { # If logFC is not available or all NA
      if (any(possible_group_label_cols %in% names(df_to_process))) {
        group_label_col_name <- intersect(possible_group_label_cols, names(df_to_process))[1]
        # Check if group_label has the expected values and is not all NA
        group_labels_in_df <- unique(df_to_process[[group_label_col_name]])
        if (all(c("Older", "Younger") %in% group_labels_in_df) && !all(is.na(df_to_process[[group_label_col_name]]))) {
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
    
    # Filter for upregulated genes and downregulated genes (dynamic directionality)
    up_genes <- character(0)
    dn_genes <- character(0)
    
    if (has_logFC) { # Use logFC for direction
      up_genes <- significant_genes %>% dplyr::filter(!!sym(possible_logfc_cols[1]) > 0) %>% dplyr::pull(!!sym(required_id_col)) %>% unique()
      dn_genes <- significant_genes %>% dplyr::filter(!!sym(possible_logfc_cols[1]) < 0) %>% dplyr::pull(!!sym(required_id_col)) %>% unique()
    } else if (use_group_label_for_direction) { # Use group_label for direction
      up_genes <- significant_genes %>% dplyr::filter(!!sym(group_label_col_name) == "Older") %>% dplyr::pull(!!sym(required_id_col)) %>% unique()
      dn_genes <- significant_genes %>% dplyr::filter(!!sym(group_label_col_name) == "Younger") %>% dplyr::pull(!!sym(required_id_col)) %>% unique()
    }
    
    # Add to gene sets list
    if (length(up_genes) > 0) {
      up_gene_sets[[paste0(sig_name, "_UP")]] <- up_genes
    }
    if (length(dn_genes) > 0) {
      dn_gene_sets[[paste0(sig_name, "_DN")]] <- dn_genes
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
#' @return A combined data frame of fgsea results.
perform_fgsea_and_combine <- function(ranked_lists, gene_sets_up, gene_sets_dn, ranked_source_name, geneset_source_name) {
  all_fgsea_results <- list()
  
  message(paste0("Initiating fgsea analysis for '", ranked_source_name, "' ranked lists vs. '", geneset_source_name, "' gene sets."))
  
  # Run for UP gene sets
  if (length(gene_sets_up) > 0) {
    message(paste0("  Processing with UP gene sets from '", geneset_source_name, "'..."))
    for (rl_name in names(ranked_lists)) {
      message(paste0("    Analyzing ranked list: ", rl_name, "..."))
      fg_up <- fgsea(pathways = gene_sets_up,
                     stats    = ranked_lists[[rl_name]],
                     minSize  = fgsea_min_size,
                     maxSize  = fgsea_max_size)
      fg_up$ranked_list_name <- rl_name
      fg_up$geneset_source_name <- geneset_source_name
      fg_up$geneset_direction <- "UP"
      all_fgsea_results[[paste0(rl_name, "_UP_gene_sets")]] <- fg_up
    }
  } else {
    message(paste("No UP gene sets available from '", geneset_source_name, "' for fgsea."))
  }
  
  # Run for DN gene sets
  if (length(gene_sets_dn) > 0) {
    message(paste0("  Processing with DN gene sets from '", geneset_source_name, "'..."))
    for (rl_name in names(ranked_lists)) {
      message(paste0("    Analyzing ranked list: ", rl_name, "..."))
      fg_dn <- fgsea(pathways = gene_sets_dn,
                     stats    = ranked_lists[[rl_name]],
                     minSize  = fgsea_min_size,
                     maxSize  = fgsea_max_size)
      fg_dn$ranked_list_name <- rl_name
      fg_dn$geneset_source_name <- geneset_source_name
      fg_dn$geneset_direction <- "DN"
      all_fgsea_results[[paste0(rl_name, "_DN_gene_sets")]] <- fg_dn
    }
  } else {
    message(paste("No DN gene sets available from '", geneset_source_name, "' for fgsea."))
  }
  
  combined_results_df <- do.call(rbind, all_fgsea_results)
  if (is.null(combined_results_df) || nrow(combined_results_df) == 0) {
    warning(paste("No fgsea results generated for ", ranked_source_name, " vs ", geneset_source_name, " analysis type."))
    return(NULL)
  }
  message(paste0("fgsea analysis complete for '", ranked_source_name, "' vs. '", geneset_source_name, "'. Total results: ", nrow(combined_results_df), " rows."))
  return(combined_results_df)
}

#' Generates and saves ggplot visualizations for fgsea results.
#'
#' @param fgsea_df A data frame of combined fgsea results.
#' @param analysis_title_prefix A string for plot titles, e.g., "Age-Centered Analysis".
#' @param output_dir Path to save the plots.
#' @param top_n Numeric, number of top results to plot.
plot_fgsea_results <- function(fgsea_df, analysis_title_prefix, output_dir, top_n = 20) {
  
  if (is.null(fgsea_df) || nrow(fgsea_df) == 0) {
    message(paste("No fgsea results for", analysis_title_prefix, "to display or plot."))
    return(invisible(NULL))
  }
  
  message(paste0("Generating plots for: ", analysis_title_prefix))
  
  # First, filter for all significant results
  all_significant_results <- fgsea_df %>%
    filter(padj < 0.05)
  
  if (nrow(all_significant_results) == 0) {
    message(paste("No significant interactions (padj < 0.05) found for", analysis_title_prefix, "to plot."))
    return(invisible(NULL))
  } else {
    message(paste0("Found ", nrow(all_significant_results), " significant interactions (padj < 0.05) for ", analysis_title_prefix, "."))
    
    # --- Plot for Top Upregulated Signatures (Positive NES) ---
    top_upregulated_results <- all_significant_results %>%
      filter(NES > 0) %>% # Select positive NES
      arrange(dplyr::desc(NES)) %>% # Order by highest NES
      slice(1:min(dplyr::n(), top_n)) # Take top N, or fewer if less than N exist
    
    if (nrow(top_upregulated_results) > 0) {
      message(paste0("  Generating plot for top ", top_n, " upregulated results (Positive NES)..."))
      plot_upregulated <- ggplot(top_upregulated_results, aes(x = reorder(pathway, NES), y = ranked_list_name)) +
        geom_point(aes(size  = size, color = NES)) +
        scale_size_continuous(name = "Gene Set Size") +
        scale_color_gradient(low = "yellow", high = "red", name = "NES") +
        coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(face = "bold", hjust = 0.5),
              axis.title.y = element_text(face = "bold")) +
        labs(x = "Gene Set (Pathway)",
             y = paste0(analysis_title_prefix, " Ranked List"),
             title = paste0("Top ", top_n, " Upregulated Gene Sets in ", analysis_title_prefix, " (Positive NES)"))
      
      plot_filename <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_top_", top_n, "_upregulated_enrichment.pdf"))
      ggsave(plot_filename, plot_upregulated, width = 12, height = 8)
      message("  Plot generated and saved to: ", plot_filename)
      
      message(paste0("\nTop ", top_n, " Upregulated Results (pathway, ranked_list_name, NES, pval, padj, size):"))
      print(top_upregulated_results %>% dplyr::select(pathway, ranked_list_name, NES, pval, padj, size))
    } else {
      message(paste("  No significant upregulated interactions (NES > 0, padj < 0.05) found for", analysis_title_prefix, "to plot."))
    }
    
    # --- Plot for Top Downregulated Signatures (Negative NES) ---
    top_downregulated_results <- all_significant_results %>%
      filter(NES < 0) %>% # Select negative NES
      arrange(NES) %>% # Order by lowest (most negative) NES
      slice(1:min(dplyr::n(), top_n)) # Take top N, or fewer if less than N exist
    
    if (nrow(top_downregulated_results) > 0) {
      message(paste0("  Generating plot for top ", top_n, " downregulated results (Negative NES)..."))
      plot_downregulated <- ggplot(top_downregulated_results, aes(x = reorder(pathway, NES), y = ranked_list_name)) +
        geom_point(aes(size  = size, color = NES)) +
        scale_size_continuous(name = "Gene Set Size") +
        scale_color_gradient(low = "darkblue", high = "lightblue", name = "NES") +
        coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(face = "bold", hjust = 0.5),
              axis.title.y = element_text(face = "bold")) +
        labs(x = "Gene Set (Pathway)",
             y = paste0(analysis_title_prefix, " Ranked List"),
             title = paste0("Top ", top_n, " Downregulated Gene Sets in ", analysis_title_prefix, " (Negative NES)"))
      
      plot_filename <- file.path(output_dir, paste0(gsub(" ", "_", analysis_title_prefix), "_top_", top_n, "_downregulated_enrichment.pdf"))
      ggsave(plot_filename, plot_downregulated, width = 12, height = 8)
      message("  Plot generated and saved to: ", plot_filename)
      
      message(paste0("\nTop ", top_n, " Downregulated Results (pathway, ranked_list_name, NES, pval, padj, size):"))
      print(top_downregulated_results %>% dplyr::select(pathway, ranked_list_name, NES, pval, padj, size))
    } else {
      message(paste("  No significant downregulated interactions (NES < 0, padj < 0.05) found for", analysis_title_prefix, "to plot."))
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
  
  # Reshape data for heatmap: pathways as rows, ranked lists as columns, NES as values
  heatmap_data <- significant_results %>%
    dplyr::select(pathway, ranked_list_name, NES) %>%
    tidyr::pivot_wider(names_from = ranked_list_name, values_from = NES, values_fill = 0) # Fill non-significant/missing with 0
  
  # Convert to matrix, setting row names
  mat <- as.matrix(heatmap_data %>% dplyr::select(-pathway))
  rownames(mat) <- heatmap_data$pathway
  
  # Determine a symmetric color range based on max absolute NES
  max_abs_nes <- max(abs(mat), na.rm = TRUE)
  # Define the color function using circlize::colorRamp2 for a diverging palette
  col_fun <- colorRamp2(c(-max_abs_nes, -max_abs_nes/2, 0, max_abs_nes/2, max_abs_nes), 
                        c("darkblue", "lightblue", "white", "pink2", "darkred")) # Matches previous pheatmap color scheme
  
  # Construct the plot title and filename
  hm_title <- paste0("Clustered Heatmap: ", analysis_title_prefix, "\n(", geneset_source_name, " ", geneset_direction, " Gene Sets vs. Ranked Lists)")
  file_name_base <- paste0(gsub(" ", "_", analysis_title_prefix), "_", geneset_source_name, "_", geneset_direction, "_heatmap")
  heatmap_output_path <- file.path(output_dir, paste0(file_name_base, ".png"))
  
  # Generate heatmap using ComplexHeatmap
  hm <- Heatmap(
    mat,
    name = "NES", # Legend name
    col = col_fun,
    na_col = "grey90", # Color for NA values
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE, # Remove row names
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 90
  )
  
  # Save heatmap to PNG using png() and draw()
  message("  Saving heatmap to: ", heatmap_output_path)
  png(heatmap_output_path, width = 2200, height = 1800, res = 300)
  draw(hm, column_title = hm_title)
  dev.off()
  
  message("  Clustered heatmap generated and saved to: ", heatmap_filename)
}



# --- Main Script Execution ---

# 1. Load OmicSignature Collections
perturb_collection <- readRDS(perturbation_collection_file)
message(paste0("Loaded ", length(perturb_collection$OmicSigList), " perturbation signatures."))

gtex_collection <- readRDS(gtex_collection_file)
message(paste0("Loaded ", length(gtex_collection$OmicSigList), " GTEX aging signatures."))


# 2. Age-Centered Analysis
message("\n--- Running Age-Centered Analysis ---")
message("Goal: Compare perturbation gene sets against aging ranked lists.")

# Prepare ranked lists from GTEX aging data
ranked_lists_age <- extract_ranked_lists(gtex_collection, "Aging")

# Prepare UP/DN gene sets from perturbation data
perturb_gene_sets_up_dn <- extract_gene_sets_up_dn(perturb_collection, logFC_filter_val, adj_pval_filter_val, geneset_top_n, "Perturbation")

# Perform fgsea and combine results
fgsea_res_age_centered <- perform_fgsea_and_combine(
  ranked_lists = ranked_lists_age,
  gene_sets_up = perturb_gene_sets_up_dn$up_gene_sets,
  gene_sets_dn = perturb_gene_sets_up_dn$dn_gene_sets,
  ranked_source_name = "Aging",
  geneset_source_name = "Perturbation"
)

# Visualize results
plot_fgsea_results(fgsea_res_age_centered, "Age-Centered Analysis", output_dir, top_n_genesets_to_plot)

fgsea_res_age_centered_up_gs <- fgsea_res_age_centered %>% filter(geneset_direction == "UP")
plot_clustered_heatmap(fgsea_res_age_centered_up_gs, "Age-Centered Analysis", "Perturbation", "UP", output_dir)

fgsea_res_age_centered_dn_gs <- fgsea_res_age_centered %>% filter(geneset_direction == "DN")
plot_clustered_heatmap(fgsea_res_age_centered_dn_gs, "Age-Centered Analysis", "Perturbation", "DN", output_dir)


# Save the full results table
if (!is.null(fgsea_res_age_centered)) {
  write_tsv(fgsea_res_age_centered, file.path(output_dir, "fgsea_results_age_centered.tsv"))
  message("Full Age-Centered fgsea results saved to: ", file.path(output_dir, "fgsea_results_age_centered.tsv"))
}


# 3. Perturbation-Centered Analysis
message("\n--- Running Perturbation-Centered Analysis ---")
message("Goal: Compare aging gene sets against perturbation ranked lists.")

# Prepare ranked lists from perturbation data
ranked_lists_perturb <- extract_ranked_lists(perturb_collection, "Perturbation")

# Prepare UP/DN gene sets from GTEX aging data
age_gene_sets_up_dn <- extract_gene_sets_up_dn(gtex_collection, logFC_filter_val, adj_pval_filter_val, geneset_top_n, "Aging")

# Perform fgsea and combine results
fgsea_res_perturb_centered <- perform_fgsea_and_combine(
  ranked_lists = ranked_lists_perturb,
  gene_sets_up = age_gene_sets_up_dn$up_gene_sets,
  gene_sets_dn = age_gene_sets_up_dn$dn_gene_sets,
  ranked_source_name = "Perturbation",
  geneset_source_name = "Aging"
)

# Visualize results
plot_fgsea_results(fgsea_res_perturb_centered, "Perturbation-Centered Analysis", output_dir, top_n_genesets_to_plot)

fgsea_res_perturb_centered_up_gs <- fgsea_res_perturb_centered %>% filter(geneset_direction == "UP")
plot_clustered_heatmap(fgsea_res_perturb_centered_up_gs, "Perturbation-Centered Analysis", "Aging", "UP", output_dir)

fgsea_res_perturb_centered_dn_gs <- fgsea_res_perturb_centered %>% filter(geneset_direction == "DN")
plot_clustered_heatmap(fgsea_res_perturb_centered_dn_gs, "Perturbation-Centered Analysis", "Aging", "DN", output_dir)


# Save the full results table
if (!is.null(fgsea_res_perturb_centered)) {
  write_tsv(fgsea_res_perturb_centered, file.path(output_dir, "fgsea_results_perturbation_centered.tsv"))
  message("Full Perturbation-Centered fgsea results saved to: ", file.path(output_dir, "fgsea_results_perturbation_centered.tsv"))
}

message("\n--- GSEA Combined Analysis Complete ---")