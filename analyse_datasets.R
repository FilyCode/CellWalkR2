library(dplyr)
library(tidyr)
library(ggplot2)
library(anndata)
library(stats)
library(energy)
library(Matrix)
library(irlba)


#' Analyze Perturbation Signatures from AnnData Objects
#'
#' This function processes OmicSignatureCollection objects in conjunction with AnnData
#' objects to compute and visualize perturbation signature metrics. It can compute
#' PCA from AnnData$X if `pca_key` is absent, extract differential expression (DE)
#' statistics, generate density plots of DE statistics, and compute Energy Distance
#' between perturbation groups and control groups using PCA components. It also
#' computes a simple TRADE-like approximate metric.
#'
#' @param omic_collection An object of class 'OmicSignatureCollection' containing
#'   multiple OmicSignature objects.
#' @param ann_data_object_paths A named list of file paths to AnnData (.h5ad) objects.
#'   Names of the list elements will be used as identifiers for the AnnData objects.
#' @param ann_data_obs_column_for_perturbation The name of the column in `adata$obs`
#'   that contains the perturbation labels (e.g., gene names, sgRNA IDs).
#' @param control_group_label The label within `ann_data_obs_column_for_perturbation`
#'   that designates the control group.
#' @param adj_p_threshold The adjusted p-value threshold to determine significant
#'   differential expression.
#' @param de_stat_for_plot The column name within the differential expression table
#'   to use for plotting and analysis (e.g., "logFC", "score").
#' @param pca_key If AnnData objects have pre-computed PCA stored in `adata$obsm`,
#'   specify the key (e.g., "X_pca"). If `NULL` and `e_distance_pc > 0`, PCA will
#'   be computed from `adata$X`.
#' @param e_distance_pc The number of principal components to use for Energy Distance
#'   calculation. If `0`, E-distance calculation is skipped.
#' @param min_cells_for_e_distance The minimum number of cells required in a group
#'   (perturbation or control) to compute Energy Distance.
#' @param compute_trade_approx Whether to compute the simple TRADE-like approximate
#'   metric (proportion of genes with large effect sizes).
#'
#' @return A list containing:
#'   \item{e_distances}{A data frame with computed energy distances for each
#'     perturbation group against the control group, or NULL if skipped.}
#'   \item{num_affected_genes}{A data frame summarizing the total number of genes,
#'     number of significantly affected genes, and TRADE approximation for each signature.}
#'   \item{de_stat_plots}{A list of ggplot objects: pooled DE statistic density,
#'     DE statistic density envelope, representative signature densities, and
#'     distribution of significant gene counts. NULL if no plots were generated.}
#'   \item{all_de_results_summary}{A data frame with comprehensive summary metrics
#'     for each processed signature.}
#'
#' @importFrom dplyr %>% arrange desc
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_density geom_line geom_ribbon geom_violin geom_boxplot geom_jitter labs theme_bw scale_color_viridis_d theme element_blank element_text
#' @importFrom anndata read_h5ad .adata.obs .adata.obsm .adata.X .adata.obsm_keys .adata.obs_keys
#' @importFrom stats density median sd approx quantile prcomp
#' @importFrom energy eudist energy.h.test
#' @importFrom Matrix colMeans
#' @importFrom irlba irlba
#'
#' @export
analyze_perturbation_signatures_from_anndata <- function(
    omic_collection,
    ann_data_object_paths = NULL,
    ann_data_obs_column_for_perturbation,
    control_group_label,
    adj_p_threshold = 0.05,
    de_stat_for_plot = "logFC",
    pca_key = "X_pca",
    e_distance_pc = 30,
    min_cells_for_e_distance = 10,
    compute_trade_approx = TRUE,
    max_pcs_compute = 50
) {
  # Basic checks on omic_collection
  if (is.null(omic_collection) || !("OmicSignatureCollection" %in% class(omic_collection)) && !("R6" %in% class(omic_collection))) {
    stop("`omic_collection` must be a non-empty OmicSignatureCollection object.")
  }
  if (length(omic_collection$OmicSigList) == 0) {
    stop("`omic_collection` contains no OmicSignatures.")
  }
  
  # Helper to extract DE table from an OmicSignature object robustly
  get_difexp_table <- function(sig_obj) {
    # Many OmicSignature objects store DE results in $difexp or slots; adapt as needed
    if (!is.null(sig_obj$difexp) && is.data.frame(sig_obj$difexp)) {
      return(sig_obj$difexp)
    }
    # try other common names
    if (!is.null(sig_obj$de) && is.data.frame(sig_obj$de)) return(sig_obj$de)
    if (!is.null(sig_obj$metadata) && !is.null(sig_obj$metadata$difexp) && is.data.frame(sig_obj$metadata$difexp)) return(sig_obj$metadata$difexp)
    return(NULL)
  }
  
  # Determine which DE stat column we'll use
  preferred_de_cols <- c(de_stat_for_plot, "logFC", "score", "t", "stat")
  
  # We'll collect outputs
  de_stats_for_plotting_data <- list()
  per_signature_summary <- list()
  e_dist_results <- list()
  
  # Load AnnData objects if provided
  loaded_ann_data <- list()
  if (!is.null(ann_data_object_paths)) {
    if (!is.list(ann_data_object_paths) || length(ann_data_object_paths) == 0) {
      stop("`ann_data_object_paths` must be a non-empty named list of file paths (or NULL).")
    }
    for (nm in names(ann_data_object_paths)) {
      path <- ann_data_object_paths[[nm]]
      if (!file.exists(path)) {
        warning(sprintf("AnnData file '%s' not found; skipping '%s'.", path, nm))
        next
      }
      ad <- tryCatch(anndata::read_h5ad(path), error = function(e) {
        warning(sprintf("Failed to read %s (%s): %s", nm, path, e$message)); NULL
      })
      if (!is.null(ad)) loaded_ann_data[[nm]] <- ad
    }
  }
  
  # --- Iterate OmicSignature objects: extract DE statistics and counts ---
  for (sig_name in names(omic_collection$OmicSigList)) {
    sig_obj <- omic_collection$OmicSigList[[sig_name]]
    difexp <- get_difexp_table(sig_obj)
    
    if (is.null(difexp)) {
      # collect blank summary but continue
      per_signature_summary[[sig_name]] <- data.frame(signature = sig_name,
                                                      n_genes = NA_integer_,
                                                      n_sig_genes = NA_integer_,
                                                      de_stat_col = NA_character_,
                                                      trade_approx = NA_real_,
                                                      stringsAsFactors = FALSE)
      next
    }
    
    # Find which column to use for plotting / stats
    de_col_found <- intersect(preferred_de_cols, colnames(difexp))
    if (length(de_col_found) == 0) {
      # attempt to find any numeric column that looks like a stat
      numeric_cols <- colnames(difexp)[sapply(difexp, is.numeric)]
      if (length(numeric_cols) == 0) {
        warning(sprintf("No numeric DE-stat columns found for signature '%s'; skipping plotting for this signature.", sig_name))
        per_signature_summary[[sig_name]] <- data.frame(signature = sig_name,
                                                        n_genes = nrow(difexp),
                                                        n_sig_genes = NA_integer_,
                                                        de_stat_col = NA_character_,
                                                        trade_approx = NA_real_,
                                                        stringsAsFactors = FALSE)
        next
      } else {
        chosen_col <- numeric_cols[1]
      }
    } else {
      chosen_col <- de_col_found[1]
    }
    
    # Ensure adjusted p-value column exists (common names)
    pval_col <- intersect(c("adj.P.Val", "padj", "adj_p", "p_adj", "p.adjust", "FDR"), colnames(difexp))
    pval_col <- if (length(pval_col) > 0) pval_col[1] else NULL
    
    # Collect DE stat vector and prepare for plotting
    de_vec <- difexp[[chosen_col]]
    # filter NAs
    de_vec <- de_vec[!is.na(de_vec)]
    
    # Store DE statistics for density plotting if sufficient data exists
    if (length(de_vec) > 2) {
      de_stats_for_plotting_data[[sig_name]] <- data.frame(signature = sig_name, de_statistic = de_vec, stringsAsFactors = FALSE)
    }
    
    # Count affected genes using adjusted p-value if available, else use effect-size threshold
    if (!is.null(pval_col)) {
      sig_count <- sum(!is.na(difexp[[pval_col]]) & difexp[[pval_col]] <= adj_p_threshold)
    } else {
      # fallback threshold: abs(effect) >= 1 (loosish)
      sig_count <- sum(!is.na(difexp[[chosen_col]]) & abs(difexp[[chosen_col]]) >= 1)
    }
    
    # TRADE-like approximate metric: proportion of genes with |effect| > median(|effect|) + 2*sd(|effect|)
    trade_approx_val <- NA_real_
    if (compute_trade_approx && length(de_vec) >= 10) {
      abs_effects <- abs(de_vec)
      thr <- stats::median(abs_effects, na.rm = TRUE) + 2 * stats::sd(abs_effects, na.rm = TRUE)
      trade_approx_val <- mean(abs_effects > thr, na.rm = TRUE) # proportion of genes with large effect
    }
    
    per_signature_summary[[sig_name]] <- data.frame(signature = sig_name,
                                                    n_genes = nrow(difexp),
                                                    n_sig_genes = sig_count,
                                                    de_stat_col = chosen_col,
                                                    trade_approx = trade_approx_val,
                                                    stringsAsFactors = FALSE)
  } # end signature loop
  
  # Combine signature summaries into a data.frame
  all_de_results_summary_df <- if (length(per_signature_summary) > 0) {
    dplyr::bind_rows(per_signature_summary)
  } else {
    data.frame()
  }
  
  # --- DE Statistic Density Plots Generation ---
  final_de_stat_plots <- NULL
  if (length(de_stats_for_plotting_data) > 0) {
    # Combine all per-signature vectors into one pooled data.frame
    plot_data_all <- dplyr::bind_rows(de_stats_for_plotting_data)
    
    # Basic pooled density plot
    plot_pooled <- ggplot2::ggplot(plot_data_all, ggplot2::aes(x = de_statistic)) +
      ggplot2::geom_density(fill = "grey70", alpha = 0.6, color = "black") +
      ggplot2::labs(title = paste0("Pooled distribution of ", de_stat_for_plot),
                    x = de_stat_for_plot, y = "Density") +
      ggplot2::theme_bw()
    
    # Compute density estimates for each signature on a common grid
    all_vals <- plot_data_all$de_statistic
    grid_x <- seq(stats::quantile(all_vals, 0.005, na.rm = TRUE),
                  stats::quantile(all_vals, 0.995, na.rm = TRUE),
                  length.out = 512)
    
    # Helper to get density on grid for a df
    get_density_on_grid <- function(df, xcol = "de_statistic", grid = grid_x) {
      v <- df[[xcol]]
      v <- v[!is.na(v)]
      if (length(v) < 5) return(rep(NA_real_, length(grid)))
      # Using bw = "nrd0" is a common default bandwidth selection
      d <- stats::density(v, from = min(grid), to = max(grid), n = length(grid), bw = "nrd0")
      # Interpolate to ensure values match the grid exactly
      density_vals <- approx(d$x, d$y, xout = grid, rule = 2)$y
      return(density_vals)
    }
    
    # Build matrix: rows = signatures, cols = grid points
    sig_names <- names(de_stats_for_plotting_data)
    dens_mat <- t(sapply(sig_names, function(s) get_density_on_grid(de_stats_for_plotting_data[[s]], "de_statistic", grid_x)))
    rownames(dens_mat) <- sig_names
    
    # Compute summary stats across signatures at each x: mean, 10th/90th percentiles
    mean_density <- apply(dens_mat, 2, mean, na.rm = TRUE)
    q10_density <- apply(dens_mat, 2, stats::quantile, probs = 0.10, na.rm = TRUE)
    q90_density <- apply(dens_mat, 2, stats::quantile, probs = 0.90, na.rm = TRUE)
    
    density_summary_df <- data.frame(x = grid_x,
                                     mean_density = mean_density,
                                     q10 = q10_density,
                                     q90 = q90_density)
    
    # Plot: mean density with 10-90% ribbon (envelope)
    plot_envelope <- ggplot2::ggplot(density_summary_df, ggplot2::aes(x = x)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q10, ymax = q90), fill = "grey80", alpha = 0.6) +
      ggplot2::geom_line(ggplot2::aes(y = mean_density), color = "firebrick", linewidth = 0.8) +
      ggplot2::labs(title = paste0("Signature density envelope (mean ± 10-90%): ", de_stat_for_plot),
                    x = de_stat_for_plot, y = "Density") +
      ggplot2::theme_bw()
    
    # Determine representative signatures: best, median, worst by chosen metric
    metric_df <- all_de_results_summary_df
    if (!("signature" %in% colnames(metric_df))) metric_df$signature <- rownames(metric_df)
    # ensure metric exists
    if (!"n_sig_genes" %in% colnames(metric_df)) metric_df$n_sig_genes <- NA_integer_
    if (!"trade_approx" %in% colnames(metric_df)) metric_df$trade_approx <- NA_real_
    
    # Compute fallback median_abs_effect if needed
    median_abs_effect <- sapply(names(de_stats_for_plotting_data), function(s) {
      vec <- de_stats_for_plotting_data[[s]]$de_statistic
      if (length(vec) == 0) return(NA_real_)
      stats::median(abs(vec), na.rm = TRUE)
    })
    med_df <- data.frame(signature = names(de_stats_for_plotting_data), median_abs_effect = median_abs_effect, stringsAsFactors = FALSE)
    # Use merge to ensure correct alignment and handling of signatures not in med_df (though they should match)
    metric_df <- merge(metric_df, med_df, by = "signature", all.x = TRUE)
    
    
    # Rank by n_sig_genes (descending), then trade_approx (descending), then median_abs_effect (descending)
    # Handle NAs gracefully for ranking
    rank_cols <- c("n_sig_genes", "trade_approx", "median_abs_effect")
    metric_df_ranked <- metric_df %>%
      dplyr::mutate(dplyr::across(
        dplyr::all_of(rank_cols),
        ~ ifelse(is.na(.), -Inf, .)
      )) %>%
      dplyr::arrange(dplyr::desc(n_sig_genes), dplyr::desc(trade_approx), dplyr::desc(median_abs_effect))
    
    # Select representative signatures
    n_sigs_available <- nrow(metric_df_ranked)
    example_sigs <- character(0)
    if (n_sigs_available > 0) {
      # Top signature
      example_sigs <- c(example_sigs, metric_df_ranked$signature[1])
      # Median signature (based on rank)
      median_idx <- max(1, ceiling(n_sigs_available / 2))
      example_sigs <- c(example_sigs, metric_df_ranked$signature[median_idx])
      # Bottom signature
      if (n_sigs_available > 1) {
        example_sigs <- c(example_sigs, metric_df_ranked$signature[n_sigs_available])
      }
    }
    example_sigs <- unique(na.omit(example_sigs)) # Remove duplicates and NAs
    
    # Plot their densities on the same axes for direct comparison
    example_plot_df <- do.call(rbind, lapply(example_sigs, function(s) {
      df <- de_stats_for_plotting_data[[s]]
      if (is.null(df)) return(NULL)
      # Use stats::density on the actual data points
      density_obj <- stats::density(df$de_statistic, na.rm = TRUE)
      data.frame(signature = s, x = density_obj$x, y = density_obj$y, stringsAsFactors = FALSE)
    }))
    
    plot_examples <- NULL
    if (!is.null(example_plot_df) && nrow(example_plot_df) > 0) {
      plot_examples <- ggplot2::ggplot(example_plot_df, ggplot2::aes(x = x, y = y, color = signature)) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::labs(title = "Representative signature density curves (best / median / worst)",
                      x = de_stat_for_plot, y = "Density") +
        ggplot2::theme_bw() +
        ggplot2::scale_color_viridis_d()
    }
    
    # Summary metrics plot: violin/boxplot of n_sig_genes across signatures
    plot_summary_metrics <- NULL
    if (nrow(metric_df) > 0) {
      # convert n_sig_genes to numeric and handle NAs for plotting
      metric_df$n_sig_genes_plot <- as.numeric(metric_df$n_sig_genes)
      plot_summary_metrics <- ggplot2::ggplot(metric_df, ggplot2::aes(x = 1, y = n_sig_genes_plot)) +
        ggplot2::geom_violin(fill = "lightblue", alpha = 0.6, trim = FALSE) + # trim=FALSE to show full distribution
        ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.8, fill = "white") +
        ggplot2::geom_jitter(width = 0.15, size = 0.6, alpha = 0.6) +
        ggplot2::labs(title = "Distribution of number of significant genes per signature",
                      x = "", y = "Number significant genes (adj p threshold)") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    }
    
    # Collect plots
    final_de_stat_plots <- list(
      pooled = plot_pooled,
      envelope = plot_envelope,
      examples = plot_examples,
      summary_metrics = plot_summary_metrics
    )
  } else {
    message("No DE-statistic vectors found to plot.")
  }
  
  # --- E-distance calculation ---
  # If e_distance_pc > 0 and AnnData objects provided, compute or extract PCA and compute energy-distance
  if (!is.null(loaded_ann_data) && length(loaded_ann_data) > 0 && e_distance_pc > 0) {
    for (adnm in names(loaded_ann_data)) {
      adata <- loaded_ann_data[[adnm]]
      
      # Ensure obs column exists
      if (!(ann_data_obs_column_for_perturbation %in% colnames(adata$obs))) {
        warning(sprintf("Obs column '%s' not in AnnData '%s'; skipping E-distance for this object.", ann_data_obs_column_for_perturbation, adnm))
        next
      }
      obs_col <- adata$obs[[ann_data_obs_column_for_perturbation]]
      
      # groups with sufficient cells
      grp_tab <- table(obs_col)
      valid_groups <- names(grp_tab[grp_tab >= min_cells_for_e_distance])
      if (!(control_group_label %in% names(grp_tab))) {
        warning(sprintf("Control label '%s' not found in AnnData '%s'. E-distance may be impossible.", control_group_label, adnm))
      }
      
      # Prepare PCA matrix:
      pca_mat <- NULL
      if (!is.null(pca_key) && pca_key %in% adata$obsm_keys()) {
        pca_mat <- adata$obsm[[pca_key]]
        message(sprintf("Using pre-computed PCA '%s' from AnnData '%s'.", pca_key, adnm))
        if (ncol(pca_mat) > e_distance_pc) {
          pca_mat <- pca_mat[, 1:e_distance_pc, drop = FALSE]
        }
      } else if (!is.null(adata$X)) {
        message(sprintf("PCA key '%s' not found or NULL in AnnData '%s'. Attempting to compute PCA from adata$X.", pca_key, adnm))
        
        target_pcs_to_compute <- min(e_distance_pc, max_pcs_compute)
        
        X <- adata$X
        pca_mat <- NULL
        
        tryCatch({
          if (inherits(X, "python.builtin.object")) {
            X_input <- as.matrix(X)
          } else if (is.matrix(X) || inherits(X, "sparseMatrix")) {
            X_input <- X
          } else {
            X_input <- as.matrix(X)
          }
          
          is_sparse_r <- inherits(X_input, "sparseMatrix")
          
          if (is_sparse_r) {
            message(sprintf("adata$X in '%s' is sparse. Using irlba for PCA.", adnm))
            if (nrow(X_input) < 2) {
              warning("Cannot compute PCA with fewer than 2 cells.")
              next
            }
            col_means <- Matrix::colMeans(X_input)
            X_centered <- sweep(X_input, 2, col_means, "-")
            
            actual_npc_to_compute <- min(nrow(X_input), ncol(X_input), target_pcs_to_compute)
            if (actual_npc_to_compute < 1) {
              warning(sprintf("Cannot compute PCA for AnnData '%s' (insufficient dimensions or requested PCs < 1); skipping E-distance for this object.", adnm))
              next
            }
            
            svd_res <- irlba::irlba(X_centered, nv = actual_npc_to_compute, maxiter = 1000, tol = 1e-6)
            pca_mat <- X_centered %*% svd_res$v[, 1:actual_npc_to_compute, drop = FALSE]
            message(sprintf("Computed PCA from sparse adata$X for '%s' (used %d PCs, requested %d).", adnm, actual_npc_to_compute, target_pcs_to_compute))
            
          } else {
            sparsity_threshold <- 0.95 
            num_zeros <- sum(X_input == 0)
            sparsity <- num_zeros / length(X_input)
            
            if (sparsity > sparsity_threshold) {
              message(sprintf("Dense adata$X in '%s' is highly sparse (%.2f%% zeros). Converting to sparse format and using irlba.", adnm, sparsity * 100))
              X_sparse <- Matrix::Matrix(X_input, sparse = TRUE)
              if (nrow(X_sparse) < 2) {
                warning("Cannot compute PCA with fewer than 2 cells.")
                next
              }
              col_means_sparse <- Matrix::colMeans(X_sparse)
              X_centered_sparse <- sweep(X_sparse, 2, col_means_sparse, "-")
              
              actual_npc_to_compute <- min(nrow(X_sparse), ncol(X_sparse), target_pcs_to_compute)
              if (actual_npc_to_compute < 1) {
                warning(sprintf("Cannot compute PCA for AnnData '%s' (insufficient dimensions or requested PCs < 1); skipping E-distance for this object.", adnm))
                next
              }
              
              svd_res <- irlba::irlba(X_centered_sparse, nv = actual_npc_to_compute, maxiter = 1000, tol = 1e-6)
              pca_mat <- X_centered_sparse %*% svd_res$v[, 1:actual_npc_to_compute, drop = FALSE]
              message(sprintf("Computed PCA from converted-sparse adata$X for '%s' (used %d PCs, requested %d).", adnm, actual_npc_to_compute, target_pcs_to_compute))
              
            } else {
              message(sprintf("Dense adata$X in '%s' (sparsity %.2f%%) is not highly sparse. Using irlba directly on dense matrix for PCA.", adnm, sparsity * 100))
              if (nrow(X_input) < 2) {
                warning("Cannot compute PCA with fewer than 2 cells.")
                next
              }
              col_means_dense <- colMeans(X_input)
              X_centered_dense <- X_input - matrix(col_means_dense, nrow = nrow(X_input), ncol = ncol(X_input), byrow = TRUE)
              
              actual_npc_to_compute <- min(ncol(X_input), nrow(X_input), target_pcs_to_compute)
              if (actual_npc_to_compute < 1) {
                warning(sprintf("Cannot compute PCA for AnnData '%s' (too few dimensions or requested PCs < 1); skipping E-distance for this object.", adnm))
                next
              }
              
              svd_res <- irlba::irlba(X_centered_dense, nv = actual_npc_to_compute, maxiter = 1000, tol = 1e-6)
              pca_mat <- X_centered_dense %*% svd_res$v[, 1:actual_npc_to_compute, drop = FALSE]
              message(sprintf("Computed PCA from dense adata$X for '%s' (used %d PCs, requested %d).", adnm, actual_npc_to_compute, target_pcs_to_compute))
            }
          }
        }, error = function(e) {
          warning(sprintf("Failed to compute PCA from AnnData$X for '%s': %s", adnm, e$message))
          pca_mat <- NULL
        })
      } else {
        warning(sprintf("No PCA key '%s' and no AnnData$X present for '%s' — cannot compute E-distance.", pca_key, adnm))
      }
      
      if (is.null(pca_mat)) {
        message(sprintf("Skipping E-distance for '%s' due to missing or failed PCA matrix computation.", adnm))
        next
      }
      
      pcs_to_use_for_edistance <- min(ncol(pca_mat), e_distance_pc) 
      if (pcs_to_use_for_edistance < 1) {
        warning(sprintf("Not enough PCA components (%d available) to satisfy e_distance_pc requirement (%d) for '%s'. Skipping E-distance.", ncol(pca_mat), e_distance_pc, adnm))
        next
      }
      pca_mat_subset <- pca_mat[, 1:pcs_to_use_for_edistance, drop = FALSE]
      
      if (nrow(pca_mat_subset) != nrow(adata$obs)) {
        warning(sprintf("PCA matrix rowcount (%d) doesn't match AnnData obs rows (%d) for '%s'. Skipping E-distance for this object.", nrow(pca_mat_subset), nrow(adata$obs), adnm))
        next
      }
      
      # Compute E-distance using energy::eudist + formula:
      # E(X, Y) = 2 * mean(||X - Y||) - mean(||X - X'||) - mean(||Y - Y'||)
      obs_col_chr <- as.character(obs_col)
      candidate_groups <- intersect(unique(obs_col_chr), valid_groups)
      
      for (grp in setdiff(candidate_groups, control_group_label)) {
        pert_idx <- which(obs_col_chr == grp)
        ctrl_idx <- which(obs_col_chr == control_group_label)
        
        if (length(pert_idx) < min_cells_for_e_distance || length(ctrl_idx) < min_cells_for_e_distance) {
          next
        }
        
        pert_pcs <- pca_mat_subset[pert_idx, , drop = FALSE]
        ctrl_pcs <- pca_mat_subset[ctrl_idx, , drop = FALSE]
        
        ed_val <- tryCatch({
          combined_pcs <- rbind(pert_pcs, ctrl_pcs)
          
          # Compute the full Euclidean distance matrix for the combined data
          dist_matrix_combined <- as.matrix(stats::dist(combined_pcs))
          
          # Extract the relevant sub-matrices:
          n_pert <- nrow(pert_pcs)
          n_ctrl <- nrow(ctrl_pcs)
          
          # d_xy: pairwise distances between perturbation group (first n_pert rows)
          #       and control group (next n_ctrl rows).
          d_xy <- dist_matrix_combined[1:n_pert, (n_pert + 1):(n_pert + n_ctrl), drop = FALSE]
          
          # d_xx: distances within the perturbation group.
          d_xx <- dist_matrix_combined[1:n_pert, 1:n_pert, drop = FALSE]
          
          # d_yy: distances within the control group.
          d_yy <- dist_matrix_combined[(n_pert + 1):(n_pert + n_ctrl), (n_pert + 1):(n_pert + n_ctrl), drop = FALSE]

          # off-diagonal means for within distances
          # mean() on a matrix calculates the mean of all its elements.
          # We need to exclude the diagonal elements for d_xx and d_yy before averaging.
          mean_xy <- mean(d_xy)
          mean_xx <- mean(d_xx[row(d_xx) != col(d_xx)]) # Excludes diagonal elements
          mean_yy <- mean(d_yy[row(d_yy) != col(d_yy)]) # Excludes diagonal elements
          
          # Calculate the energy distance component using the formula
          2 * mean_xy - mean_xx - mean_yy
          
        }, error = function(e) {
          warning(sprintf("Error computing energy distance for group '%s' in '%s': %s", grp, adnm, e$message))
          NA_real_ # Return NA if computation fails
        })
        
        e_dist_results[[paste0(grp, "_", adnm)]] <- data.frame(
          perturbation = grp,
          annobj = adnm,
          energy_distance = ed_val,
          stringsAsFactors = FALSE
        )
      }
    }
  } else {
    message("E-distance skipped: either no AnnData objects loaded or e_distance_pc <= 0.")
  }
  
  
  e_distances_df <- if (length(e_dist_results) > 0) {
    dplyr::bind_rows(e_dist_results)
  } else {
    NULL
  }
  
  # Number of affected genes summary is already in all_de_results_summary_df (per-signature)
  num_affected_genes_df <- all_de_results_summary_df[, c("signature", "n_genes", "n_sig_genes", "trade_approx")]
  
  return(list(
    e_distances = e_distances_df,
    num_affected_genes = num_affected_genes_df,
    de_stat_plots = final_de_stat_plots,
    all_de_results_summary = all_de_results_summary_df
  ))
}


e_dist_fig <- function(e_distances_df) {
  
  e_distance_plots <- NULL
  if (!is.null(e_distances_df) && nrow(e_distances_df) > 0) {
    # Order perturbations by energy distance within each AnnData object
    e_distances_df <- e_distances_df %>%
      dplyr::group_by(annobj) %>%
      dplyr::arrange(dplyr::desc(energy_distance), .by_group = TRUE) %>%
      dplyr::mutate(
        rank_within_ann = dplyr::row_number(),
        perturbation = factor(
          perturbation,
          levels = unique(perturbation[order(annobj, -energy_distance)])
        )
      ) %>%
      dplyr::ungroup()
    
    # 1) Dot plot: energy distance per perturbation, colored by annobj
    edist_dot <- ggplot2::ggplot(
      e_distances_df,
      ggplot2::aes(
        x = perturbation,
        y = energy_distance,
        color = annobj
      )
    ) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Energy distance to control per perturbation",
        x = "Perturbation",
        y = "Energy distance"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6)
      )
    
    # 2) Distribution plot: overall distribution per annobj
    edist_violin <- ggplot2::ggplot(
      e_distances_df,
      ggplot2::aes(x = annobj, y = energy_distance, fill = annobj)
    ) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.6) +
      ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white") +
      ggplot2::geom_jitter(width = 0.1, alpha = 0.5, size = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(
        title = "Distribution of energy distances to control",
        x = "AnnData object",
        y = "Energy distance"
      ) +
      ggplot2::theme_bw()
    
    # 3) If you want “top N” strongest perturbations per AnnData object
    top_n <- 30L
    e_distances_top <- e_distances_df %>%
      dplyr::group_by(annobj) %>%
      dplyr::slice_max(order_by = energy_distance, n = top_n, with_ties = FALSE) %>%
      dplyr::ungroup()
    
    edist_top_bar <- ggplot2::ggplot(
      e_distances_top,
      ggplot2::aes(
        x = reorder(perturbation, energy_distance),
        y = energy_distance,
        fill = annobj
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste0("Top ", top_n, " perturbations by energy distance to control"),
        x = "Perturbation",
        y = "Energy distance"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6)
      )
    
    e_distance_plots <- list(
      dot = edist_dot,
      violin = edist_violin,
      top_bar = edist_top_bar
    )
  }
  return(e_distance_plots)
}


# --- Main Script Execution ---

# Load the OmicSignatureCollection
# Ensure the path is correct and the file exists.
omic_collection <- readRDS("/restricted/projectnb/agedisease/projects/challenge2025/results/perturbational_omic_sigs/replogle_2022/Replogle_Perturb_Combined_OmicSignatureCollection.rds")

# Define paths to AnnData objects
k562_h5ad_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/replogle_2022/K562_essential_raw_singlecell_01.h5ad"
rpe1_h5ad_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/replogle_2022/rpe1_raw_singlecell_01.h5ad"
ann_data_input_paths <- list(K562 = k562_h5ad_path, RPE1 = rpe1_h5ad_path)

# --- Function Call ---
# If you want E-distance computed from AnnData$X (since PCA might be missing), ensure `pca_key = NULL` and `e_distance_pc > 0`.
# The function will attempt to compute PCA using irlba for sparse matrices or prcomp for dense.
print("Running analyze_perturbation_signatures_from_anndata...")
results <- analyze_perturbation_signatures_from_anndata(
  omic_collection = omic_collection,
  ann_data_object_paths = ann_data_input_paths,
  ann_data_obs_column_for_perturbation = "gene",
  control_group_label = "non-targeting",
  adj_p_threshold = 0.05,
  de_stat_for_plot = "logFC",
  pca_key = NULL,        # Set to NULL to allow function to compute PCA from adata$X if needed
  e_distance_pc = 30,    # Set to 0 to skip E-distance calculation entirely
  min_cells_for_e_distance = 10,
  compute_trade_approx = TRUE
)
print("Function execution complete.")

# --- Explore Results ---

# 1. Energy Distances
if (!is.null(results$e_distances) && nrow(results$e_distances) > 0) {
  print("--- Energy Distances (first few rows) ---")
  print(head(na.omit(results$e_distances)))
  
  e_dist_plots <- e_dist_fig(results$e_distances)
  if (!is.null(e_dist_plots)) {
    print("Energy distance dot plot:")
    print(e_dist_plots$dot)
    
    print("Energy distance distribution per AnnData object:")
    print(e_dist_plots$violin)
    
    print("Top perturbations by energy distance:")
    print(e_dist_plots$top_bar)
  }
  
} else {
  print("--- Energy Distance calculation was skipped or yielded no results. Check warnings above. ---")
}

# 2. Number of Affected Genes Summary
if (!is.null(results$num_affected_genes) && nrow(results$num_affected_genes) > 0) {
  print("--- Affected Genes Summary (first few rows) ---")
  print(head(results$num_affected_genes))
} else {
  print("--- No Affected Genes summary generated. ---")
}

# 3. DE Statistic Density Plots
if (!is.null(results$de_stat_plots) && length(results$de_stat_plots) > 0) {
  print("--- Displaying DE Statistic Density Plots ---")
  
  # Display pooled density plot
  if (!is.null(results$de_stat_plots$pooled)) {
    print("Pooled DE statistic density:")
    plot(results$de_stat_plots$pooled)
  }
  
  # Display density envelope plot
  if (!is.null(results$de_stat_plots$envelope)) {
    print("DE statistic density envelope:")
    plot(results$de_stat_plots$envelope)
  }
  
  # Display representative signature density plots
  if (!is.null(results$de_stat_plots$examples)) {
    print("Representative signature density curves:")
    plot(results$de_stat_plots$examples)
  }
  
  # Display summary metrics plot
  if (!is.null(results$de_stat_plots$summary_metrics)) {
    print("Distribution of number of significant genes per signature:")
    plot(results$de_stat_plots$summary_metrics)
  }
  
} else {
  print("--- No DE statistic density plots were generated. ---")
}

# 4. Comprehensive Summary of DE Results
if (!is.null(results$all_de_results_summary) && nrow(results$all_de_results_summary) > 0) {
  print("--- Comprehensive DE Results Summary (first few rows) ---")
  print(head(results$all_de_results_summary))
} else {
  print("--- No comprehensive DE results summary available. ---")
}

