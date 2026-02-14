# Updated analyze_perturbation_signatures_from_anndata
# - Computes PCA from AnnData$X if `pca_key` absent and e_distance_pc > 0
# - Always extracts DE statistics & counts affected genes from OmicSignature objects
# - Produces DE-statistic density plots (combined + per-signature) when possible
# - Computes Energy Distance when possible (requires cell counts >= min_cells_for_e_distance)
# - Adds a simple TRADE-like approximate metric (optional, conservative)
#
# Dependencies: dplyr, tidyr, ggplot2, anndata, energy, Matrix (if sparse matrices), stats
#
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
  
  # --- checks / required packages ---
  required_packages <- c("dplyr","tidyr","ggplot2","anndata","energy","stats")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' required but not installed.", pkg))
    }
  }
  
  # Basic checks on omic_collection
  if (is.null(omic_collection) || !"OmicSignatureCollection" %in% class(omic_collection) && !("R6" %in% class(omic_collection))) {
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
  
  # Iterate OmicSignature objects: extract de-statistics and counts
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
    
    # Collect de stat vector and prepare for plotting
    de_vec <- difexp[[chosen_col]]
    # filter NAs
    de_vec <- de_vec[!is.na(de_vec)]
    if (length(de_vec) > 2) {
      de_stats_for_plotting_data[[sig_name]] <- data.frame(signature = sig_name, de_statistic = de_vec, stringsAsFactors = FALSE)
    }
    
    # Count affected genes using adjusted p-value if available, else use effect-size threshold
    if (!is.null(pval_col)) {
      sig_count <- sum(!is.na(difexp[[pval_col]]) & difexp[[pval_col]] <= adj_p_threshold)
    } else {
      # fallback threshold: abs(effect) > 1 (loosish)
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
  all_de_results_summary_df <- do.call(rbind, per_signature_summary)
  
  
  final_de_stat_plots <- NULL
  if (length(de_stats_for_plotting_data) > 0) {
    # Combine all per-signature vectors into one pooled data.frame
    plot_data_all <- do.call(rbind, de_stats_for_plotting_data)
    # basic pooled density plot
    plot_pooled <- ggplot2::ggplot(plot_data_all, ggplot2::aes(x = de_statistic)) +
      ggplot2::geom_density(fill = "grey70", alpha = 0.6, color = "black") +
      ggplot2::labs(title = paste0("Pooled distribution of ", de_stat_for_plot),
                    x = de_stat_for_plot, y = "Density") +
      ggplot2::theme_bw()
    
    # Compute density estimates for each signature on a common grid so we can compute mean/quantiles
    # Define grid
    all_vals <- plot_data_all$de_statistic
    grid_x <- seq(quantile(all_vals, 0.005, na.rm = TRUE),
                  quantile(all_vals, 0.995, na.rm = TRUE),
                  length.out = 512)
    
    # Helper to get density on grid for a df
    get_density_on_grid <- function(df, xcol = "de_statistic", grid = grid_x) {
      v <- df[[xcol]]
      v <- v[!is.na(v)]
      if (length(v) < 5) return(rep(NA_real_, length(grid)))
      d <- stats::density(v, from = min(grid), to = max(grid), n = length(grid), bw = "nrd0")
      # d$x should match grid if from/to/n are set, but to be robust we'll interpolate
      density_vals <- approx(d$x, d$y, xout = grid, rule = 2)$y
      return(density_vals)
    }
    
    # Build matrix: rows = signatures, cols = grid points
    sig_names <- names(de_stats_for_plotting_data)
    dens_mat <- t(sapply(sig_names, function(s) get_density_on_grid(de_stats_for_plotting_data[[s]], "de_statistic", grid_x)))
    rownames(dens_mat) <- sig_names
    
    # compute summary stats across signatures at each x: mean, 10th/90th percentiles
    mean_density <- apply(dens_mat, 2, function(x) mean(x, na.rm = TRUE))
    q10_density <- apply(dens_mat, 2, function(x) quantile(x, probs = 0.10, na.rm = TRUE))
    q90_density <- apply(dens_mat, 2, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
    
    density_summary_df <- data.frame(x = grid_x,
                                     mean_density = mean_density,
                                     q10 = q10_density,
                                     q90 = q90_density)
    
    # Plot: mean density with 10-90% ribbon (envelope)
    plot_density_envelope <- ggplot2::ggplot(density_summary_df, ggplot2::aes(x = x)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q10, ymax = q90), fill = "grey80", alpha = 0.6) +
      ggplot2::geom_line(ggplot2::aes(y = mean_density), color = "firebrick", linewidth = 0.8) +
      ggplot2::labs(title = paste0("Signature density envelope (mean ± 10-90%): ", de_stat_for_plot),
                    x = de_stat_for_plot, y = "Density") +
      ggplot2::theme_bw()
    
    # Choose representative signatures: best, median, worst by chosen metric
    # Metric: prefer n_sig_genes; fallback to trade_approx; fallback to median absolute effect
    metric_df <- all_de_results_summary_df
    if (!("signature" %in% colnames(metric_df))) metric_df$signature <- rownames(metric_df)
    # ensure metric exists
    if (!"n_sig_genes" %in% colnames(metric_df)) metric_df$n_sig_genes <- NA_integer_
    if (!"trade_approx" %in% colnames(metric_df)) metric_df$trade_approx <- NA_real_
    # compute fallback median_abs_effect if needed
    median_abs_effect <- sapply(names(de_stats_for_plotting_data), function(s) median(abs(de_stats_for_plotting_data[[s]]$de_statistic), na.rm = TRUE))
    med_df <- data.frame(signature = names(de_stats_for_plotting_data), median_abs_effect = median_abs_effect, stringsAsFactors = FALSE)
    metric_df <- merge(metric_df, med_df, by = "signature", all.y = TRUE)
    
    # rank by n_sig_genes (descending), pick top (best) and bottom (worst) and median
    metric_df <- metric_df[order(-as.numeric(metric_df$n_sig_genes), decreasing = FALSE), ] # keep stable ordering
    # remove NA signatures if needed
    metric_df$rank_metric <- rank(-as.numeric(ifelse(is.na(metric_df$n_sig_genes), -Inf, as.numeric(metric_df$n_sig_genes))), ties.method = "first")
    # sort by that rank
    metric_df <- metric_df[order(metric_df$rank_metric), ]
    
    # define examples:
    # best = highest n_sig_genes (if present), worst = lowest, median = median rank
    valid_metric_idx <- which(!is.na(metric_df$n_sig_genes) & is.finite(metric_df$n_sig_genes))
    if (length(valid_metric_idx) >= 3) {
      best_sig <- metric_df$signature[which.max(metric_df$n_sig_genes)]
      worst_sig <- metric_df$signature[which.min(metric_df$n_sig_genes)]
      median_sig <- metric_df$signature[order(metric_df$n_sig_genes)][ceiling(length(valid_metric_idx)/2)]
    } else {
      # fallback to median_abs_effect
      metric_df2 <- metric_df[order(metric_df$median_abs_effect, decreasing = TRUE), ]
      best_sig <- metric_df2$signature[1]
      worst_sig <- metric_df2$signature[nrow(metric_df2)]
      median_sig <- metric_df2$signature[max(1, ceiling(nrow(metric_df2)/2))]
    }
    
    example_sigs <- unique(na.omit(c(best_sig, median_sig, worst_sig)))
    # Plot their densities on the same axes for direct comparison
    example_plot_df <- do.call(rbind, lapply(example_sigs, function(s) {
      df <- de_stats_for_plotting_data[[s]]
      if (is.null(df)) return(NULL)
      d <- stats::density(df$de_statistic, na.rm = TRUE)
      data.frame(signature = s, x = d$x, y = d$y, stringsAsFactors = FALSE)
    }))
    plot_examples <- NULL
    if (!is.null(example_plot_df) && nrow(example_plot_df) > 0) {
      plot_examples <- ggplot2::ggplot(example_plot_df, ggplot2::aes(x = x, y = y, color = signature)) +
        ggplot2::geom_line(size = 0.9) +
        ggplot2::labs(title = "Representative signature density curves (best / median / worst)",
                      x = de_stat_for_plot, y = "Density") +
        ggplot2::theme_bw() +
        ggplot2::scale_color_viridis_d()
    }
    
    # Additionally, a compact summary scatter/violin of per-signature metrics:
    # We'll use n_sig_genes and trade_approx (if available)
    summary_metrics_df <- metric_df
    # plot: violin/boxplot of n_sig_genes across signatures
    plot_summary_metrics <- NULL
    if (nrow(summary_metrics_df) >= 5) {
      # convert n_sig_genes to numeric and handle NAs
      summary_metrics_df$n_sig_genes_num <- as.numeric(summary_metrics_df$n_sig_genes)
      plot_summary_metrics <- ggplot2::ggplot(summary_metrics_df, ggplot2::aes(x = 1, y = n_sig_genes_num)) +
        ggplot2::geom_violin(fill = "lightblue", alpha = 0.6) +
        ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.8) +
        ggplot2::geom_jitter(width = 0.15, size = 0.6, alpha = 0.6) +
        ggplot2::labs(title = "Distribution of number of significant genes per signature",
                      x = "", y = "Number significant genes (adj p threshold)") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    }
    
    # Collect plots
    final_de_stat_plots <- list(
      pooled = plot_pooled,
      envelope = plot_density_envelope,
      examples = plot_examples,
      summary_metrics = plot_summary_metrics
    )
  } else {
    message("No DE-statistic vectors found to plot.")
  }
  }
  
  # E-distance calculation:
  # If e_distance_pc > 0 and AnnData objects provided, compute or extract PCA and compute energy_distance
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
      
      # prepare PCA matrix:
      pca_mat <- NULL
      if (!is.null(pca_key) && pca_key %in% adata$obsm_keys()) {
        pca_mat <- adata$obsm[[pca_key]]
      } else {
        # attempt to compute PCA from adata$X (if present)
        if (!is.null(adata$X)) {
          X <- adata$X
          # if sparse from reticulate, convert to R Matrix or numeric matrix carefully
          # try to coerce to dense if small; otherwise use prcomp on sparse (via svds if matrix package available)
          # For simplicity, try prcomp on dense (may be memory heavy)
          tryCatch({
            if ("dgCMatrix" %in% class(X) || inherits(X, "matrix")) {
              Xdense <- as.matrix(X)
            } else {
              # if python sparse matrix object, try converting via reticulate interface
              Xdense <- tryCatch({
                as.matrix(X)
              }, error = function(e) {
                stop("Cannot coerce AnnData$X to R matrix for PCA. Provide PCA in .obsm or set e_distance_pc = 0.")
              })
            }
            # center & scale (prcomp does centering)
            npc <- min(ncol(Xdense), nrow(Xdense), max_pcs_compute, e_distance_pc)
            if (npc < 2) {
              warning(sprintf("Cannot compute PCA for AnnData '%s' (too few dimensions); skipping E-distance for this object.", adnm))
              next
            }
            pr <- stats::prcomp(Xdense, center = TRUE, scale. = FALSE, retx = TRUE)
            pca_mat <- pr$x[, 1:min(ncol(pr$x), e_distance_pc), drop = FALSE]
            message(sprintf("Computed PCA from AnnData$X for '%s' (using %d PCs).", adnm, ncol(pca_mat)))
          }, error = function(e) {
            warning(sprintf("Failed to compute PCA from AnnData$X for '%s': %s", adnm, e$message))
            pca_mat <- NULL
          })
        } else {
          warning(sprintf("No PCA key and no AnnData$X present for '%s' — cannot compute E-distance.", adnm))
        }
      }
      
      if (is.null(pca_mat)) next
      
      # iterate perturbation groups that correspond to signatures
      # We'll attempt E-distance for any group name that appears in omic signatures parsed names
      # Create vector of candidate group names (intersect with valid_groups)
      candidate_groups <- intersect(unique(as.character(obs_col)), valid_groups)
      # For each candidate (except control) compute energy_distance between its cells and control cells
      for (grp in setdiff(candidate_groups, control_group_label)) {
        pert_idx <- which(obs_col == grp)
        ctrl_idx <- which(obs_col == control_group_label)
        if (length(pert_idx) < min_cells_for_e_distance || length(ctrl_idx) < min_cells_for_e_distance) {
          next
        }
        # ensure pca_mat rows correspond to adata cells (AnnData read_h5ad should preserve order)
        # Subset rows
        # safety: ensure pca_mat has same number of rows as nrow(adat$obs)
        if (nrow(pca_mat) != nrow(adata$obs)) {
          warning(sprintf("PCA matrix rowcount (%d) doesn't match AnnData obs rows (%d) for '%s'. Skipping E-distance for this object.", nrow(pca_mat), nrow(adata$obs), adnm))
          next
        }
        pcs_to_use <- min(ncol(pca_mat), e_distance_pc)
        pert_pcs <- pca_mat[pert_idx, 1:pcs_to_use, drop = FALSE]
        ctrl_pcs <- pca_mat[ctrl_idx, 1:pcs_to_use, drop = FALSE]
        
        # compute energy distance
        ed_val <- tryCatch({
          energy::energy_distance(pert_pcs, ctrl_pcs)
        }, error = function(e) {
          warning(sprintf("Error computing energy_distance for group '%s' in '%s': %s", grp, adnm, e$message))
          NA_real_
        })
        e_dist_results[[paste0(grp, "_", adnm)]] <- data.frame(perturbation = grp, annobj = adnm, energy_distance = ed_val, stringsAsFactors = FALSE)
      } # end group loop
    } # end annData loop
  } else {
    message("E-distance skipped: either no AnnData objects loaded or e_distance_pc <= 0.")
  }
  
  e_distances_df <- if (length(e_dist_results) > 0) do.call(rbind, e_dist_results) else NULL
  
  # Count affected genes summary is already in all_de_results_summary_df (per-signature)
  num_affected_genes_df <- all_de_results_summary_df[, c("signature", "n_genes", "n_sig_genes", "trade_approx")]
  
  return(list(
    e_distances = e_distances_df,
    num_affected_genes = num_affected_genes_df,
    de_stat_plots = final_de_stat_plots,
    pca_plots = NULL, # we didn't create visual PCA objects here; the function can be extended to create them if desired
    all_de_results_summary = all_de_results_summary_df
  ))
}




omic_collection <- readRDS("/restricted/projectnb/agedisease/projects/challenge2025/results/perturbational_omic_sigs/replogle_2022/Replogle_Perturb_Combined_OmicSignatureCollection.rds")

k562_h5ad_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/replogle_2022/K562_essential_raw_singlecell_01.h5ad"
rpe1_h5ad_path <- "/restricted/projectnb/agedisease/CBMrepositoryData/replogle_2022/rpe1_raw_singlecell_01.h5ad"
ann_data_input_paths <- list(K562 = k562_h5ad_path, RPE1 = rpe1_h5ad_path)

ann_data_obs_perturb_col <- 
seurat_metadata_control_label <- 

# If you want E-distance computed from AnnData$X (since PCA is missing), set pca_key = NULL and e_distance_pc > 0
results <- analyze_perturbation_signatures_from_anndata(
  omic_collection = omic_collection,
  ann_data_object_paths = ann_data_input_paths,
  ann_data_obs_column_for_perturbation = "gene",
  control_group_label = "non-targeting",
  adj_p_threshold = 0.05,
  de_stat_for_plot = "logFC",
  pca_key = NULL,        # allow function to compute PCA from adata$X if needed
  e_distance_pc = 30,    # set to 0 to skip E-distance entirely
  min_cells_for_e_distance = 10,
  compute_trade_approx = TRUE
)

# Explore results
if (!is.null(results$e_distances)) {
  print("Energy distances (subset):")
  print(head(results$e_distances))
} else {
  print("No E-distance results.")
}
print("Affected genes summary:")
print(head(results$num_affected_genes))
if (!is.null(results$de_stat_plots)) {
  # display combined density
  print(results$de_stat_plots$all_combined)
}
print("DE results summary (first rows):")
print(head(results$all_de_results_summary))


# --- Explore the results ---
# 1. Energy Distances (will be NULL as E-distance is skipped)
if (!is.null(results$e_distances)) {
  print("Energy Distances:")
  print(results$e_distances)
} else {
  print("Energy Distance calculation was skipped as PCA was not available or e_distance_pc was 0. Check warnings above.")
}

# 2. Number of Affected Genes
print("Number of Affected Genes:")
print(results$num_affected_genes)

# 3. DE Statistic Density Plots
if (length(results$de_stat_plots) > 0) {
  print("Differential Expression Statistic Density Plots:")
  # You can view these plots:
  # print(results$de_stat_plots[["all_combined"]]) # Combined plot for all perturbations
  # plot(results$de_stat_plots[["NAF1"]])         # Individual plot for NAF1 perturbation (if found)
} else {
  print("No DE statistic density plots generated.")
}

# 4. PCA Plots (will be NULL as PCA is absent)
if (!is.null(results$pca_plots) && length(results$pca_plots) > 0) {
  print("PCA Plots:")
  # Example: print(results$pca_plots[["K562"]])
} else {
  print("PCA plots were not generated as PCA data was absent.")
}

# 5. Summary of DE results for all processed signatures
print("Summary of Differential Expression Results:")
print(results$all_de_results_summary)

