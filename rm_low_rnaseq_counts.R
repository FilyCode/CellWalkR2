#' Remove those genes without at least 'min_thresh' (1) reads 
#' per 'reads_per' (million) in at least 'min_samples'
#' @author Stefano Monti (smonti at bu.edu)
#' @author Andrew Chen (andrewdr at bu.edu)
#' 
#' @export
rm_low_rnaseq_counts <- function(
    eset,
    min_samples = NULL,
    class_id = NULL,
    assay_id = NULL,
    reads_per = 1000000,
    min_thresh = 1) 
{
  ## BEGIN input checks
  stopifnot(methods::is(eset, "ExpressionSet") || 
              methods::is(eset, "SummarizedExperiment") || 
                methods::is(eset, "Seurat"))
  stopifnot(xor(is.null(min_samples), is.null(class_id)))
  ## END input checks

  if (methods::is(eset, "ExpressionSet")) 
  {
    if (!is.null(class_id)) {
      stopifnot(class_id %in% colnames(Biobase::pData(eset)))
      groups <- Biobase::pData(eset)[, class_id]
      min_samples <- max(min_thresh, table(groups))
    }
    rpm <- colSums(Biobase::exprs(eset)) / reads_per
    filter_ind <- t(apply(Biobase::exprs(eset), 1, function(x) {
      x > rpm
    }))
    filter_ind_rowsums <- apply(filter_ind, 1, sum)
    return(eset[filter_ind_rowsums >= min_samples, ])
  } 
  else if (methods::is(eset, "SummarizedExperiment")) 
  {
    if (is.null(assay_id)) {
      assay_id <- names(SummarizedExperiment::assays(eset))[1]
    }
    stopifnot(assay_id %in% names(SummarizedExperiment::assays(eset)))

    if (!is.null(class_id)) {
      stopifnot(class_id %in% colnames(SummarizedExperiment::colData(eset)))
      groups <- SummarizedExperiment::colData(eset)[, class_id]
      min_samples <- max(min_thresh, table(groups))
    }
    counts <- SummarizedExperiment::assays(eset)[[assay_id]]
    if (methods::is(counts, "DelayedMatrix")) {
      counts <- as.matrix(counts) # Convert to a dense matrix
    }
    # --- DEBUGGING START ---
    # Add a check here to ensure counts is numeric and doesn't have NAs/Infs
    if (!is.numeric(counts)) {
      stop("Counts matrix is not numeric after conversion. Current type: ", typeof(counts))
    }
    if (any(is.na(counts))) {
      message("Warning: Counts matrix contains NA values.")
    }
    if (any(!is.finite(counts))) {
      message("Warning: Counts matrix contains non-finite values (Inf, -Inf, NaN).")
    }
    
    
    rpm <- colSums(counts) / reads_per
    filter_ind <- t(apply(counts, 1, function(x) {
      x > rpm
    }))
    filter_ind_rowsums <- apply(filter_ind, 1, sum)
    return(eset[filter_ind_rowsums >= min_samples, ])
  } 
  else if (methods::is(eset, "Seurat")) 
  {
    # Default to RNA assay if not specified
    if (is.null(assay_id)) assay_id <- "RNA"
    stopifnot(assay_id %in% names(eset@assays))
    counts <- Seurat::GetAssayData(eset, assay = assay_id, slot = "counts")
    meta <- eset@meta.data

    # If class_id is specified, use it to determine min_samples
    if (!is.null(class_id)) {
      stopifnot(class_id %in% colnames(meta))
      groups <- meta[, class_id]
      min_samples <- max(min_thresh, table(groups))
    }
    rpm <- colSums(counts) / reads_per
    filter_ind <- t(apply(counts, 1, function(x) {
      x > rpm
    }))
    filter_ind_rowsums <- rowSums(filter_ind)
    keep_genes <- which(filter_ind_rowsums >= min_samples)
    # Subset Seurat object to keep only the filtered genes
    return(subset(eset, features = rownames(counts)[keep_genes]))
  } 
  else {
    stop("unrecognized oject type: ", class(eset))
  }
}
