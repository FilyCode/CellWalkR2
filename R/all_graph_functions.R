library(CellWalkR)
library(ArchR) # Only for ArchR-specific functions if used; data extraction will use zellkonverter/Signac
library(Seurat)
library(Signac)
library(zellkonverter) # For reading H5AD files
library(Matrix)
library(matrixStats)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(FNN)
library(irlba)
library(pheatmap)
library(igraph)
library(cowplot) # For plot_grid
library(dplyr)   # For %>% and filter
library(future.apply)
set.seed(123) # For reproducibility



#' Create Spatial KNN Graph
#'
#' \code{createSpatialKNNGraph()} computes a cell-to-cell similarity matrix based on spatial proximity
#' using a K-nearest neighbors approach. This function is now intended for *validation purposes*
#' or for scenarios where spatial proximity is explicitly used as a modality.
#'
#' @param spatialCoords matrix or data.frame with columns for spatial x and y coordinates,
#'                      rownames should be cell barcodes.
#' @param k_neighbors integer, the number of nearest neighbors to consider for each cell.
#' @return a sparseMatrix of cell-to-cell similarity, with values normalized between 0 and 1.
#' @export
createSpatialKNNGraph <- function(spatialCoords, k_neighbors = 10) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required for spatial KNN graph. Please install it with install.packages('FNN')")
  }
  if (missing(spatialCoords) || !(is.matrix(spatialCoords) || is.data.frame(spatialCoords))) {
    stop("Must provide 'spatialCoords' as a matrix or data.frame with cell barcodes as rownames.")
  }
  if (ncol(spatialCoords) < 2) {
    stop("'spatialCoords' must have at least two columns (e.g., x and y coordinates).")
  }
  if (is.null(rownames(spatialCoords))) {
    stop("'spatialCoords' must have rownames (cell barcodes).")
  }
  if (nrow(spatialCoords) <= 1) {
    message("Warning: Too few cells (<= 1) for spatial KNN graph computation. Returning empty or identity matrix.")
    if (nrow(spatialCoords) == 0) return(Matrix::Matrix(0, 0, 0, sparse = TRUE, dimnames = list(NULL, NULL)))
    return(Matrix::Matrix(1, 1, 1, sparse = TRUE, dimnames = list(rownames(spatialCoords), rownames(spatialCoords))))
  }

  if (k_neighbors >= nrow(spatialCoords)) {
    warning("'k_neighbors' is greater than or equal to the number of cells. Setting k_neighbors to num_cells - 1.")
    k_neighbors <- nrow(spatialCoords) - 1
  }
  if (k_neighbors < 1) {
    stop("'k_neighbors' must be at least 1 after adjustments.")
  }

  message("Computing spatial distances and KNN graph...")

  knn_result <- FNN::get.knn(spatialCoords, k = k_neighbors)

  adj_list_i <- rep(1:nrow(spatialCoords), each = k_neighbors)
  adj_list_j <- as.vector(knn_result$index)
  distances <- as.vector(knn_result$distance)

  max_dist_in_knn <- max(distances)
  if (max_dist_in_knn == 0) {
    similarity <- rep(1, length(distances))
  } else {
    similarity <- 1 - (distances / max_dist_in_knn)
  }

  valid_indices <- adj_list_j <= nrow(spatialCoords) & adj_list_j >= 1
  adj_list_i <- adj_list_i[valid_indices]
  adj_list_j <- adj_list_j[valid_indices]
  similarity <- similarity[valid_indices]

  spatial_graph <- Matrix::sparseMatrix(
    i = adj_list_i,
    j = adj_list_j,
    x = similarity,
    dims = c(nrow(spatialCoords), nrow(spatialCoords)),
    dimnames = list(rownames(spatialCoords), rownames(spatialCoords))
  )

  spatial_graph <- pmax(spatial_graph, t(spatial_graph))
  diag(spatial_graph) <- 1

  message("Spatial KNN graph created.")
  return(spatial_graph)
}

#' Combine Multiple Cell-to-Cell Graphs
#'
#' \code{combineMultiModalCellGraphs()} takes a list of cell-to-cell similarity matrices
#' and combines them into a single matrix using a weighted sum. This allows integrating
#' similarity information from different modalities (e.g., RNA, ATAC).
#'
#' @param listOfCellGraphs list of sparseMatrix or matrix objects, each representing
#'                         cell-to-cell similarity. Rownames and colnames must be consistent
#'                         across all matrices and represent cell barcodes.
#' @param weights numeric vector of weights, corresponding to the order of matrices
#'                in `listOfCellGraphs`. Weights will be normalized to sum to 1.
#' @return a sparseMatrix representing the combined cell-to-cell similarity graph.
#' @export
combineMultiModalCellGraphs <- function(listOfCellGraphs, weights) {
  if (missing(listOfCellGraphs) || !is.list(listOfCellGraphs) || length(listOfCellGraphs) == 0) {
    stop("Must provide 'listOfCellGraphs' as a non-empty list of matrices.")
  }
  if (missing(weights) || !is.numeric(weights) || length(weights) != length(listOfCellGraphs)) {
    stop("Must provide 'weights' as a numeric vector of the same length as 'listOfCellGraphs'.")
  }
  if (sum(weights) != 1) {
    message("Weights do not sum to 1. Normalizing weights to sum to 1.")
    weights <- weights / sum(weights)
  }

  first_graph <- listOfCellGraphs[[1]]
  if (is.null(first_graph)) {
    stop("The first graph in 'listOfCellGraphs' is NULL or empty. Ensure all graphs are valid matrices and non-empty.")
  }
  num_cells <- nrow(first_graph)
  cell_names <- rownames(first_graph)

  if (is.null(cell_names) || length(cell_names) == 0) {
    stop("The first graph in 'listOfCellGraphs' has no rownames (cell barcodes).")
  }

  for (i in 2:length(listOfCellGraphs)) {
    current_graph <- listOfCellGraphs[[i]]
    if (is.null(current_graph) || nrow(current_graph) == 0) {
      stop(paste0("Graph ", i, " in 'listOfCellGraphs' is NULL or empty. Ensure all graphs are valid matrices and non-empty."))
    }

    # Use isTRUE(all.equal(...)) for robust comparison of dimensions
    if (!isTRUE(all.equal(dim(current_graph), c(num_cells, num_cells)))) {
      stop(paste0("Graph ", i, " has inconsistent dimensions (", paste(dim(current_graph), collapse="x"),
                  ") compared to the first graph (", num_cells, "x", num_cells, "). All graphs must have the same number of cells. ",
                  "Current logic expects all input graphs to share the exact same cell identities and order. ",
                  "If you have non-overlapping cell sets, consider a spatial integration/imputation strategy first to unify the cell universe."))
    }
    # Check if cell names are exactly the same and in the same order
    if (!identical(rownames(current_graph), cell_names) || !identical(colnames(current_graph), cell_names)) {
      stop(paste0("Graph ", i, " has inconsistent cell names or order. All graphs must have the same rownames/colnames in the same order. Please reorder matrices if necessary."))
    }
  }

  message("Combining cell graphs using weighted sum...")
  combined_graph <- Matrix::Matrix(0, nrow = num_cells, ncol = num_cells, sparse = TRUE, dimnames = list(cell_names, cell_names))

  for (i in 1:length(listOfCellGraphs)) {
    combined_graph <- combined_graph + (listOfCellGraphs[[i]] * weights[i])
  }

  message("Multi-modal cell graph created.")
  return(combined_graph)
}


#' Read 10x Xenium Spatial RNA-seq Data
#'
#' \code{readXeniumData()} reads 10x Xenium H5 file for gene counts and
#' `cells.csv` for spatial metadata into a Seurat object.
#'
#' @param rna_h5_path Path to the Xenium cell_feature_matrix.h5 file.
#' @param cells_csv_path Path to the Xenium cells.csv.gz file.
#' @param project_name Project name for the Seurat object.
#' @param min_cells Include features (genes) detected in at least this many cells.
#' @param min_features Include cells with at least this many features (genes).
#' @return A Seurat object containing RNA counts and spatial metadata.
#' @export
readXeniumData <- function(rna_h5_path, cells_csv_path, project_name = "Xenium_RNA", min_cells = 3, min_features = 100) {
  message("Reading Xenium spatial RNA-seq data from: ", rna_h5_path)

  raw_mat_list <- Read10X_h5(rna_h5_path)

  # Explicitly select the gene expression matrix, as 10x Xenium H5 contains multiple types
  if (!("Gene Expression" %in% names(raw_mat_list))) {
    stop("Expected 'Gene Expression' matrix not found in the H5 file. Found: ", paste(names(raw_mat_list), collapse = ", "))
  }
  mat <- raw_mat_list[["Gene Expression"]]

  # Check for empty matrix or zero dimensions/entries
  if (prod(dim(mat)) == 0 || nrow(mat) == 0 || ncol(mat) == 0) {
    stop("`mat` is an empty or zero-dimension matrix after extracting 'Gene Expression'. Cannot create Seurat object.")
  }
  if (sum(mat > 0) == 0) {
    warning("`mat` contains only zero entries after extraction. This might cause issues downstream.")
  }
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("`mat` is missing rownames (features/genes) or colnames (cells) after extraction. Cannot create Seurat object.")
  }

  # Create Seurat object with initial filtering
  seurat_rna <- CreateSeuratObject(
    counts = mat,
    project = project_name,
    min.cells = min_cells,    # Filter genes (features)
    min.features = min_features # Filter cells
  )

  message("Seurat object created with initial filtering. Original cells: ", ncol(mat), ", original features: ", nrow(mat),
          ". After filtering (min.cells=", min_cells, ", min.features=", min_features, "): ", ncol(seurat_rna), " cells, ", nrow(seurat_rna), " features.")

  # Add cell metadata (centroids etc.)
  meta <- read.csv(cells_csv_path)
  rownames(meta) <- meta$cell_id

  common_cells_meta_rna <- intersect(colnames(seurat_rna), rownames(meta))
  if (length(common_cells_meta_rna) == 0) {
    stop("No common cells found between filtered Seurat object and cells.csv. Check data integrity.")
  }
  seurat_rna <- subset(seurat_rna, cells = common_cells_meta_rna)
  seurat_rna <- AddMetaData(seurat_rna, metadata = meta[common_cells_meta_rna, , drop = FALSE])

  message("Xenium data loaded into Seurat object. Final object contains ", ncol(seurat_rna), " cells and ", nrow(seurat_rna), " features.")
  return(seurat_rna)
}


#' Read Spatial ATAC-seq Data from AnnData (.h5ad)
#'
#' \code{readAnnDataATAC()} reads spatial scATAC-seq data from an AnnData (.h5ad) file
#' using `zellkonverter`, extracting the peak matrix, peak GRanges, and spatial coordinates.
#'
#' @param atac_h5ad_path Path to the AnnData .h5ad file containing ATAC-seq data.
#' @param x_col Name of the column in `colData` for X-coordinates (e.g., "xcor").
#' @param y_col Name of the column in `colData` for Y-coordinates (e.g., "ycor").
#' @param genome Genome build (e.g., "hg38"). Required for `GRanges`.
#' @return A list containing:
#'   \itemize{
#'     \item `counts_matrix`: A sparse matrix of peak counts (cells x peaks).
#'     \item `peaks_gr`: A `GRanges` object of peak genomic regions.
#'     \item `spatial_coords`: A data.frame of spatial x, y coordinates.
#'     \item `cell_names`: A character vector of cell barcodes.
#'   }
#' @export
readAnnDataATAC <- function(atac_h5ad_path, x_col = "xcor", y_col = "ycor", genome = "hg38") {
  message("Reading spatial ATAC-seq data from AnnData (.h5ad)...")
  if (!file.exists(atac_h5ad_path)) {
    stop(paste0("AnnData file not found at: ", atac_h5ad_path))
  }

  sce <- readH5AD(atac_h5ad_path)

  # SingleCellExperiment assays are always FEATURES x CELLS.
  # CellWalker2 expects count matrices to be CELLS x FEATURES.
  # Therefore, we MUST transpose the assay matrix.
  counts_atac_feature_by_cell <- as(assay(sce, "X"), "dgCMatrix")
  counts_atac <- t(counts_atac_feature_by_cell) # THIS IS THE NECESSARY TRANSPOSE
  message("ATAC counts matrix transposed from feature-by-cell to cell-by-feature (", dim(counts_atac)[1], " cells x ", dim(counts_atac)[2], " features/peaks).")

  # Parse peak names into GRanges object (rowData refers to features)
  peaks_gr <- tryCatch({
    # Check if rowData contains 'chrom', 'start', 'end' explicitly
    rd_df <- as.data.frame(rowData(sce)) # Convert DataFrame to data.frame for easier column access
    if (!is.null(rd_df) && "chrom" %in% colnames(rd_df) && "start" %in% colnames(rd_df) && "end" %in% colnames(rd_df)) {
      GRanges(seqnames = rd_df$chrom,
              ranges = IRanges(start = as.integer(rd_df$start), end = as.integer(rd_df$end)),
              seqinfo = Seqinfo(genome = genome))
    } else {
      # Fallback to parsing SCE rownames (which are feature names)
      Signac::StringToGRanges(rownames(sce), sep = c(":", "-"))
    }
  }, error = function(e) {
    stop(paste("Could not parse peak genomic ranges. Ensure rowData contains 'chrom', 'start', 'end' or original SCE rownames are in 'chr:start-end' format.", e$message))
  })

  # Extract spatial coordinates from `colData(sce)` (cell metadata)
  if (!all(c(x_col, y_col) %in% colnames(colData(sce)))) {
    stop(paste0("Spatial coordinate columns '", x_col, "' and '", y_col, "' not found in AnnData colData. Available columns: ", paste(colnames(colData(sce)), collapse = ", ")))
  }
  spatial_coords <- as.data.frame(colData(sce)[, c(x_col, y_col), drop = FALSE])
  colnames(spatial_coords) <- c("x", "y") # Standardize column names
  rownames(spatial_coords) <- colnames(sce) # Original cell names from SCE

  # Ensure all data aligns by cell names (rownames of counts_atac are cells, rownames of spatial_coords are cells)
  common_cells_final <- intersect(rownames(counts_atac), rownames(spatial_coords))
  if (length(common_cells_final) == 0) {
    stop("No common cells found between ATAC counts and spatial coordinates. Check data integrity and cell ID consistency.")
  }

  counts_atac <- counts_atac[common_cells_final, , drop = FALSE]
  spatial_coords <- spatial_coords[common_cells_final, , drop = FALSE]

  message("AnnData ATAC-seq data loaded. Final dimensions (cells x peaks): ", paste(dim(counts_atac), collapse = "x"), ".")
  return(list(
    counts_matrix = counts_atac, # Now correctly cells x peaks
    peaks_gr = peaks_gr,
    spatial_coords = spatial_coords,
    cell_names = rownames(counts_atac)
  ))
}


#' Mirror Spatial Coordinates
#'
#' \code{mirrorSpatialCoordinates()} mirrors the specified spatial coordinate axis.
#' This is useful for aligning spatial datasets that might be imaged in different orientations.
#'
#' @param spatialCoords A data.frame or matrix with columns for spatial x and y coordinates.
#' @param mirror_axis Character string, either "x" or "y", indicating which axis to mirror.
#' @return A data.frame or matrix with the specified axis mirrored.
#' @export
mirrorSpatialCoordinates <- function(spatialCoords, mirror_axis = "y") {
  if (!(mirror_axis %in% c("x", "y"))) {
    stop("mirror_axis must be 'x' or 'y'.")
  }

  mirrored_coords <- spatialCoords
  if (mirror_axis == "y") {
    # Mirror y-axis: new_y = max_y - old_y
    max_y <- max(spatialCoords$y, na.rm = TRUE)
    mirrored_coords$y <- max_y - spatialCoords$y
    message("Spatial y-coordinates mirrored.")
  } else { # mirror_axis == "x"
    # Mirror x-axis: new_x = max_x - old_x
    max_x <- max(spatialCoords$x, na.rm = TRUE)
    mirrored_coords$x <- max_x - spatialCoords$x
    message("Spatial x-coordinates mirrored.")
  }
  return(mirrored_coords)
}

#' Apply Spatial Filters (Bounding Box and Proximity)
#'
#' \code{applySpatialFilters()} applies bounding box and proximity filters to spatial data.
#'
#' @param spatial_rna_meta A data.frame of metadata for spatial RNA-seq, including spatial x/y columns.
#' @param spatial_atac_meta A data.frame of metadata for spatial ATAC-seq, including spatial x/y columns.
#' @param x_col_rna Name of the X-coordinate column in `spatial_rna_meta`.
#' @param y_col_rna Name of the Y-coordinate column in `spatial_rna_meta`.
#' @param x_col_atac Name of the X-coordinate column in `spatial_atac_meta`.
#' @param y_col_atac Name of the Y-coordinate column in `spatial_atac_meta`.
#' @param x_min_rna Numeric, minimum X-coordinate for RNA bounding box filter.
#' @param x_max_rna Numeric, maximum X-coordinate for RNA bounding box filter.
#' @param x_min_atac Numeric, minimum X-coordinate for ATAC bounding box filter.
#' @param x_max_atac Numeric, maximum X-coordinate for ATAC bounding box filter.
#' @param y_min_rna Numeric, minimum Y-coordinate for RNA bounding box filter.
#' @param y_max_rna Numeric, maximum Y-coordinate for RNA bounding box filter.
#' @param y_min_atac Numeric, minimum Y-coordinate for ATAC bounding box filter.
#' @param y_max_atac Numeric, maximum Y-coordinate for ATAC bounding box filter.
#' @param proximity_radius Numeric, radius for proximity filter.
#' @param min_neighbors_rna Integer, minimum neighbors for RNA proximity filter.
#' @param min_neighbors_atac Integer, minimum neighbors for ATAC proximity filter.
#' @return A list containing `filtered_rna_meta` and `filtered_atac_meta` (data.frames).
#' @export
applySpatialFilters <- function(
    spatial_rna_meta, spatial_atac_meta,
    x_col_rna = "x", y_col_rna = "y",
    x_col_atac = "x", y_col_atac = "y",
    x_min_rna = -Inf, x_max_rna = Inf,
    x_min_atac = -Inf, x_max_atac = Inf,
    y_min_rna = -Inf, y_max_rna = Inf,
    y_min_atac = -Inf, y_max_atac = Inf,
    proximity_radius = 50,
    min_neighbors_rna = 10, min_neighbors_atac = 5
) {
  message("Applying spatial bounding box filter...")
  filtered_rna_meta <- spatial_rna_meta %>%
    filter(!is.na(.data[[x_col_rna]]) & !is.na(.data[[y_col_rna]]) &
             .data[[x_col_rna]] >= x_min_rna & .data[[x_col_rna]] <= x_max_rna &
             .data[[y_col_rna]] >= y_min_rna & .data[[y_col_rna]] <= y_max_rna)

  filtered_atac_meta <- spatial_atac_meta %>%
    filter(!is.na(.data[[x_col_atac]]) & !is.na(.data[[y_col_atac]]) &
             .data[[x_col_atac]] >= x_min_atac & .data[[x_col_atac]] <= x_max_atac &
             .data[[y_col_atac]] >= y_min_atac & .data[[y_col_atac]] <= y_max_atac)

  message("Applying spatial proximity filter...")
  proximity_filter_single <- function(meta, xcol, ycol, radius, min_neighbors) {
    coords <- as.matrix(meta[, c(xcol, ycol)])
    n <- nrow(coords)
    if (n == 0) return(meta[0,,drop=FALSE])

    if (requireNamespace("RANN", quietly = TRUE)) {
      nn <- RANN::nn2(coords, coords, k = min(min_neighbors + 1, n, na.rm = TRUE))$nn.dists # Use RANN, k must be <= n
      within <- rowSums(nn <= radius) - 1 # Exclude self
      keep <- within >= min_neighbors
    } else {
      warning("RANN not installed, using slower proximity filter. Install with install.packages('RANN') for faster filtering.")
      keep <- logical(n)
      r2 <- radius * radius
      for (i in seq_len(n)) {
        dx <- coords[,1] - coords[i,1]
        dy <- coords[,2] - coords[i,2]
        d2 <- dx*dx + dy*dy
        keep[i] <- sum(d2 <= r2) - 1 >= min_neighbors
      }
    }
    meta[keep, , drop = FALSE]
  }

  filtered_rna_meta <- proximity_filter_single(filtered_rna_meta, x_col_rna, y_col_rna, proximity_radius, min_neighbors_rna)
  filtered_atac_meta <- proximity_filter_single(filtered_atac_meta, x_col_atac, y_col_atac, proximity_radius, min_neighbors_atac)

  message("Spatial filters applied. RNA cells remaining: ", nrow(filtered_rna_meta), "; ATAC cells remaining: ", nrow(filtered_atac_meta), ".")
  return(list(filtered_rna_meta = filtered_rna_meta, filtered_atac_meta = filtered_atac_meta))
}


#' Create ATAC-seq Gene Activity Matrix by Aggregating Peak Counts
#'
#' \code{createPeakBasedGeneActivityMatrix()} calculates gene activity scores (pseudo-expression)
#' from a cell-by-peak ATAC count matrix and a `GRanges` object of peaks. It identifies peaks
#' overlapping with extended gene regions (promoter + gene body) and sums their counts per cell.
#' This bypasses the need for raw fragment files, adapting to existing peak matrices.
#'
#' @param atac_counts_matrix Sparse matrix of ATAC-seq peak counts (cells x peaks).
#' @param atac_peaks_gr `GRanges` object of ATAC-seq peaks, matching columns of `atac_counts_matrix`.
#' @param genome Genome build (e.g., "hg38") for gene annotations.
#' @param extend_upstream Integer, distance upstream of TSS to extend gene regions.
#' @param extend_downstream Integer, distance downstream of gene end (TTS) to extend gene regions.
#' @return A sparse matrix of cell-by-gene activity scores.
#' @export
createPeakBasedGeneActivityMatrix <- function(atac_counts_matrix, atac_peaks_gr, genome = "hg38", extend_upstream = 2000, extend_downstream = 0) {
  if (!is(atac_counts_matrix, "Matrix") || !is(atac_peaks_gr, "GRanges")) {
    stop("`atac_counts_matrix` must be a sparse matrix and `atac_peaks_gr` must be a GRanges object.")
  }
  if (ncol(atac_counts_matrix) != length(atac_peaks_gr)) {
    stop("Number of columns in `atac_counts_matrix` (", ncol(atac_counts_matrix), ") must match the length of `atac_peaks_gr` (", length(atac_peaks_gr), ").")
  }
  if (nrow(atac_counts_matrix) == 0) {
    message("No cells in ATAC counts matrix. Returning empty gene activity matrix.")
    return(Matrix::Matrix(0, 0, 0, sparse = TRUE))
  }

  message("Creating ATAC-seq gene activity matrix by aggregating peak counts...")

  # 1. Get gene annotations (promoters + gene bodies) from CellWalkR::getRegions
  gene_regions <- CellWalkR::getRegions(geneBody = TRUE, genome = genome, names = "Entrez")
  if (length(gene_regions) == 0) {
    stop(paste0("Could not retrieve gene regions for genome '", genome, "' using CellWalkR::getRegions. Check genome spelling or implementation."))
  }

  # Ensure seqlevels styles are compatible
  seqlevelsStyle(gene_regions) <- seqlevelsStyle(atac_peaks_gr)

  # 2. Extend gene regions
  # trim() is important after extension to handle chromosome boundaries
  extended_gene_regions <- GenomicRanges::trim(
    GenomicRanges::resize(
      x = gene_regions,
      width = GenomicRanges::width(gene_regions) + extend_upstream + extend_downstream,
      fix = "start" # Extend from TSS/start of original region
    )
  )

  # Ensure gene_id is present and valid for naming the output matrix
  if (is.null(extended_gene_regions$gene_id) || any(is.na(extended_gene_regions$gene_id))) {
    stop("`gene_id` metadata column is missing or contains NA values in `extended_gene_regions`.")
  }

  # 3. Map peaks to extended gene regions
  # Find overlaps between peaks and extended gene regions
  overlaps <- GenomicRanges::findOverlaps(atac_peaks_gr, extended_gene_regions)

  if (length(overlaps) == 0) {
    message("No peaks overlap with gene regions. Returning empty gene activity matrix.")
    # Return empty matrix with correct number of cells and gene names
    gene_names <- unique(extended_gene_regions$gene_id)
    return(Matrix::Matrix(0, nrow(atac_counts_matrix), length(gene_names), sparse = TRUE, dimnames = list(rownames(atac_counts_matrix), gene_names)))
  }

  # 4. Aggregate peak counts per gene
  # `queryHits(overlaps)` are indices of peaks (columns in `atac_counts_matrix`)
  # `subjectHits(overlaps)` are indices of gene_regions (rows in `extended_gene_regions`)

  # Create a sparse matrix to map peaks to genes (features x genes)
  peak_to_gene_map <- Matrix::sparseMatrix(
    i = queryHits(overlaps), # Peak index
    j = subjectHits(overlaps), # Gene index
    x = 1, # A 1 indicates this peak overlaps this gene
    dims = c(length(atac_peaks_gr), length(extended_gene_regions))
  )
  colnames(peak_to_gene_map) <- extended_gene_regions$gene_id
  rownames(peak_to_gene_map) <- names(atac_peaks_gr) # Optional, for clarity

  # Perform matrix multiplication: (cells x peaks) %*% (peaks x genes) = (cells x genes)
  # This sums the counts of overlapping peaks for each gene, for each cell.
  gene_activity_matrix <- atac_counts_matrix %*% peak_to_gene_map

  # Remove duplicate gene_id columns if any genes have multiple entries in extended_gene_regions
  # This can happen if getRegions returns multiple ranges for the same gene_id (e.g., separate promoter and gene body regions).
  # We need to sum these.
  unique_gene_ids <- unique(colnames(gene_activity_matrix))
  if (length(unique_gene_ids) < ncol(gene_activity_matrix)) {
    message("Aggregating counts for duplicate gene IDs in gene activity matrix using vectorized sparse matrix multiplication.")

    # Get the original (potentially duplicated) column names
    original_colnames <- colnames(gene_activity_matrix)

    # Create a mapping from each original column name to its position in `unique_gene_ids`
    # e.g., if original_colnames = c("A", "B", "A"), unique_gene_ids = c("A", "B")
    # then col_map_idx = c(1, 2, 1)
    col_map_idx <- match(original_colnames, unique_gene_ids)

    # Construct a sparse "aggregation" matrix.
    # This matrix will have:
    # - Rows equal to the number of columns in `gene_activity_matrix` (original genes with duplicates).
    # - Columns equal to the number of `unique_gene_ids`.
    # An entry (i, j) will be 1 if the i-th original gene column should be summed into the j-th unique gene column.
    agg_map_matrix <- Matrix::sparseMatrix(
      i = 1:ncol(gene_activity_matrix), # Row index for each original gene column
      j = col_map_idx,                  # Column index for its corresponding unique gene
      x = 1,                            # Value to sum (1, indicating inclusion)
      dims = c(ncol(gene_activity_matrix), length(unique_gene_ids))
    )

    # Perform the aggregation using matrix multiplication.
    # (cells x original_genes) %*% (original_genes x unique_genes) = (cells x unique_genes)
    gene_activity_matrix <- gene_activity_matrix %*% agg_map_matrix

    # Assign the correct unique gene IDs as column names.
    # (Row names are automatically preserved from the original `gene_activity_matrix` during multiplication).
    colnames(gene_activity_matrix) <- unique_gene_ids
  }

  message("ATAC-seq gene activity matrix created. Dimensions (cells x genes): ", paste(dim(gene_activity_matrix), collapse = "x"))
  return(gene_activity_matrix)
}


#' Impute ATAC-seq Gene Activity Data onto RNA-seq Spots Based on Spatial Proximity
#'
#' \code{imputeATACtoRNAbySpatial()} imputes ATAC-seq gene activity scores onto RNA-seq spots
#' by finding the `k` nearest ATAC cells for each RNA spot and averaging their gene activity scores.
#' This creates an ATAC gene activity matrix where rows correspond to RNA spots.
#'
#' @param rna_spatial_coords Data.frame of RNA spatial coordinates (cell_id, x, y).
#' @param atac_spatial_coords Data.frame of ATAC spatial coordinates (cell_id, x, y).
#' @param atac_gene_activity_matrix Sparse matrix of ATAC-seq gene activity scores (cells x genes),
#'                           with rownames matching `rownames(atac_spatial_coords)`.
#' @param k_impute_spatial Integer, the number of nearest ATAC cells to consider for imputation.
#' @return A sparse matrix of imputed ATAC gene activity scores (RNA_cells x genes).
#' @export
imputeATACtoRNAbySpatial <- function(rna_spatial_coords, atac_spatial_coords, atac_gene_activity_matrix, k_impute_spatial = 5, weight_by = c("uniform", "inverse_dist", "gaussian"), gaussian_sigma = NULL, eps = 1e-6, return_weights = FALSE) {
  weight_by <- match.arg(weight_by)
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for spatial imputation. Please install it with install.packages('RANN')")
  }
  if (nrow(rna_spatial_coords) == 0 || nrow(atac_spatial_coords) == 0 || nrow(atac_gene_activity_matrix) == 0) {
    message("No cells in RNA or ATAC spatial coordinates or ATAC gene activity. Skipping imputation.")
    empty <- Matrix::Matrix(0, nrow(rna_spatial_coords), ncol(atac_gene_activity_matrix), sparse = TRUE, dimnames = list(rownames(rna_spatial_coords), colnames(atac_gene_activity_matrix)))
    if (return_weights) return(list(imputed = empty, weights = NULL))
    return(empty)
  }

  k_impute_spatial <- min(k_impute_spatial, nrow(atac_spatial_coords))
  if (k_impute_spatial == 0) {
    warning("k_impute_spatial became 0, no ATAC cells available for imputation. Returning empty matrix.")
    empty <- Matrix::Matrix(0, nrow(rna_spatial_coords), ncol(atac_gene_activity_matrix), sparse = TRUE, dimnames = list(rownames(rna_spatial_coords), colnames(atac_gene_activity_matrix)))
    if (return_weights) return(list(imputed = empty, weights = NULL))
    return(empty)
  }

  nn_result <- RANN::nn2(
    query = as.matrix(rna_spatial_coords[, c("x", "y")]),
    data = as.matrix(atac_spatial_coords[, c("x", "y")]),
    k = k_impute_spatial
  )

  # Build weight values per neighbor depending on weighting scheme
  distances_mat <- nn_result$nn.dists
  if (weight_by == "uniform") {
    weight_vals_mat <- matrix(1 / k_impute_spatial, nrow = nrow(distances_mat), ncol = ncol(distances_mat))
  } else if (weight_by == "inverse_dist") {
    weight_vals_mat <- 1 / (distances_mat + eps)
    # Normalize row-wise to sum to 1
    row_sums <- rowSums(weight_vals_mat)
    row_sums[row_sums == 0] <- 1
    weight_vals_mat <- weight_vals_mat / row_sums
  } else { # gaussian
    if (is.null(gaussian_sigma)) {
      sigma <- median(distances_mat[is.finite(distances_mat) & distances_mat > 0], na.rm = TRUE)
      if (!is.finite(sigma) || sigma <= 0) sigma <- 1
    } else {
      sigma <- gaussian_sigma
    }
    weight_vals_mat <- exp(- (distances_mat^2) / (2 * sigma^2))
    # normalize rows
    row_sums <- rowSums(weight_vals_mat)
    row_sums[row_sums == 0] <- 1
    weight_vals_mat <- weight_vals_mat / row_sums
  }

  # Construct sparse weight matrix W (N_rna x N_atac)
  row_indices_W <- rep(seq_len(nrow(rna_spatial_coords)), each = k_impute_spatial)
  col_indices_W <- as.vector(nn_result$nn.idx)
  weight_values <- as.vector(t(weight_vals_mat)) # ensure same order as flattened indices

  # Filter invalid indices if any
  valid_mask <- !is.na(col_indices_W) & col_indices_W > 0 & col_indices_W <= nrow(atac_spatial_coords)
  row_indices_W <- row_indices_W[valid_mask]
  col_indices_W <- col_indices_W[valid_mask]
  weight_values <- weight_values[valid_mask]

  W <- Matrix::sparseMatrix(
    i = row_indices_W,
    j = col_indices_W,
    x = weight_values,
    dims = c(nrow(rna_spatial_coords), nrow(atac_spatial_coords)),
    dimnames = list(rownames(rna_spatial_coords), rownames(atac_spatial_coords))
  )

  imputed_atac_matrix <- W %*% atac_gene_activity_matrix
  message("ATAC-seq gene activity data imputation complete. Imputed matrix dims: ", paste(dim(imputed_atac_matrix), collapse = "x"))
  if (return_weights) return(list(imputed = imputed_atac_matrix, weights = W))
  return(imputed_atac_matrix)
}


#' Simulate pseudo-bulk RNA by spatial neighbor aggregation
#'
#' simulatePseudoBulkRNA() creates pseudo-bulk expression matrices for RNA spots by
#' aggregating counts from each cell's spatial neighbors. This can be used to mimic
#' the multi-cell sampling / smoothing effect of ATAC measurements (2–10 cells).
#'
#' For each cell, its S nearest spatial neighbors (including itself, if returned by
#' the nearest-neighbor search) are identified and their counts are summed or
#' averaged to form a pseudo-bulk profile. A separate pseudo-bulk matrix is
#' generated for each requested neighborhood size.
#'
#' @param seurat_obj Seurat object containing RNA counts (assay = "RNA", layer = "counts").
#' @param rna_spatial_coords data.frame with rownames matching RNA cell IDs and
#'   columns \code{x}, \code{y} giving spatial coordinates.
#' @param sample_sizes Integer vector of neighborhood sizes (including self) to aggregate,
#'   e.g. \code{c(2, 3, 5, 10)}.
#' @param method Aggregation method for neighbors: \code{"sum"} (default) to sum counts,
#'   or \code{"mean"} to average counts within each neighborhood.
#'
#' @return A named list of matrices (genes x pseudoCells), one for each neighborhood size,
#'   with list elements named like \code{"S2"}, \code{"S3"}, etc. Columns are in the same
#' order as the input cells (i.e. one pseudo-cell per original cell).
#'
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix sparseMatrix colSums Diagonal
simulatePseudoBulkRNA <- function(seurat_obj,
                                  rna_spatial_coords,
                                  sample_sizes = c(2, 3, 5, 10),
                                  method = c("sum", "mean")) {
  method <- match.arg(method)

  # Extract counts matrix: genes x cells
  mat <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts"))
  cell_ids <- colnames(mat)

  # Check and align spatial coordinates
  if (!all(cell_ids %in% rownames(rna_spatial_coords))) {
    stop("Spatial coords must contain all RNA cell IDs as rownames.")
  }
  coords <- as.matrix(rna_spatial_coords[cell_ids, c("x", "y")])
  if (any(is.na(coords))) {
    stop("NA spatial coords detected for some RNA cells.")
  }

  # RANN::nn2 requires k <= n
  n_cells <- nrow(coords)
  max_k <- min(max(sample_sizes), n_cells)

  nn <- RANN::nn2(data = coords, query = coords, k = max_k)

  out <- list()

  for (S in sample_sizes) {
    S_eff <- min(S, n_cells)
    idx_mat <- nn$nn.idx[, 1:S_eff, drop = FALSE]  # Ncells x S_eff (indices into cell_ids)

    # Build a (cells x pseudoCells) sparse weight matrix Wp
    # Wp[j, i] = 1 if cell j is a neighbor of pseudo-cell i
    row_indices <- as.vector(idx_mat)                    # neighbor cell indices (row positions)
    col_indices <- rep(seq_len(n_cells), each = S_eff)   # pseudo-cell indices (column positions)
    vals <- rep(1, length(row_indices))

    Wp <- Matrix::sparseMatrix(i = row_indices,
                               j = col_indices,
                               x = vals,
                               dims = c(n_cells, n_cells))

    if (method == "mean") {
      # Normalize each column so weights sum to 1 (average over neighbors)
      col_sums <- Matrix::colSums(Wp)
      inv_col_sums <- 1 / ifelse(col_sums == 0, 1, col_sums)
      D <- Matrix::Diagonal(x = inv_col_sums)
      Wp <- Wp %*% D
    }

    # mat: genes x cells; Wp: cells x pseudoCells -> pseudos: genes x pseudoCells
    pseudos <- mat %*% Wp

    colnames(pseudos) <- cell_ids
    rownames(pseudos) <- rownames(mat)

    out[[paste0("S", S)]] <- pseudos
  }

  return(out)
}


#' Compare neighbor rank ranges vs random cells (median spatial distance distributions)
#'
#' assessRangeVsRandom() computes per-cell median spatial distances for two neighbor
#' rank ranges taken from a symmetrized feature graph (e.g., top 1-10 and 11-20),
#' and for random cell samples. Returns a density plot with all three distributions.
#'
#' @param feature_graph sparse matrix (cells x cells) with rownames = cell IDs
#' @param feature_matrix_for_similarity features x cells (only used to validate alignment)
#' @param spatial_coords_df data.frame with x,y and rownames = cell IDs
#' @param k_full integer maximum neighbor rank to consider (default 20)
#' @param range1 integer vector length 2 (e.g., c(1,10))
#' @param range2 integer vector length 2 (e.g., c(11,20))
#' @param n_random_per_cell integer, how many random cells to sample per cell for random baseline
#' @param modality_name label for plot titles
#' @param n_workers parallel workers (currently not used inside function but kept for signature)
assessRangeVsRandom <- function(feature_graph, feature_matrix_for_similarity, spatial_coords_df,
                                k_full = 20, range1 = c(1,10), range2 = c(11,20),
                                n_random_per_cell = 20, modality_name = "Modality", n_workers = 4) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install Matrix.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table.")
  if (!requireNamespace("RANN", quietly = TRUE)) stop("Install RANN.")
  if (is.null(rownames(feature_graph)) || is.null(colnames(feature_matrix_for_similarity)) || is.null(rownames(spatial_coords_df))) {
    stop("feature_graph rownames, feature_matrix_for_similarity colnames and spatial_coords_df rownames must be set and matching.")
  }
  common_cells <- Reduce(intersect, list(rownames(feature_graph), colnames(feature_matrix_for_similarity), rownames(spatial_coords_df)))
  if (length(common_cells) == 0) stop("No common cells across inputs.")
  fg <- feature_graph[common_cells, common_cells, drop = FALSE]
  coords <- spatial_coords_df[common_cells, c("x","y"), drop = FALSE]
  cell_ids <- common_cells
  num_cells <- length(cell_ids)
  # symmetrize
  fg_sym <- (fg + Matrix::t(fg)) / 2
  diag(fg_sym) <- 0
  # neighbor ranking from feature graph
  s <- Matrix::summary(fg_sym)
  neighbors_rank_list <- vector("list", num_cells)
  names(neighbors_rank_list) <- cell_ids
  if (nrow(s) > 0) {
    dt <- data.table::as.data.table(s)
    data.table::setnames(dt, c("i","j","x"))
    dt[, ':=' (i = as.integer(i), j = as.integer(j))]
    dt_sorted <- dt[order(i, -x)]
    grouped <- split(dt_sorted, by = "i", keep.by = FALSE)
    for (i_idx in names(grouped)) {
      ii <- as.integer(i_idx)
      neighbors_rank_list[[ii]] <- cell_ids[grouped[[i_idx]]$j]
    }
  } else {
    neighbors_rank_list <- lapply(neighbors_rank_list, function(x) character(0))
  }
  coords_mat <- as.matrix(coords)
  # helper to compute median spatial distance for given neighbor rank range
  get_median_spatial_for_range <- function(range_vec) {
    res <- numeric(num_cells)
    for (i in seq_len(num_cells)) {
      neighs <- neighbors_rank_list[[i]]
      if (length(neighs) == 0) { res[i] <- NA; next }
      sel_idx <- seq(range_vec[1], range_vec[2])
      sel <- neighs[sel_idx[sel_idx <= length(neighs)]]
      sel <- sel[!is.na(sel) & sel %in% cell_ids]
      if (length(sel) == 0) { res[i] <- NA; next }
      d <- sqrt((coords_mat[i,1] - coords_mat[sel,1])^2 + (coords_mat[i,2] - coords_mat[sel,2])^2)
      res[i] <- median(d, na.rm = TRUE)
    }
    res
  }
  m1 <- get_median_spatial_for_range(range1)
  m2 <- get_median_spatial_for_range(range2)
  # random per-cell medians
  set.seed(123)
  mrand <- numeric(num_cells)
  for (i in seq_len(num_cells)) {
    pool <- setdiff(cell_ids, cell_ids[i])
    if (length(pool) == 0) { mrand[i] <- NA; next }
    sel <- sample(pool, min(n_random_per_cell, length(pool)))
    d <- sqrt((coords_mat[i,1] - coords_mat[sel,1])^2 + (coords_mat[i,2] - coords_mat[sel,2])^2)
    mrand[i] <- median(d, na.rm = TRUE)
  }
  df_plot <- data.frame(
    cell_id = rep(cell_ids, 3),
    median_spatial = c(m1, m2, mrand),
    group = rep(c(paste0("NN_", range1[1], "_", range1[2]), paste0("NN_", range2[1], "_", range2[2]), "Random"), each = num_cells),
    stringsAsFactors = FALSE
  )
  df_plot$group <- factor(df_plot$group, levels = unique(df_plot$group))
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = median_spatial, color = group, fill = group)) +
    ggplot2::geom_density(alpha = 0.3, na.rm = TRUE) +
    ggplot2::labs(title = paste0(modality_name, ": Median Spatial Distance — range vs range vs random"), x = "Median Spatial Distance", y = "Density") +
    ggplot2::theme_minimal()
  list(df = df_plot, plot = p)
}


#' Preprocess Spatial Multi-omics Data for CellWalker2
#'
#' \code{preprocessSpatialMultiomics()} is a wrapper to import, process, and filter
#' spatial scRNA-seq and spatial scATAC-seq data for CellWalker2. It returns the
#' individually filtered data for each modality, without forcing common cell IDs,
#' allowing for flexible downstream integration.
#'
#' @param rna_h5_path Path to the Xenium cell_feature_matrix.h5 file.
#' @param cells_csv_path Path to the Xenium cells.csv.gz file.
#' @param atac_h5ad_path Path to the AnnData .h5ad file for ATAC-seq data.
#' @param atac_x_col Name of the X-coordinate column in ATAC AnnData `colData` (e.g., "xcor").
#' @param atac_y_col Name of the Y-coordinate column in ATAC AnnData `colData` (e.g., "ycor").
#' @param genome Genome build (e.g., "hg38") for ATAC peaks.
#' @param mirror_rna_y Boolean, if TRUE, mirrors the y-coordinates of the RNA spatial data.
#' @param rna_proj_name Project name for the RNA Seurat object.
#' @param rna_min_cells_initial For Seurat's `CreateSeuratObject`, include features (genes) detected in at least this many cells.
#' @param rna_min_features_initial For Seurat's `CreateSeuratObject`, include cells with at least this many features (genes).
#' @param rna_nfeatures For Seurat's `FindVariableFeatures`, number of features to select as top variable features.
#' @param filter_params A list of parameters for `applySpatialFilters`:
#'   \itemize{
#'     \item `x_min_rna`, `x_max_rna`, `y_min_rna`, `y_max_rna`: Bounding box coordinates for RNA.
#'     \item `x_min_atac`, `x_max_atac`, `y_min_atac`, `y_max_atac`: Bounding box coordinates for ATAC.
#'     \item `proximity_radius`: Radius for proximity filter.
#'     \item `min_neighbors_rna`, `min_neighbors_atac`: Min neighbors for proximity filter.
#'   }
#' @return A list containing preprocessed data for RNA and ATAC (individually filtered):
#'   \itemize{
#'     \item `seurat_obj_rna`: Preprocessed Seurat object for spatial RNA-seq.
#'     \item `rna_spatial_coords`: Data.frame of RNA spatial coordinates (cell_id, x, y).
#'     \item `atac_counts_matrix`: Sparse matrix of ATAC-seq peak counts (cells x peaks).
#'     \item `atac_peaks_gr`: `GRanges` object of ATAC-seq peaks.
#'     \item `atac_spatial_coords`: Data.frame of ATAC spatial coordinates (cell_id, x, y).
#'   }
#' @export
preprocessSpatialMultiomics <- function(
    rna_h5_path, cells_csv_path,
    atac_h5ad_path, atac_x_col = "xcor", atac_y_col = "ycor", genome = "hg38",
    mirror_rna_y = TRUE,
    rna_proj_name = "Xenium_RNA",
    rna_min_cells_initial = 3,
    rna_min_features_initial = 100,
    rna_nfeatures = 3000,
    filter_params = list(
      x_min_rna = -Inf, x_max_rna = Inf, x_min_atac = -Inf, x_max_atac = Inf, y_min_rna = -Inf, y_max_rna = Inf, y_min_atac = -Inf, y_max_atac = Inf,
      proximity_radius = 50, min_neighbors_rna = 10, min_neighbors_atac = 5
    )
) {
  message("Starting `preprocessSpatialMultiomics` pipeline...")
  # 1. Read Spatial RNA-seq Data
  seurat_rna <- readXeniumData(
    rna_h5_path, cells_csv_path,
    project_name = rna_proj_name,
    min_cells = rna_min_cells_initial,
    min_features = rna_min_features_initial
  )

  seurat_rna@meta.data$x <- as.numeric(seurat_rna@meta.data$x_centroid)
  seurat_rna@meta.data$y <- as.numeric(seurat_rna@meta.data$y_centroid)
  rna_spatial_coords <- seurat_rna@meta.data[, c("x", "y"), drop = FALSE]
  rownames(rna_spatial_coords) <- colnames(seurat_rna)

  # 2. Read Spatial ATAC-seq Data
  atac_data <- readAnnDataATAC(atac_h5ad_path, x_col = atac_x_col, y_col = atac_y_col, genome = genome)
  atac_counts_matrix <- atac_data$counts_matrix
  atac_peaks_gr <- atac_data$peaks_gr
  atac_spatial_coords <- atac_data$spatial_coords

  # 3. Apply Spatial Filters (Bounding Box and Proximity)
  message("Applying spatial filters (bounding box and proximity)...")
  filtered_coords_list <- applySpatialFilters(
    spatial_rna_meta = rna_spatial_coords,
    spatial_atac_meta = atac_spatial_coords,
    x_col_rna = "x", y_col_rna = "y",
    x_col_atac = "x", y_col_atac = "y",
    x_min_rna = filter_params$x_min_rna, x_max_rna = filter_params$x_max_rna,
    x_min_atac = filter_params$x_min_atac, x_max_atac = filter_params$x_max_atac,
    y_min_rna = filter_params$y_min_rna, y_max_rna = filter_params$y_max_rna,
    y_min_atac = filter_params$y_min_atac, y_max_atac = filter_params$y_max_atac,
    proximity_radius = filter_params$proximity_radius,
    min_neighbors_rna = filter_params$min_neighbors_rna,
    min_neighbors_atac = filter_params$min_neighbors_atac
  )

  # Update original data with filtered cells
  filtered_rna_cells <- rownames(filtered_coords_list$filtered_rna_meta)
  filtered_atac_cells <- rownames(filtered_coords_list$filtered_atac_meta)

  seurat_rna <- subset(seurat_rna, cells = filtered_rna_cells)
  rna_spatial_coords <- rna_spatial_coords[filtered_rna_cells, , drop = FALSE]

  atac_counts_matrix <- atac_counts_matrix[filtered_atac_cells, , drop = FALSE]
  atac_spatial_coords <- atac_spatial_coords[filtered_atac_cells, , drop = FALSE]

  message("After spatial filtering: RNA cells = ", ncol(seurat_rna), ", ATAC cells = ", nrow(atac_counts_matrix), ".")

  # 4. Mirror RNA spatial coordinates if requested (after filtering, so max_y is correct for filtered data)
  if (mirror_rna_y) {
    rna_spatial_coords <- mirrorSpatialCoordinates(rna_spatial_coords, mirror_axis = "y")
    seurat_rna@meta.data$y <- rna_spatial_coords[rownames(seurat_rna@meta.data), "y"]
    message("RNA Y-coordinates mirrored.")
  }

  # 5. Seurat standard preprocessing for RNA
  message("Performing Seurat standard preprocessing for RNA data...")
  DefaultAssay(seurat_rna) <- "RNA"

  seurat_rna <- NormalizeData(seurat_rna, verbose = FALSE)
  seurat_rna <- FindVariableFeatures(seurat_rna, nfeatures = rna_nfeatures, verbose = FALSE)

  num_var_features <- length(VariableFeatures(seurat_rna))
  message("  Found ", num_var_features, " variable features.")

  if (num_var_features == 0) {
    warning("No variable features found after FindVariableFeatures(). RunPCA may fail.")
  } else {
    seurat_rna <- ScaleData(seurat_rna, verbose = FALSE)
    seurat_rna <- RunPCA(seurat_rna, verbose = FALSE)
  }

  message("`preprocessSpatialMultiomics` pipeline complete. Returning filtered individual modality data.")
  return(list(
    seurat_obj_rna = seurat_rna,
    rna_spatial_coords = rna_spatial_coords,
    atac_counts_matrix = atac_counts_matrix,
    atac_peaks_gr = atac_peaks_gr,
    atac_spatial_coords = atac_spatial_coords
  ))
}


#' Create RNA-based Cell-to-Cell Similarity Graph from Seurat Object
#'
#' \code{createRNACellGraph()} computes a cell-to-cell similarity matrix based on
#' PCA embeddings from a preprocessed Seurat object. This function extracts the
#' PCA coordinates and computes a sparse K-nearest neighbors (KNN) graph.
#' Similarity scores are calculated using a robust Gaussian kernel based on
#' the median neighbor distance.
#'
#' @param seurat_obj A Seurat object that has already undergone `NormalizeData`,
#'                   `FindVariableFeatures`, `ScaleData`, and `RunPCA`.
#' @param dims An integer vector specifying which principal components (PCs) to use
#'             for calculating cell-cell similarity. Defaults to 1:30, but will
#'             be capped at the number of available PCs.
#' @param k_neighbors Integer, the number of nearest neighbors to consider for each cell.
#'                    This controls the sparsity of the resulting graph.
#' @return A sparseMatrix of cell-to-cell similarity, with values normalized
#'         between 0 and 1, where 1 indicates perfect similarity. Rownames and
#'         colnames are cell barcodes. Returns an empty (0x0) or 1x1 identity matrix
#'         if 0 or 1 cells remain after processing.
#' @export
createRNACellGraph <- function(seurat_obj, dims = 1:30, k_neighbors = 20) {
  # --- 1. Initial Validation ---
  if (!is(seurat_obj, "Seurat")) {
    stop("Input 'seurat_obj' must be a Seurat object.")
  }

  message("Computing RNA-based cell-cell KNN graph using Gaussian Kernel...")

  if (is.null(seurat_obj@reductions$pca)) {
    stop("Error: PCA has not been run on the provided Seurat object. Please ensure RunPCA is called.")
  }

  pca_embeddings_full <- seurat_obj@reductions$pca@cell.embeddings
  max_pcs_available <- ncol(pca_embeddings_full)

  if (max_pcs_available < 1) {
    stop("Error: No principal components found in the Seurat object. Check PCA computation or if variable features were found.")
  }

  # --- 2. Dimensionality Handling ---
  effective_dims <- dims[dims <= max_pcs_available]

  if (length(effective_dims) < 2) {
    warning("Insufficient principal components (less than 2) available or selected. Returning NULL graph.")
    return(NULL)
  }

  pca_embeddings <- pca_embeddings_full[, effective_dims, drop = FALSE]

  # Failsafe: Check for non-finite values (NA, NaN, Inf)
  if (any(!is.finite(pca_embeddings))) {
    warning("Non-finite values detected in PCA embeddings. Removing affected cells.")
    initial_nrow <- nrow(pca_embeddings)
    pca_embeddings <- na.omit(pca_embeddings)
    if (nrow(pca_embeddings) < initial_nrow) {
      message("Removed ", initial_nrow - nrow(pca_embeddings), " cells. Remaining: ", nrow(pca_embeddings))
    }
  }

  # Failsafe: Check for zero variance in PCA (collapsed data)
  pca_var <- apply(pca_embeddings, 2, var)
  if (all(pca_var == 0)) {
    stop("CRITICAL: PCA embeddings have zero variance. Check ScaleData/RunPCA steps.")
  }

  num_cells <- nrow(pca_embeddings)

  # --- 3. Handle Small Datasets ---
  if (num_cells <= 1) {
    message("Warning: Too few cells (<= 1) for graph computation.")
    if (num_cells == 0) {
      return(Matrix::Matrix(0, 0, 0, sparse = TRUE, dimnames = list(NULL, NULL)))
    } else {
      return(Matrix::Matrix(1, 1, 1, sparse = TRUE, dimnames = list(rownames(pca_embeddings), rownames(pca_embeddings))))
    }
  }

  # Adjust k_neighbors if it exceeds cell count
  if (k_neighbors >= num_cells) {
    warning("'k_neighbors' (", k_neighbors, ") >= number of cells. Setting to num_cells - 1.")
    k_neighbors <- max(1, num_cells - 1)
  }

  # --- 4. KNN Computation ---
  message("  Computing k-nearest neighbors (k = ", k_neighbors, ") using FNN::get.knn...")
  knn_result <- FNN::get.knn(pca_embeddings, k = k_neighbors)

  if (is.null(knn_result) || is.null(knn_result$nn.index) || is.null(knn_result$nn.dist)) {
    stop("FNN::get.knn returned NULL or empty results.")
  }

  adj_list_i_all <- rep(1:num_cells, each = k_neighbors)
  adj_list_j_raw <- as.vector(knn_result$nn.index)
  distances_raw <- as.vector(knn_result$nn.dist)

  # --- 5. Robust Similarity Logic (Gaussian Kernel) ---
  # We use the median distance as our bandwidth (sigma).
  # This adapts the similarity decay to the local density of the data.
  sigma <- median(distances_raw, na.rm = TRUE)
  if (sigma <= 0) sigma <- max(distances_raw, na.rm = TRUE) # Fallback for low-variance data
  if (sigma <= 0) sigma <- 1 # Final absolute fallback to prevent division by zero

  # Gaussian similarity formula: exp(-d^2 / (2 * sigma^2))
  # fix to max of 1 as that is the theoretical limit (should normally not be bigger but with weird underflows and roundings it can get over it)
  similarity <- pmin(1, exp(- (distances_raw^2) / (2 * sigma^2)))

  # --- 6. Sparse Matrix Construction ---
  # Only keep valid indices and similarities
  valid_entries <- !is.na(adj_list_j_raw) & adj_list_j_raw > 0 & is.finite(similarity)

  rna_cell_graph_init <- Matrix::sparseMatrix(
    i = adj_list_i_all[valid_entries],
    j = adj_list_j_raw[valid_entries],
    x = similarity[valid_entries],
    dims = c(num_cells, num_cells),
    dimnames = list(rownames(pca_embeddings), rownames(pca_embeddings)),
    giveCsparse = TRUE
  )

  # Ensure diagonal is 1
  diag(rna_cell_graph_init) <- 1

  # --- 7. Robust Symmetrization ---
  message("  Symmetrizing graph using max rule...")
  rna_triplets <- summary(rna_cell_graph_init)

  combined_triplets_df <- data.frame(
    i = c(rna_triplets$i, rna_triplets$j),
    j = c(rna_triplets$j, rna_triplets$i),
    x = c(rna_triplets$x, rna_triplets$x)
  )

  if (!requireNamespace("data.table", quietly = TRUE)) {
    message("  WARNING: Package 'data.table' not installed. Using base R (slower).")
    symmetric_edges <- aggregate(x ~ i + j, data = combined_triplets_df, FUN = max)
  } else {
    data.table::setDT(combined_triplets_df)
    symmetric_edges <- combined_triplets_df[, .(x = max(x)), by = .(i, j)]
  }

  rna_cell_graph <- Matrix::sparseMatrix(
    i = symmetric_edges$i,
    j = symmetric_edges$j,
    x = symmetric_edges$x,
    dims = dim(rna_cell_graph_init),
    dimnames = dimnames(rna_cell_graph_init),
    giveCsparse = TRUE
  )

  message("RNA Graph created using Gaussian Kernel. Non-zero entries: ", Matrix::nnzero(rna_cell_graph))
  return(rna_cell_graph)
}


#' Create ATAC-based Cell-to-Cell Similarity Graph from Gene Activity Matrix
#'
#' \code{createATACCellGraph()} computes a cell-to-cell similarity matrix based on
#' ATAC-seq gene activity scores. It performs sparse PCA on the gene activity matrix
#' and then constructs a K-nearest neighbors (KNN) graph from the PCA embeddings,
#' using `FNN::get.knn` with a Gaussian kernel for similarity calculation.
#'
#' @param atac_gene_activity_matrix Sparse matrix of ATAC-seq gene activity scores (cells x genes).
#' @param atac_pca_dims Integer, number of principal components to use for ATAC graph construction.
#' @param atac_k_neighbors Integer, the number of nearest neighbors to consider for each cell
#'                         in the ATAC KNN graph.
#' @return A sparseMatrix of cell-to-cell similarity, with values normalized
#'         between 0 and 1, where 1 indicates perfect similarity. Rownames and
#'         colnames are cell barcodes. Returns an empty (0x0) or 1x1 identity matrix
#'         if 0 or 1 cells remain after processing.
#' @importFrom irlba prcomp_irlba
#' @export
createATACCellGraph <- function(atac_gene_activity_matrix, atac_pca_dims = 30, atac_k_neighbors = 20) {
  message("Computing ATAC-based cell-cell graph from gene activity matrix using PCA+kNN+Gaussian Kernel...")

  if (!is(atac_gene_activity_matrix, "Matrix") && !is.matrix(atac_gene_activity_matrix)) {
    stop("`atac_gene_activity_matrix` must be a matrix or sparse Matrix (cells x genes).")
  }

  num_cells <- nrow(atac_gene_activity_matrix)
  cell_ids <- rownames(atac_gene_activity_matrix)
  if (is.null(cell_ids)) cell_ids <- as.character(seq_len(num_cells))

  if (num_cells <= 1) {
    message("Warning: Too few cells (<= 1) remaining for ATAC-based cell graph computation. Returning empty or identity matrix.")
    if (num_cells == 0) {
      return(list(graph = Matrix::Matrix(0, 0, 0, sparse = TRUE, dimnames = list(NULL, NULL)), pca = NULL))
    } else { # num_cells == 1
      return(list(graph = Matrix::Matrix(1, 1, 1, sparse = TRUE, dimnames = list(cell_ids, cell_ids)), pca = NULL))
    }
  }

  non_zero_genes <- Matrix::colSums(atac_gene_activity_matrix) > 0
  if (sum(non_zero_genes) == 0) {
    warning("ATAC gene activity matrix contains no active genes. Cannot create ATAC feature graph. Returning empty matrix.")
    empty_mat <- Matrix::Matrix(0, num_cells, num_cells, sparse = TRUE, dimnames = list(cell_ids, cell_ids))
    return(list(graph = empty_mat, pca = NULL))
  }

  atac_mat_for_dimred <- atac_gene_activity_matrix[, non_zero_genes, drop = FALSE]
  message("  Performing sparse PCA on ATAC gene activity matrix...")

  num_features_for_pca <- ncol(atac_mat_for_dimred)
  if (num_features_for_pca < 2) {
    warning("Too few features (", num_features_for_pca, ") for ATAC dimensionality reduction. Returning empty graph.")
    empty_mat <- Matrix::Matrix(0, num_cells, num_cells, sparse = TRUE, dimnames = list(cell_ids, cell_ids))
    return(list(graph = empty_mat, pca = NULL))
  }

  effective_atac_pca_dims <- min(atac_pca_dims, num_features_for_pca)
  if (effective_atac_pca_dims < 2 && num_features_for_pca >= 2) effective_atac_pca_dims <- 2
  if (effective_atac_pca_dims < 1) {
    warning("Not enough dimensions for ATAC PCA. Returning empty graph.")
    empty_mat <- Matrix::Matrix(0, num_cells, num_cells, sparse = TRUE, dimnames = list(cell_ids, cell_ids))
    return(list(graph = empty_mat, pca = NULL))
  }

  # Use irlba on sparse input if available, avoiding as.matrix conversion when possible
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for sparse PCA. Install with install.packages('irlba').")
  }

  # irlba::prcomp_irlba works well for sparse matrices in many cases
  atac_pca_res <- tryCatch({
    irlba::prcomp_irlba(atac_mat_for_dimred, n = effective_atac_pca_dims, center = TRUE, scale. = FALSE)
  }, error = function(e) {
    # Fallback: try converting to dense if small enough
    warning("irlba failed on sparse matrix; attempting dense prcomp_irlba conversion as fallback (may be memory intensive): ", e$message)
    dense_mat <- as.matrix(atac_mat_for_dimred)
    irlba::prcomp_irlba(dense_mat, n = effective_atac_pca_dims, center = TRUE, scale. = FALSE)
  })

  if (is.null(atac_pca_res) || is.null(atac_pca_res$x)) {
    stop("ATAC PCA failed to produce embeddings. Aborting ATAC graph construction.")
  }

  pca_embeddings_atac <- atac_pca_res$x
  rownames(pca_embeddings_atac) <- cell_ids

  # Adjust k_neighbors
  if (atac_k_neighbors >= num_cells) {
    warning("'atac_k_neighbors' (", atac_k_neighbors, ") >= number of cells (", num_cells, "). Setting k_neighbors to num_cells - 1.")
    atac_k_neighbors <- max(1, num_cells - 1)
  }
  if (atac_k_neighbors < 1) stop("Error: 'atac_k_neighbors' must be at least 1 for KNN graph computation after adjustments.")

  message("  Computing k-nearest neighbors (k = ", atac_k_neighbors, ") on ATAC PCs using FNN::get.knn...")
  knn_result <- FNN::get.knn(pca_embeddings_atac, k = atac_k_neighbors)
  if (is.null(knn_result) || is.null(knn_result$nn.index) || is.null(knn_result$nn.dist)) {
    stop("FNN::get.knn returned NULL or empty results (nn.index/nn.dist) for ATAC graph.")
  }

  adj_list_j_raw <- as.vector(knn_result$nn.index)
  distances_raw <- as.vector(knn_result$nn.dist)
  adj_list_i_all <- rep(1:num_cells, each = atac_k_neighbors)

  valid_entries <- !is.na(adj_list_j_raw) & adj_list_j_raw > 0 & adj_list_j_raw <= num_cells &
    is.finite(distances_raw) & distances_raw >= 0
  if (sum(valid_entries) == 0) {
    warning("No valid neighbors found after filtering. Returning empty ATAC graph.")
    empty_mat <- Matrix::Matrix(0, num_cells, num_cells, sparse = TRUE, dimnames = list(cell_ids, cell_ids))
    return(list(graph = empty_mat, pca = atac_pca_res))
  }

  adj_list_i <- adj_list_i_all[valid_entries]
  adj_list_j <- adj_list_j_raw[valid_entries]
  distances <- distances_raw[valid_entries]

  sigma <- median(distances, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) sigma <- max(distances, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1

  similarity <- pmin(1, exp(- (distances^2) / (2 * sigma^2)))
  similarity[is.na(similarity) | !is.finite(similarity)] <- 0
  similarity[similarity < 0] <- 0

  final_idx <- which(similarity > 0)
  if (length(final_idx) == 0) {
    warning("No positive similarities after Gaussian kernel. Returning empty ATAC graph.")
    empty_mat <- Matrix::Matrix(0, num_cells, num_cells, sparse = TRUE, dimnames = list(cell_ids, cell_ids))
    return(list(graph = empty_mat, pca = atac_pca_res))
  }

  adj_list_i <- as.integer(adj_list_i[final_idx])
  adj_list_j <- as.integer(adj_list_j[final_idx])
  similarity <- as.numeric(similarity[final_idx])

  atac_cell_graph_init <- Matrix::sparseMatrix(
    i = adj_list_i,
    j = adj_list_j,
    x = similarity,
    dims = c(num_cells, num_cells),
    dimnames = list(cell_ids, cell_ids),
    giveCsparse = TRUE
  )

  diag(atac_cell_graph_init) <- 1

  message("  Symmetrizing graph (max rule)...")
  atac_triplets <- Matrix::summary(atac_cell_graph_init)
  combined_triplets_df <- data.frame(
    i = c(atac_triplets$i, atac_triplets$j),
    j = c(atac_triplets$j, atac_triplets$i),
    x = c(atac_triplets$x, atac_triplets$x)
  )
  if (!requireNamespace("data.table", quietly = TRUE)) {
    symmetric_edges <- aggregate(x ~ i + j, data = combined_triplets_df, FUN = max)
  } else {
    data.table::setDT(combined_triplets_df)
    symmetric_edges <- combined_triplets_df[, .(x = max(x)), by = .(i, j)]
  }

  atac_cell_graph <- Matrix::sparseMatrix(
    i = symmetric_edges$i,
    j = symmetric_edges$j,
    x = symmetric_edges$x,
    dims = dim(atac_cell_graph_init),
    dimnames = dimnames(atac_cell_graph_init),
    giveCsparse = TRUE
  )

  message("ATAC-based cell graph created. Dimensions: ", paste(dim(atac_cell_graph), collapse = "x"))
  return(list(graph = atac_cell_graph, pca = atac_pca_res))
}


# #' Construct Combined RNA+ImputedATAC Cell Graph
# #'
# #' \code{constructRNATACGraph()} creates a unified cell-to-cell similarity graph
# #' where nodes are RNA spots, and ATAC feature information is imputed onto these
# #' RNA spots based on spatial proximity. It then combines RNA feature similarities
# #' with these imputed ATAC feature similarities.
# #'
# #' @param seurat_obj_rna_preprocessed Preprocessed Seurat object for spatial RNA-seq.
# #' @param rna_spatial_coords_df Data.frame of RNA spatial coordinates (cell_id, x, y).
# #' @param atac_counts_matrix_filtered Sparse matrix of ATAC-seq peak counts (cells x peaks).
# #' @param atac_peaks_gr_filtered `GRanges` object of ATAC-seq peaks.
# #' @param atac_spatial_coords_filtered Data.frame of ATAC spatial coordinates (cell_id, x, y).
# #' @param genome Genome build (e.g., "hg38") for ATAC gene activity.
# #' @param rna_feature_dims Integer vector of PCs for RNA feature graph (e.g., 1:30).
# #' @param rna_feature_k_neighbors Integer, k for RNA feature KNN graph.
# #' @param atac_impute_k_neighbors Integer, k for imputing ATAC to RNA spatially.
# #' @param atac_feature_sim_method Character, method for ATAC feature similarity (e.g., "sparseJaccard").
# #' @param modality_weights Numeric vector of length 2, weights for RNA Feature
# #'                         and Imputed ATAC Feature graphs, respectively. Will be normalized.
# #' @return A sparseMatrix representing the combined RNA+ImputedATAC cell-to-cell graph.
# #' @export
# constructRNATACGraph <- function(
    #     seurat_obj_rna_preprocessed,
#     rna_spatial_coords_df,
#     atac_counts_matrix_filtered,
#     atac_peaks_gr_filtered,
#     atac_spatial_coords_filtered,
#     genome = "hg38",
#     rna_feature_dims = 1:30,
#     rna_feature_k_neighbors = 20,
#     atac_impute_k_neighbors = 5,
#     atac_feature_sim_method = "sparseJaccard",
#     modality_weights = c(0.5, 0.5)
# ) {
#   message("Starting `constructRNATACGraph` pipeline...")
#   rna_cell_ids <- colnames(seurat_obj_rna_preprocessed)
#   if (length(rna_cell_ids) == 0) {
#     stop("No RNA cells remaining for graph construction. Please check preprocessing steps.")
#   }
#
#   # 1. Create RNA Feature Graph
#   rna_feature_graph <- createRNACellGraph(seurat_obj_rna_preprocessed, dims = rna_feature_dims, k_neighbors = rna_feature_k_neighbors)
#   if (is.null(rna_feature_graph) || nrow(rna_feature_graph) == 0) {
#     warning("RNA Feature Graph could not be created meaningfully (empty or NULL). It will be represented as zeros in the combined graph.")
#     rna_feature_graph <- Matrix::Matrix(0, length(rna_cell_ids), length(rna_cell_ids), sparse = TRUE, dimnames = list(rna_cell_ids, rna_cell_ids))
#   }
#
#   # 2. Create ATAC Gene Activity Matrix (using custom peak-based aggregation)
#   atac_gene_activity_matrix <- createPeakBasedGeneActivityMatrix(
#     atac_counts_matrix = atac_counts_matrix_filtered,
#     atac_peaks_gr = atac_peaks_gr_filtered,
#     genome = genome
#   )
#
#   # 3. Impute ATAC Gene Activity onto RNA spots
#   imputed_atac_gene_activity_rna_cells <- imputeATACtoRNAbySpatial(
#     rna_spatial_coords = rna_spatial_coords_df,
#     atac_spatial_coords = atac_spatial_coords_filtered,
#     atac_gene_activity_matrix = atac_gene_activity_matrix,
#     k_impute_spatial = atac_impute_k_neighbors
#   )
#
#   # 4. Create Imputed ATAC Feature Graph (on RNA spots)
#   if (nrow(imputed_atac_gene_activity_rna_cells) > 0 && ncol(imputed_atac_gene_activity_rna_cells) > 0) {
#     imputed_atac_feature_graph <- computeCellSim(imputed_atac_gene_activity_rna_cells, method = atac_feature_sim_method)
#     # Ensure graph is aligned to RNA cell IDs
#     rownames(imputed_atac_feature_graph) <- colnames(imputed_atac_feature_graph) <- rna_cell_ids
#   } else {
#     warning("Imputed ATAC gene activity matrix is empty, cannot create imputed ATAC feature graph. It will be represented as zeros.")
#     imputed_atac_feature_graph <- Matrix::Matrix(0, length(rna_cell_ids), length(rna_cell_ids), sparse = TRUE, dimnames = list(rna_cell_ids, rna_cell_ids))
#   }
#
#   # Ensure all graphs have the correct cell IDs and order for combination
#   rna_feature_graph_ordered <- rna_feature_graph[rna_cell_ids, rna_cell_ids]
#   imputed_atac_feature_graph_ordered <- imputed_atac_feature_graph[rna_cell_ids, rna_cell_ids]
#
#   # 5. Combine RNA Feature and Imputed ATAC Feature Graphs
#   list_of_graphs_to_combine <- list(
#     rna_feature_graph_ordered,
#     imputed_atac_feature_graph_ordered
#   )
#
#   combined_graph <- combineMultiModalCellGraphs(list_of_graphs_to_combine, modality_weights)
#   message("`constructRNATACGraph` pipeline complete. Final combined graph dimensions: ", paste(dim(combined_graph), collapse = "x"))
#   return(combined_graph)
# }


#' Generate Separate RNA and ATAC Cell-to-Cell Similarity Graphs
#'
#' \code{generateSeparateModalGraphs()} creates independent cell-to-cell similarity graphs
#' for spatial RNA-seq and spatial ATAC-seq data. The RNA graph is based on PCA embeddings,
#' and the ATAC graph is based on gene activity scores derived from peak aggregation,
#' followed by sparse PCA and KNN graph construction.
#'
#' @param seurat_obj_rna_preprocessed Preprocessed Seurat object for spatial RNA-seq.
#' @param atac_counts_matrix_filtered Sparse matrix of ATAC-seq peak counts (cells x peaks).
#' @param atac_peaks_gr_filtered `GRanges` object of ATAC-seq peaks.
#' @param genome Genome build (e.g., "hg38") for ATAC gene activity.
#' @param rna_feature_dims Integer vector of PCs for RNA feature graph (e.g., 1:30).
#' @param rna_feature_k_neighbors Integer, k for RNA feature KNN graph.
#' @param atac_pca_dims Integer, number of principal components to use for ATAC graph construction.
#' @param atac_k_neighbors Integer, the number of nearest neighbors to consider for each cell
#'                         in the ATAC KNN graph.
#' @return A list containing:
#'   \itemize{
#'     \item `rna_graph`: Sparse matrix of RNA-based cell-to-cell similarity.
#'     \item `atac_graph`: Sparse matrix of ATAC-based cell-to-cell similarity.
#'     \item `atac_gene_activity_matrix`: The gene activity matrix derived from ATAC peaks (useful for ATAC-specific cell typing or visualization).
#'   }
#' @export
generateSeparateModalGraphs <- function(
    seurat_obj_rna_preprocessed,
    atac_counts_matrix_filtered,
    atac_peaks_gr_filtered,
    genome = "hg38",
    rna_feature_dims = 1:30,
    rna_feature_k_neighbors = 20,
    atac_pca_dims = 30,
    atac_k_neighbors = 20
) {
  message("Starting `generateSeparateModalGraphs` pipeline...")

  # 1. Create RNA Feature Graph
  rna_graph <- createRNACellGraph(seurat_obj_rna_preprocessed, dims = rna_feature_dims, k_neighbors = rna_feature_k_neighbors)
  rna_cell_ids <- colnames(seurat_obj_rna_preprocessed)
  if (is.null(rna_graph) || nrow(rna_graph) == 0) {
    warning("RNA Feature Graph could not be created meaningfully (empty or NULL). Returning a zero matrix for RNA graph.")
    rna_graph <- Matrix::Matrix(0, length(rna_cell_ids), length(rna_cell_ids), sparse = TRUE, dimnames = list(rna_cell_ids, rna_cell_ids))
  }

  # 2. Create ATAC Gene Activity Matrix
  atac_gene_activity_matrix <- createPeakBasedGeneActivityMatrix(
    atac_counts_matrix = atac_counts_matrix_filtered,
    atac_peaks_gr = atac_peaks_gr_filtered,
    genome = genome
  )

  # 3. Create ATAC Feature Graph (using new PCA+KNN approach)
  atac_graph_res <- createATACCellGraph(
    atac_gene_activity_matrix,
    atac_pca_dims = atac_pca_dims,
    atac_k_neighbors = atac_k_neighbors
  )
  if (is.list(atac_graph_res) && !is.null(atac_graph_res$graph)) {
    atac_graph <- atac_graph_res$graph
    atac_pca_res <- atac_graph_res$pca
  } else if (is(atac_graph_res, "Matrix")) {
    atac_graph <- atac_graph_res
    atac_pca_res <- NULL
  } else {
    atac_graph <- Matrix::Matrix(0, nrow(atac_gene_activity_matrix), nrow(atac_gene_activity_matrix),
                                 sparse = TRUE, dimnames = list(rownames(atac_gene_activity_matrix), rownames(atac_gene_activity_matrix)))
    atac_pca_res <- NULL
  }
  atac_cell_ids <- rownames(atac_gene_activity_matrix) # Use IDs from the gene activity matrix
  if (is.null(atac_cell_ids)) atac_cell_ids <- seq_len(nrow(atac_gene_activity_matrix))
  if (is.null(atac_graph) || nrow(atac_graph) == 0) {
    warning("ATAC Feature Graph could not be created meaningfully (empty or NULL). Returning a zero matrix for ATAC graph.")
    atac_graph <- Matrix::Matrix(0, length(atac_cell_ids), length(atac_cell_ids), sparse = TRUE, dimnames = list(atac_cell_ids, atac_cell_ids))
  }

  message("`generateSeparateModalGraphs` pipeline complete. Returning separate RNA and ATAC graphs.")
  return(list(
    rna_graph = rna_graph,
    atac_graph = atac_graph,
    atac_gene_activity_matrix = atac_gene_activity_matrix,
    atac_pca_res = if (exists("atac_pca_res")) atac_pca_res else NULL
  ))
}


#' Compute Label Edges
#'
#' \code{computeTypeEdges} generates a matrix of edges from each cell type label to each cell from gene expression.
#' The edge weight is normalized gene expression of markers weighted by the log2FC.
#'
#' @param exprMat_norm a gene by cell data.frame or matrix, rownames are genes, colnames are cell barcodes
#' @param markers marker genes for each cell type, columns are gene, cluster, p_val_adj (optional) and avg_log2FC (optional).
#' @param pval.cutoff select markers with adjust pvalue (\code{p_val_adj}) < pval.cutoff. Default: 0.05
#' @param log2FC.cutoff select markers with abs(log2FC.cutoff)> 0.5 if only.pos = F or log2FC.cutoff>0.5 if only.pos = T. Default: 0.5
#' @param only.pos only include positive markers
#' @param force if true (default), this function will use marker genes that matched the name in exprMat_norm to compute edge weight;
#'              if false, this function will raise an error if some marker genes are not in exprMat_norm.
#' @return a matrix of edges from each cell type label to each cell
#' @export
#' @import data.table
changedCellWalkercomputeTypeEdges <- function(exprMat_norm, markers, pval.cutoff = 0.05, log2FC.cutoff = 0.5, only.pos = FALSE, force = TRUE)
{
  if(!requireNamespace("data.table", quietly = TRUE)){
    stop("Must install data.table")
  }
  if(missing(exprMat_norm) || (!is(exprMat_norm, "data.frame") & !is(exprMat_norm, "matrix")  & !is(exprMat_norm, "Matrix"))){
    stop("Must provide a dataframe or matrix of RNA data")
  }
  if(is.null(colnames(exprMat_norm)) | is.null(rownames(exprMat_norm)))
  {
    stop('exprMat_norm must have column and row names')
  }
  if(missing(markers) || !is(markers, "data.frame")){
    stop("Must provide a dataframe of markers")
  }
  if(is.null(markers$gene) || is.null(markers$cluster)){
    stop("markers must have 'gene' and 'cluster' columns")
  }


  if(is.null(markers$avg_log2FC)) {
    warning('avg_log2FC is not present in the columns of markers, so markers will not be filtered by avg_log2FC')
    markers$avg_log2FC = log2FC.cutoff # Assign a default, but note it's not a true value
  }
  if(is.null(markers$p_val_adj)) {
    warning('p_val_adj is not present in the columns of markers, so markers will not be filtered by p_val_adj')
    markers$p_val_adj =  pval.cutoff # Assign a default, but note it's not a true value
  }

  markers = data.table(markers)
  if(only.pos)
  {
    markers = markers[avg_log2FC >= log2FC.cutoff & p_val_adj <= pval.cutoff]
  }else{
    markers = markers[abs(avg_log2FC) >= log2FC.cutoff & p_val_adj <= pval.cutoff]
  }

  if(force == FALSE)
  {
    stopifnot(all(markers$gene %in% rownames(exprMat_norm)))
  }

  # Filter markers to only include those present in the expression matrix
  markers_inter <- markers[markers$gene %in% rownames(exprMat_norm)]

  if(length(unique(markers_inter$cluster)) < length(unique(markers$cluster)))
  {
    # This check needs to be more careful if any cluster completely loses all its markers
    # It might be better to check each unique cluster from the original markers.
    original_clusters <- unique(markers$cluster)
    clusters_with_markers_after_filter <- unique(markers_inter$cluster)
    missing_clusters <- setdiff(original_clusters, clusters_with_markers_after_filter)
    if (length(missing_clusters) > 0) {
      stop(paste("The following clusters lost all their markers after filtering and intersection with expression matrix:",
                 paste(missing_clusters, collapse = ", "),
                 ". Need to lower the pval/logFC thresholds or check if expression matrix contains the markers."))
    }
  }

  # Ensure genes in exprMat_norm are ordered consistently with markers_inter for the calculation
  # This might not be strictly necessary if .SD[['gene']] is always correctly indexed,
  # but it's a good practice to be aware of potential indexing mismatches with sparse matrices.

  # Use drop = FALSE to ensure subsetting always returns a matrix-like object
  labelEdges = markers_inter[, {
    # Extract genes for the current cluster
    current_genes <- .SD[['gene']]

    # Subset exprMat_norm, ensuring it remains a matrix even for single genes
    subset_expr_mat <- exprMat_norm[current_genes, , drop = FALSE]

    # Extract log2FC values, ensuring they align with the subset_expr_mat rows
    # This step is crucial if current_genes is not guaranteed to be in the same order as rownames(subset_expr_mat)
    # However, since we filter markers_inter by exprMat_norm rownames, and .SD[['gene']] comes from markers_inter,
    # the order should be preserved for the multiplication.
    # We explicitly match just to be extra safe, though for direct multiplication,
    # the vector is recycled against the rows of the matrix, which usually works if lengths match.
    # To be absolutely sure, it might be safer to ensure avg_log2FC is aligned.

    # A simpler approach, relying on vector recycling, is usually sufficient and less code:
    gene_log2FC <- .SD[['avg_log2FC']]

    # Perform element-wise multiplication and then colSums
    # The matrix subset_expr_mat will have genes as rows, cells as columns
    # gene_log2FC will be a vector of length = number of genes in current_genes
    # R will recycle gene_log2FC down the columns of subset_expr_mat for multiplication
    weighted_expr <- subset_expr_mat * gene_log2FC

    # Check for zero sum of absolute log2FC to avoid division by zero
    sum_abs_log2FC <- sum(abs(gene_log2FC))
    score_val <- if (sum_abs_log2FC > 0) {
      colSums(weighted_expr) / sum_abs_log2FC
    } else {
      rep(0, ncol(exprMat_norm)) # If no valid log2FC, score is 0
    }

    list("score" = score_val, "cell" = colnames(exprMat_norm))
  }, by = cluster]

  labelEdges = reshape2::acast(labelEdges, cell~cluster, value.var = 'score')
  labelEdges[labelEdges <0] = 0

  labelEdges
}


# --- Plotting functions for quick visualization and validation ---

# (User provided this, including here for completeness of testing)
plot_spatial_rna_atac <- function(rna_spatial_coords_df, atac_spatial_coords_df,
                                  color_by_rna = NULL, color_by_atac = NULL,
                                  pt_size = 0.6, pt_alpha = 0.8,
                                  same_limits = FALSE,
                                  rna_title = "RNA (spatial)",
                                  atac_title = "ATAC (spatial)") {

  if (nrow(rna_spatial_coords_df) == 0 && nrow(atac_spatial_coords_df) == 0) {
    message("No spatial coordinates to plot.")
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data to plot") + theme_void())
  }

  # Common limits
  if (same_limits) {
    xlim <- range(c(rna_spatial_coords_df$x, atac_spatial_coords_df$x), na.rm = TRUE)
    ylim <- range(c(rna_spatial_coords_df$y, atac_spatial_coords_df$y), na.rm = TRUE)
  } else { xlim <- ylim <- NULL }

  make_plot <- function(df, title, color_col) {
    p <- ggplot(df, aes(x = x, y = y)) +
      coord_fixed() + ggtitle(title) + theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    if (!is.null(xlim)) p <- p + xlim(xlim) + ylim(ylim)
    if (!is.null(color_col)) {
      if (is.numeric(df[[color_col]])) {
        p <- p + geom_point(aes_string(color = color_col), size = pt_size, alpha = pt_alpha) +
          scale_color_viridis_c()
      } else {
        p <- p + geom_point(aes_string(color = color_col), size = pt_size, alpha = pt_alpha) +
          scale_color_brewer(palette = "Set1", na.value = "grey70")
      }
    } else {
      p <- p + geom_point(color = "#2c3e50", size = pt_size, alpha = pt_alpha)
    }
    p + xlab("x") + ylab("y")
  }

  # Only create plots if data exists
  p_rna <- if (nrow(rna_spatial_coords_df) > 0) make_plot(rna_spatial_coords_df, rna_title, color_by_rna) else NULL
  p_atac <- if (nrow(atac_spatial_coords_df) > 0) make_plot(atac_spatial_coords_df, atac_title, color_by_atac) else NULL

  if (is.null(p_rna) && is.null(p_atac)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data to plot") + theme_void())
  } else if (is.null(p_rna)) {
    return(p_atac)
  } else if (is.null(p_atac)) {
    return(p_rna)
  } else {
    plot_grid(p_rna, p_atac, ncol = 1, align = "v", rel_heights = c(1, 1))
  }
}


#' Plot Spatial and Graph Coherence for Validation
#'
#' \code{plotSpatialGraphCoherence()} visualizes the relationship between graph-based
#' neighbors and spatially proximate neighbors for a subset of cells.
#' This helps to validate if the non-spatial (RNA+ATAC) graph preserves spatial patterns.
#'
#' @param combined_rna_atac_graph The combined cell-to-cell graph (sparse matrix, RNA spots as nodes).
#' @param rna_spatial_coords_df Data.frame of RNA spatial coordinates (cell_id, x, y).
#' @param n_query_cells Integer, number of random query cells to pick for plotting.
#' @param k_neighbors_to_show Integer, number of graph/spatial neighbors to highlight for each query cell.
#' @param background_alpha Numeric, alpha value for 'Other' cells (0 to 1). Lower means more transparent.
#' @param neighbor_pt_size Numeric, point size for neighbor cells.
#' @param query_pt_size Numeric, point size for query cells.
#' @return A ggplot object visualizing query cells, their graph neighbors, and their spatial neighbors.
#' @export
plotSpatialGraphCoherence <- function(combined_rna_atac_graph,
                                      rna_spatial_coords_df,
                                      n_query_cells = 4,
                                      k_top_neighbors = 10,
                                      k_next_neighbors = 10,
                                      background_alpha = 0.05,
                                      neighbor_pt_size = 1.8,
                                      query_pt_size = 3.5) {
  # Distinguish Top vs Next neighbor groups, and to visualize overlaps between Graph and Spatial neighbor ranks.

  if (nrow(rna_spatial_coords_df) == 0 || nrow(combined_rna_atac_graph) == 0) {
    message("No cells available for plotting spatial graph coherence.")
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data to plot") + theme_void())
  }

  num_cells <- nrow(rna_spatial_coords_df)
  cell_ids <- rownames(rna_spatial_coords_df)

  k_total_neighbors_to_retrieve <- k_top_neighbors + k_next_neighbors
  if (k_total_neighbors_to_retrieve >= num_cells) {
    warning("Total requested neighbors (k_top + k_next = ", k_total_neighbors_to_retrieve, ") is >= number of cells (", num_cells, "). Adjusting to retrieve num_cells - 1 total neighbors and splitting them as best as possible.")
    k_total_neighbors_to_retrieve <- max(1, num_cells - 1)
    if (k_total_neighbors_to_retrieve <= k_top_neighbors) {
      k_next_neighbors <- 0
    } else {
      k_next_neighbors <- k_total_neighbors_to_retrieve - k_top_neighbors
    }
    k_top_neighbors <- min(k_top_neighbors, k_total_neighbors_to_retrieve)
  }
  if (k_top_neighbors == 0 && k_next_neighbors == 0) {
    message("No neighbors can be shown (k_top_neighbors and k_next_neighbors are 0). Only query cells will be plotted with 'Other' cells in background.")
  }

  query_cells <- sample(cell_ids, min(n_query_cells, num_cells))
  message("Analyzing graph and spatial coherence for ", length(query_cells), " query cells (top=", k_top_neighbors, ", next=", k_next_neighbors, ")...")
  cells_to_plot_per_query_list <- list()

  for (q_id in query_cells) {
    graph_sims <- combined_rna_atac_graph[q_id, ]
    graph_neighbors_indices_all_sorted <- order(graph_sims, decreasing = TRUE)
    q_id_idx_in_graph <- match(q_id, cell_ids)
    graph_neighbors_indices_all_sorted <- graph_neighbors_indices_all_sorted[graph_neighbors_indices_all_sorted != q_id_idx_in_graph]
    actual_neighbors_available <- length(graph_neighbors_indices_all_sorted)
    k_total_actual <- min(k_total_neighbors_to_retrieve, actual_neighbors_available)
    selected_graph_neighbor_indices <- graph_neighbors_indices_all_sorted[1:k_total_actual]
    selected_graph_neighbor_ids <- cell_ids[selected_graph_neighbor_indices]

    spatial_coords_matrix <- as.matrix(rna_spatial_coords_df[, c("x", "y")])
    spatial_nn_result <- RANN::nn2(
      query = spatial_coords_matrix[q_id, , drop = FALSE],
      data = spatial_coords_matrix,
      k = min(k_total_actual + 1, num_cells)
    )
    spatial_neighbors_indices_raw <- spatial_nn_result$nn.idx[1, ]
    self_idx <- match(q_id, rownames(spatial_coords_matrix))
    spatial_neighbors_indices_raw <- spatial_neighbors_indices_raw[spatial_neighbors_indices_raw != self_idx & spatial_neighbors_indices_raw > 0]
    actual_spatial_available <- length(spatial_neighbors_indices_raw)
    k_total_spatial_actual <- min(k_total_actual, actual_spatial_available)
    selected_spatial_neighbor_indices <- spatial_neighbors_indices_raw[1:k_total_spatial_actual]
    selected_spatial_neighbor_ids <- cell_ids[selected_spatial_neighbor_indices]

    # Assign neighbor types (Top vs Next) for both graph and spatial selections
    graph_neighbor_types <- rep("Other", length(selected_graph_neighbor_ids))
    spatial_neighbor_types <- rep("Other", length(selected_spatial_neighbor_ids))

    if (k_top_neighbors > 0 && length(selected_graph_neighbor_indices) >= k_top_neighbors) {
      top_graph_indices_in_selection <- 1:k_top_neighbors
      graph_neighbor_types[top_graph_indices_in_selection] <- "Graph (Top)"
    }
    if (k_next_neighbors > 0 && length(selected_graph_neighbor_indices) > k_top_neighbors) {
      next_graph_indices_in_selection <- (k_top_neighbors + 1):min(k_top_neighbors + k_next_neighbors, actual_neighbors_available)
      graph_neighbor_types[next_graph_indices_in_selection] <- "Graph (Next)"
    }

    if (k_top_neighbors > 0 && length(selected_spatial_neighbor_indices) >= k_top_neighbors) {
      top_spatial_indices_in_selection <- 1:k_top_neighbors
      spatial_neighbor_types[top_spatial_indices_in_selection] <- "Spatial (Top)"
    }
    if (k_next_neighbors > 0 && length(selected_spatial_neighbor_indices) > k_top_neighbors) {
      next_spatial_indices_in_selection <- (k_top_neighbors + 1):min(k_top_neighbors + k_next_neighbors, actual_spatial_available)
      spatial_neighbor_types[next_spatial_indices_in_selection] <- "Spatial (Next)"
    }

    current_query_cells_df_base <- data.frame(
      cell      = cell_ids,
      x         = rna_spatial_coords_df[cell_ids, "x"],
      y         = rna_spatial_coords_df[cell_ids, "y"],
      query_id  = q_id,
      type_final = "Other",
      stringsAsFactors = FALSE
    )
    current_query_cells_df_base$type_final[current_query_cells_df_base$cell == q_id] <- "Query"

    graph_neighbor_data <- data.frame(cell = selected_graph_neighbor_ids, neighbor_type = graph_neighbor_types, stringsAsFactors = FALSE)
    spatial_neighbor_data <- data.frame(cell = selected_spatial_neighbor_ids, neighbor_type = spatial_neighbor_types, stringsAsFactors = FALSE)
    combined_neighbors <- unique(rbind(graph_neighbor_data, spatial_neighbor_data))

    for(i in 1:nrow(combined_neighbors)) {
      cell <- combined_neighbors$cell[i]
      is_graph_top <- cell %in% selected_graph_neighbor_ids[graph_neighbor_types == "Graph (Top)"]
      is_graph_next <- cell %in% selected_graph_neighbor_ids[graph_neighbor_types == "Graph (Next)"]
      is_spatial_top <- cell %in% selected_spatial_neighbor_ids[spatial_neighbor_types == "Spatial (Top)"]
      is_spatial_next <- cell %in% selected_spatial_neighbor_ids[spatial_neighbor_types == "Spatial (Next)"]
      final_type <- "Other"
      if (cell == q_id) final_type <- "Query"
      else if (is_graph_top && is_spatial_top) final_type <- "Both (Top)"
      else if (is_graph_next && is_spatial_next) final_type <- "Both (Next)"
      else if (is_graph_top && is_spatial_next) final_type <- "Overlap (Graph Top/Spatial Next)"
      else if (is_graph_next && is_spatial_top) final_type <- "Overlap (Graph Next/Spatial Top)"
      else if (is_graph_top) final_type <- "Graph (Top)"
      else if (is_graph_next) final_type <- "Graph (Next)"
      else if (is_spatial_top) final_type <- "Spatial (Top)"
      else if (is_spatial_next) final_type <- "Spatial (Next)"
      current_query_cells_df_base$type_final[current_query_cells_df_base$cell == cell] <- final_type
    }
    cells_to_plot_per_query_list[[q_id]] <- current_query_cells_df_base
  }

  if (length(cells_to_plot_per_query_list) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No query cells or neighbors to plot") + theme_void())
  }
  final_plot_data <- do.call(rbind, cells_to_plot_per_query_list)

  defined_plot_type_levels <- c(
    "Query",
    "Both (Top)",
    "Graph (Top)",
    "Spatial (Top)",
    "Both (Next)",
    "Graph (Next)",
    "Spatial (Next)",
    "Overlap (Graph Top/Spatial Next)",
    "Overlap (Graph Next/Spatial Top)",
    "Other"
  )
  final_plot_data$type_final <- factor(final_plot_data$type_final, levels = defined_plot_type_levels)

  color_values <- c(
    "Query" = "red",
    "Both (Top)" = "#d73027",
    "Graph (Top)" = "#fdae61",
    "Spatial (Top)" = "#a6d96a",
    "Both (Next)" = "#7570b3",
    "Graph (Next)" = "#1b9e77",
    "Spatial (Next)" = "#e7298a",
    "Overlap (Graph Top/Spatial Next)" = "#b2df8a",
    "Overlap (Graph Next/Spatial Top)" = "#CAB2D6",
    "Other" = "grey"
  )
  present_levels <- intersect(names(color_values), levels(final_plot_data$type_final))
  color_values_present <- color_values[present_levels]

  p <- ggplot2::ggplot(final_plot_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(data = dplyr::filter(final_plot_data, type_final == "Other"),
                        ggplot2::aes(color = type_final), size = 0.5, alpha = background_alpha) +
    ggplot2::geom_point(data = dplyr::filter(final_plot_data, type_final %in% c("Spatial (Top)", "Graph (Top)")),
                        ggplot2::aes(color = type_final), size = neighbor_pt_size, alpha = 0.8, shape = 11) +
    ggplot2::geom_point(data = dplyr::filter(final_plot_data, type_final %in% c("Spatial (Next)", "Graph (Next)")),
                        ggplot2::aes(color = type_final), size = neighbor_pt_size - 0.3, alpha = 0.7, shape = 17) +
    ggplot2::geom_point(data = dplyr::filter(final_plot_data, type_final %in% c("Both (Top)", "Both (Next)", "Overlap (Graph Top/Spatial Next)", "Overlap (Graph Next/Spatial Top)")),
                        ggplot2::aes(color = type_final), size = neighbor_pt_size + 0.2, alpha = 0.9, shape = 1) +
    ggplot2::geom_point(data = dplyr::filter(final_plot_data, type_final == "Query"),
                        ggplot2::aes(color = type_final), size = query_pt_size, shape = 8, stroke = 1.2) +
    ggplot2::scale_color_manual(values = color_values_present, na.value = "grey") +
    ggplot2::facet_wrap(~query_id) +
    ggplot2::labs(title = paste0("Graph & Spatial Neighborhood Comparison (k_top=", k_top_neighbors, ", k_next=", k_next_neighbors, ")"), color = "Cell Category") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "bottom",
                   panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   plot.background = ggplot2::element_rect(fill = "white", color = NA),
                   axis.title = ggplot2::element_text(size = 8),
                   strip.text = ggplot2::element_text(size = 8)) +
    ggplot2::coord_fixed()

  message("Spatial graph coherence plot generated comparing neighbor groups.")
  return(p)
}


#' Assess Overall Spatial-Graph Coherence (Revised)
#'
#' \code{assessSpatialGraphCoherenceOverall()} quantifies and visualizes the global coherence
#' between feature-based cell similarity graphs and spatial proximity for a given modality.
#' This version emphasizes speed and robustness for large sparse graphs (~50k cells) and:
#' - Computes Median Spatial Distance to Feature-Graph Neighbors (MSD-FGN) using the top-k
#'   neighbors from the provided feature_graph (symmetrized).
#' - Computes Median Feature Similarity to Spatial Neighbors (MFS-SN) using cosine similarity
#'   on the provided feature matrix (fast per-cell dot-products with precomputed norms) with PCA embedding.
#' - Computes Median Feature Similarity to Feature-Graph Neighbors (MFS-FGN) as a sanity check.
#' - The "Feature-Graph Distance to Spatial Neighbors" metric (MFD-SN) has been removed (better calculation takes too long, simple calculation gives no valid results).
#' - The Combined Coherence Score (CCS) averages normalized MSD-FGN and MFS-SN.
#' - All per-cell heavy work runs in parallel batches via future.apply.
#'
#' @param feature_graph A sparse matrix (cells x cells) representing feature-based
#'   cell-to-cell similarities (e.g., RNA graph or ATAC graph). Rownames/Colnames must be cell IDs.
#' @param feature_matrix_for_similarity A sparse matrix (features x cells) used for calculating
#'   cosine similarity between cells for the "spatial neighbors" metric. Rownames are features,
#'   colnames are cell IDs, and must match `feature_graph`.
#' @param spatial_coords_df A data.frame (cells x 2) of spatial x, y coordinates.
#'   Rownames must be cell IDs and match those in `feature_graph`.
#' @param k_neighbors Integer, the number of nearest neighbors to consider for coherence assessment.
#' @param modality_name Character string (e.g., "RNA", "ATAC") for labeling results and plots.
#' @param n_workers Integer, number of parallel workers. NULL defaults to detectCores() - 1.
#' @param batch_size Integer, batch size for parallel processing.
#' @param n_pca_embed Integer, number of PCA dimensions for embedding features. Set NA or 0 to skip.
#' @return A list containing:
#'   \itemize{
#'     \item `coherence_df`: A data.frame with metrics for each cell.
#'     \item `plot_dist_spatial_for_graph_nn`: ggplot object for MSD-FGN density.
#'     \item `plot_sim_graph_for_spatial_nn`: ggplot object for MFS-SN density.
#'     \item `plot_sim_graph_for_graph_nn`: ggplot object for MFS-FGN density.
#'     \item `plot_spatial_coherence_metrics`: ggplot object with spatial subplots for selected metrics.
#'     \item `diagnostics`: List containing computational diagnostics.
#'   }
#' @export
assessSpatialGraphCoherenceOverall <- function(
    feature_graph,
    feature_matrix_for_similarity,
    spatial_coords_df,
    k_neighbors = 10,
    modality_name = "Modality",
    n_workers = NULL,        # number of parallel workers; NULL -> auto (detectCores - 1)
    batch_size = 1000L,      # batch size for parallel processing
    n_pca_embed = 50L,      # set NA to skip PCA embedding
    neighbor_rank_range = c(1, 10)
) {

  # --- Dependencies check ---
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install 'Matrix'.")
  if (!requireNamespace("RANN", quietly = TRUE)) stop("Install 'RANN'.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install 'data.table'.")
  if (!requireNamespace("future.apply", quietly = TRUE)) stop("Install 'future.apply'.")
  # irlba optional for PCA embedding
  use_irlba <- requireNamespace("irlba", quietly = TRUE)

  # --- Input validation ---
  if (!is(feature_graph, "Matrix")) stop("`feature_graph` must be a sparse Matrix (e.g., dgCMatrix).")
  if (!is(feature_matrix_for_similarity, "Matrix")) stop("`feature_matrix_for_similarity` must be a sparse Matrix (e.g., dgCMatrix).")
  if (!is.data.frame(spatial_coords_df) || ncol(spatial_coords_df) < 2) stop("`spatial_coords_df` must be a data.frame with at least two columns (x,y).")
  if (is.null(rownames(feature_graph)) || is.null(rownames(spatial_coords_df)) || is.null(colnames(feature_matrix_for_similarity))) {
    stop("All inputs must have matching cell IDs (rownames for graph/spatial, colnames for feature matrix).")
  }
  if (!is.numeric(neighbor_rank_range) || length(neighbor_rank_range) != 2) {
    stop("`neighbor_rank_range` must be an integer-like numeric vector of length 2, e.g. c(1,10) or c(11,20).")
  }
  neighbor_rank_range <- as.integer(neighbor_rank_range)
  if (any(neighbor_rank_range < 1) || neighbor_rank_range[1] > neighbor_rank_range[2]) {
    stop("`neighbor_rank_range` must be positive with neighbor_rank_range[1] <= neighbor_rank_range[2].")
  }

  # Derive effective k_neighbors based on neighbor_rank_range upper bound for any internal nn computations
  required_k_for_internal_nn <- max(k_neighbors, neighbor_rank_range[2])

  # Align cell sets
  graph_cell_ids <- rownames(feature_graph)
  spatial_cell_ids <- rownames(spatial_coords_df)
  feature_cell_ids <- colnames(feature_matrix_for_similarity)
  common_cells <- Reduce(intersect, list(graph_cell_ids, spatial_cell_ids, feature_cell_ids))
  if (length(common_cells) == 0) stop("No common cell IDs found across inputs.")
  if (length(common_cells) < nrow(feature_graph) || length(common_cells) < nrow(spatial_coords_df) || length(common_cells) < ncol(feature_matrix_for_similarity)) {
    message("Warning: Subsetting to ", length(common_cells), " common cells.")
  }
  feature_graph <- feature_graph[common_cells, common_cells, drop = FALSE]
  spatial_coords_df <- spatial_coords_df[common_cells, , drop = FALSE]
  feature_matrix_for_similarity <- feature_matrix_for_similarity[, common_cells, drop = FALSE]
  num_cells <- nrow(feature_graph)
  cell_ids <- rownames(feature_graph)

  if (num_cells <= 1) {
    message("Too few cells to compute coherence. Returning empty result.")
    empty_plots <- function() ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data") + ggplot2::theme_void()
    return(list(
      coherence_df = data.frame(cell_id = character(0), modality = character(0)),
      plot_dist_spatial_for_graph_nn = empty_plots(),
      plot_sim_graph_for_spatial_nn = empty_plots(),
      plot_sim_graph_for_graph_nn = empty_plots(),
      # plot_dist_graph_for_spatial_nn is removed
      plot_spatial_coherence_metrics = empty_plots(),
      diagnostics = list(note = "Too few cells")
    ))
  }

  if (k_neighbors >= num_cells) {
    warning("k_neighbors >= num_cells; adjusting to num_cells - 1")
    k_neighbors <- max(1, num_cells - 1)
  }
  if (k_neighbors == 0) {
    stop("k_neighbors must be at least 1 after adjustments for coherence assessment.")
  }

  # --- Parallel plan ---
  if (is.null(n_workers)) {
    n_workers <- max(1L, parallel::detectCores() - 1L)
  } else {
    n_workers <- as.integer(n_workers)
    if (n_workers < 1L) n_workers <- 1L
  }
  old_plan <- future::plan()
  options(future.globals.maxSize = +Inf)
  future::plan(future::multisession, workers = n_workers)
  on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)

  message("Assessing spatial-graph coherence for ", modality_name, " (n_cells=", num_cells, ", k=", k_neighbors, ", workers=", n_workers, ")")

  # --- Prepare spatial matrix and precompute norms for cosine ---
  spatial_coords_mat_raw <- as.matrix(spatial_coords_df[, c("x", "y")])

  # Normalize spatial coordinates to a common scale (0 to 1, based on the largest dimension range)
  # This makes distances comparable across different tissues/slides/modalities.
  x_range <- range(spatial_coords_mat_raw[, 1])
  y_range <- range(spatial_coords_mat_raw[, 2])

  max_dim_range <- max(diff(x_range), diff(y_range))

  if (max_dim_range == 0) {
    # All points are identical, cannot normalize distances meaningfully.
    # Use raw coordinates, distances will be 0.
    warning("Spatial coordinates have zero range for modality '", modality_name, "'. Using raw coordinates for distance calculations. Metric 1 will likely be 0.")
    spatial_coords_mat_normalized <- spatial_coords_mat_raw
  } else {
    # Translate coordinates so the minimum is 0, then scale so the maximum range is 1.
    spatial_coords_mat_normalized <- sweep(spatial_coords_mat_raw, 2, c(min(x_range), min(y_range)), "-")
    spatial_coords_mat_normalized <- spatial_coords_mat_normalized / max_dim_range
  }

  # Precompute norms for feature matrix (columns are cells)
  feature_matrix_norms <- sqrt(Matrix::colSums(feature_matrix_for_similarity * feature_matrix_for_similarity))

  # --- Optional PCA embedding to speed up cosine similarity ---
  use_pca_embedding <- !is.null(n_pca_embed) && !is.na(n_pca_embed) && n_pca_embed > 0 && use_irlba
  feature_matrix_for_similarity_to_use <- feature_matrix_for_similarity
  if (use_pca_embedding) {
    message("Computing PCA embedding (n_pca_embed=", n_pca_embed, ") via irlba for faster cosine computations...")
    # Convert to dense matrix for irlba; if too large, user can set n_pca_embed=NA to skip
    dense_mat <- as.matrix(feature_matrix_for_similarity) # This is typically Features x Cells (F x C)
    sv <- irlba::irlba(dense_mat, nv = n_pca_embed)

    # sv$v is C x k (k=n_pca_embed)
    # Transposing gives k x C (PCs x Cells)
    feature_matrix_embed <- t(sv$v %*% diag(sv$d))

    # Assign cell IDs as column names so subsetting works.
    if (ncol(feature_matrix_embed) == length(cell_ids)) {
      colnames(feature_matrix_embed) <- cell_ids # Assign cell IDs as column names
      rownames(feature_matrix_embed) <- paste0("PC", seq_len(n_pca_embed)) # Assign PC names for rows
      feature_matrix_for_similarity_to_use <- feature_matrix_embed
      # Recalculate norms for the embedded matrix (columns are cells)
      feature_matrix_norms <- sqrt(Matrix::colSums(feature_matrix_for_similarity_to_use * feature_matrix_for_similarity_to_use))
      rm(dense_mat, sv); gc() # Clean up intermediate objects
    } else {
      # If dimensions don't match, revert to using the original matrix and log a warning.
      warning("Dimension mismatch during PCA embedding: Number of columns in embedded matrix (", ncol(feature_matrix_embed), ") does not match number of cells (", length(cell_ids), "). Using original feature matrix for similarity calculations.")
      # feature_matrix_for_similarity_to_use remains the original matrix
      rm(dense_mat, sv); gc()
    }
  }

  # --- Symmetrize the feature_graph using averaging ---
  message("Symmetrizing feature_graph using averaging...")
  feature_graph_sym <- (feature_graph + Matrix::t(feature_graph)) / 2
  diag(feature_graph_sym) <- 0 # Remove self-loops

  gs_check <- Matrix::summary(feature_graph_sym)
  if (nrow(gs_check) == 0) {
    warning("Symmetrized feature_graph contains no edges. Metrics 1 & 3 will be NA. Plots for MSD-FGN and MFS-FGN will likely be white.")
  } else {
    message("Symmetrized graph has ", nrow(gs_check), " edges (after removing self-loops).")
  }

  # --- Prepare vector for Metric 1 & 3 (vectorized) ---
  median_spatial_dist_of_graph_neighbors_vec <- rep(NA_real_, num_cells) # Metric 1
  median_graph_sim_of_graph_neighbors_vec <- rep(NA_real_, num_cells)   # Metric 3

  gs <- data.table::as.data.table(Matrix::summary(feature_graph_sym))
  data.table::setnames(gs, c("i","j","x"))
  setnames(gs, c("i","j"), c("query_cell_idx", "neighbor_cell_idx"))
  gs <- gs[query_cell_idx != neighbor_cell_idx] # Safeguard

  if (nrow(gs) > 0) {
    message("  Calculating Metric 1 (MSD-FGN) and Metric 3 (MFS-FGN)...")
    # Calculate spatial distances for all edges using NORMALIZED coordinates
    x1 <- spatial_coords_mat_normalized[gs$query_cell_idx, 1]
    y1 <- spatial_coords_mat_normalized[gs$query_cell_idx, 2]
    x2 <- spatial_coords_mat_normalized[gs$neighbor_cell_idx, 1]
    y2 <- spatial_coords_mat_normalized[gs$neighbor_cell_idx, 2]
    gs[, spatial_distance := sqrt((x1 - x2)^2 + (y1 - y2)^2)]

    # Ensure k_neighbors doesn't exceed number of available neighbors for a cell
    # Rank neighbors globally, then select the subset corresponding to neighbor_rank_range
    graph_neighbors_ranked_all <- gs[, .SD[order(-x)], by = query_cell_idx]
    # For each query cell, take rows from neighbor_rank_range[1]..neighbor_rank_range[2] (bounded by available neighbors)
    graph_neighbors_ranked <- graph_neighbors_ranked_all[, {
      start_r <- neighbor_rank_range[1]
      end_r   <- neighbor_rank_range[2]
      available <- .N
      if (available < start_r) {
        # no neighbors in this rank range -> return zero rows
        .SD[0]
      } else {
        end_sel <- min(end_r, available)
        .SD[start_r:end_sel]
      }
    }, by = query_cell_idx]


    if (nrow(graph_neighbors_ranked) > 0) {
      median_spatial_results <- graph_neighbors_ranked[, .(median_dist = median(spatial_distance, na.rm = TRUE)), by = query_cell_idx]
      median_graph_sim_results_of_graph_nn <- graph_neighbors_ranked[, .(median_graph_sim = median(x, na.rm = TRUE)), by = query_cell_idx]

      median_spatial_dist_of_graph_neighbors_vec[median_spatial_results$query_cell_idx] <- median_spatial_results$median_dist
      median_graph_sim_of_graph_neighbors_vec[median_graph_sim_results_of_graph_nn$query_cell_idx] <- median_graph_sim_results_of_graph_nn$median_graph_sim
    } else {
      message("Warning: No graph neighbors found for any cell to rank. Metrics 1 & 3 will be NA.")
    }
  } else {
    warning("Symmetrized feature_graph contains no edges (after removing self-loops). Metrics 1 & 3 will be NA.")
  }

  # --- Prepare spatial nearest neighbors (RANN) ---
  all_spatial_nn_results <- RANN::nn2(data = spatial_coords_mat_normalized, k = k_neighbors + 1)
  all_nn_idx <- all_spatial_nn_results$nn.idx

  # Pre-allocate result vectors for metrics calculated in batches
  median_graph_sim_of_spatial_neighbors_vec <- rep(NA_real_, num_cells) # Metric 2
  # Metric 4 (median_graph_dist_of_spatial_neighbors) is REMOVED.

  # --- Batch processing function (per-batch) ---
  # This function calculates Metric 2 (MFS-SN)
  process_batch_fast <- function(sources) {
    nb <- length(sources)
    res_metric2 <- numeric(nb) # Median Feature Similarity to Spatial Neighbors (MFS-SN)

    for (ii in seq_along(sources)) {
      i <- sources[ii] # Current cell index

      # Get spatial neighbors for cell 'i'
      # use the same neighbor_rank_range to select spatial neighbors (by rank)
      raw_targets_full <- all_nn_idx[i, ]
      raw_targets_full <- raw_targets_full[raw_targets_full != i & raw_targets_full > 0]
      # select neighbors by the requested rank range
      sel_start <- neighbor_rank_range[1]
      sel_end   <- neighbor_rank_range[2]
      if (length(raw_targets_full) < sel_start) {
        raw_targets <- integer(0)
      } else {
        raw_targets <- raw_targets_full[sel_start:min(sel_end, length(raw_targets_full))]
      }

      if (length(raw_targets) == 0) {
        res_metric2[ii] <- NA_real_
        next # Skip if no spatial neighbors found
      }

      # --- Metric 2: Median Feature Similarity to Spatial Neighbors (MFS-SN) ---

      # Get the cell NAMES for the spatial neighbors using the indices from raw_targets
      neighbor_cell_names <- cell_ids[raw_targets]

      # Get the column names of the matrix we are indexing into
      matrix_colnames <- colnames(feature_matrix_for_similarity_to_use)

      # Find the NUMERICAL COLUMN INDICES in matrix_colnames that match neighbor_cell_names.
      # This is crucial for robustly avoiding "subscript out of bounds" if name resolution is faulty.
      neighbor_col_indices_in_matrix <- match(neighbor_cell_names, matrix_colnames)

      # Check if any names were not found. This would result in NA in neighbor_col_indices_in_matrix.
      # If NAs are present, we need to filter them out to avoid errors.
      if (any(is.na(neighbor_col_indices_in_matrix))) {
        # Identify which names were not found and their corresponding original raw_targets indices.
        missing_names_mask <- is.na(neighbor_col_indices_in_matrix)
        missing_names <- neighbor_cell_names[missing_names_mask]
        corresponding_raw_targets <- raw_targets[missing_names_mask]

        warning("Cell names from spatial neighbors (indices: ",
                paste(corresponding_raw_targets, collapse=", "),
                "; names: ", paste(missing_names, collapse=", "),
                ") not found in column names of feature_matrix_for_similarity_to_use. Setting metrics to NA for these neighbors.")

        # Filter out NA indices to only use valid numerical indices for matrix subsetting.
        valid_neighbor_col_indices <- na.omit(neighbor_col_indices_in_matrix)

        if (length(valid_neighbor_col_indices) == 0) {
          # If no valid neighbors could be found after filtering, return NA for metrics
          res_metric2[ii] <- NA_real_
          next # Skip to the next source cell
        }
        # Use the valid numerical column indices for sparse matrix column selection.
        mat2 <- feature_matrix_for_similarity_to_use[, valid_neighbor_col_indices, drop = FALSE]

      } else {
        # All neighbor cell names were found as column names.
        # Use the numerical indices derived from match() to subset the matrix.
        mat2 <- feature_matrix_for_similarity_to_use[, neighbor_col_indices_in_matrix, drop = FALSE]
        valid_neighbor_col_indices = neighbor_col_indices_in_matrix
      }

      # Indexing for the query cell `i` (which is a numerical index into cell_ids and columns of the matrix)
      vec1 <- feature_matrix_for_similarity_to_use[, i, drop = FALSE]

      # Calculate dot products and cosine similarity
      dot_products <- as.vector(crossprod(vec1, mat2))

      # Norms for the query cell `i`
      norm1 <- feature_matrix_norms[i] # `i` is a 1-based numerical index

      # Norms for the neighbor cells (using the valid numerical column indices)
      norm2 <- feature_matrix_norms[valid_neighbor_col_indices] # Use indices found by match

      denom <- norm1 * norm2
      cos_sim <- ifelse(denom == 0, 0, dot_products / denom)

      # Clip cosine similarity values to be within [0, 1] range.
      # This addresses cases where values might exceed 1 due to numerical precision or unusual data characteristics.
      cos_sim <- pmax(0, pmin(1, cos_sim))

      res_metric2[ii] <- median(cos_sim, na.rm = TRUE)

      # --- Metric 4 (median_graph_dist_of_spatial_neighbors) is REMOVED. ---
    }
    list(idx = sources, metric2_med = res_metric2)
  }

  # --- Run batches in parallel ---
  seq_starts <- seq(1L, num_cells, by = batch_size)
  batch_idx_list <- lapply(seq_starts, function(s) {
    e <- min(s + batch_size - 1L, num_cells)
    seq.int(s, e)
  })

  message("  Calculating Metric 2 (MFS-SN) in parallel batches...")
  batch_results <- future.apply::future_lapply(batch_idx_list, FUN = process_batch_fast)

  # Collect results into the preallocated vectors
  for (br in batch_results) {
    median_graph_sim_of_spatial_neighbors_vec[br$idx] <- br$metric2_med
  }

  # --- Assemble coherence dataframe ---
  # Define which metrics to include in the coherence_df and for availability count
  metrics_to_include_in_df <- c(
    "median_spatial_dist_of_graph_neighbors", # Metric 1
    "median_graph_sim_of_spatial_neighbors",     # Metric 2
    "median_graph_sim_of_graph_neighbors"        # Metric 3
  )

  coherence_df <- data.frame(
    cell_id = cell_ids,
    modality = modality_name,
    median_spatial_dist_of_graph_neighbors = median_spatial_dist_of_graph_neighbors_vec,
    median_graph_sim_of_spatial_neighbors = median_graph_sim_of_spatial_neighbors_vec,
    median_graph_sim_of_graph_neighbors = median_graph_sim_of_graph_neighbors_vec
    # metric 4 (median_graph_dist_of_spatial_neighbors) is REMOVED
  )
  data.table::setDT(coherence_df)

  # Calculate number of available metrics per cell, excluding the removed one
  coherence_df$num_metrics_available <- rowSums(!is.na(coherence_df[, metrics_to_include_in_df, with = FALSE]))
  message("Per-cell metrics availability summary:")
  print(table(coherence_df$num_metrics_available))

  # --- Remove cells with NA in any of the core metrics ---
  initial_nrow_coherence_df <- nrow(coherence_df)
  coherence_df <- na.omit(coherence_df)

  if (nrow(coherence_df) < initial_nrow_coherence_df) {
    message(paste0("Info: Removed ", initial_nrow_coherence_df - nrow(coherence_df), " cells with NA coherence metrics via na.omit()."))
  }

  # Return empty plots if no data remains after NA removal
  if (nrow(coherence_df) == 0) {
    message("CRITICAL WARNING: coherence_df is empty after removing NAs. All plots will be empty.")
    empty_plots <- function() ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data") +
      ggplot2::theme_void()
    return(list(
      coherence_df = data.frame(
        cell_id = character(0),
        modality = character(0),
        median_spatial_dist_of_graph_neighbors = numeric(0),
        median_graph_sim_of_spatial_neighbors = numeric(0),
        median_graph_sim_of_graph_neighbors = numeric(0),
        combined_coherence_score = numeric(0)
      ),
      plot_dist_spatial_for_graph_nn = empty_plots(),
      plot_sim_graph_for_spatial_nn = empty_plots(),
      plot_sim_graph_for_graph_nn = empty_plots(),
      plot_spatial_coherence_metrics = empty_plots(),
      diagnostics = list(note = "Empty coherence_df after NA removal")
    ))
  }


  # --- Calculate Combined Coherence Score (Metric 5) ---
  # New logic combines Metric 1 (Distance) and Metric 2 (Similarity)
  coherence_df$combined_coherence_score <- NA_real_ # Initialize

  # Define which metrics to use for the new combined score
  metric1_col <- "median_spatial_dist_of_graph_neighbors" # Distance metric (needs inversion)
  metric2_col <- "median_graph_sim_of_spatial_neighbors" # Similarity metric (no inversion)

  # Find cells where both required metrics are available
  valid_idx_for_combined <- which(!is.na(coherence_df[[metric1_col]]) & !is.na(coherence_df[[metric2_col]]))

  if (length(valid_idx_for_combined) > 0) {
    metric1_vals <- coherence_df[[metric1_col]][valid_idx_for_combined]
    metric2_vals <- coherence_df[[metric2_col]][valid_idx_for_combined]

    # Normalize and invert Metric 1 (distance): lower distance = higher score
    min_m1 <- min(metric1_vals, na.rm = TRUE)
    max_m1 <- max(metric1_vals, na.rm = TRUE)
    norm_m1 <- if (max_m1 == min_m1) rep(0.5, length(metric1_vals)) else 1 - (metric1_vals - min_m1) / (max_m1 - min_m1)

    # Normalize Metric 2 (similarity): higher similarity = higher score (no inversion needed)
    min_m2 <- min(metric2_vals, na.rm = TRUE)
    max_m2 <- max(metric2_vals, na.rm = TRUE)
    norm_m2 <- if (max_m2 == min_m2) rep(0.5, length(metric2_vals)) else (metric2_vals - min_m2) / (max_m2 - min_m2)

    # Average the normalized scores
    coherence_df$combined_coherence_score[valid_idx_for_combined] <- (norm_m1 + norm_m2) / 2
  }

  # create human-readable neighbor range label for titles
  neighbor_range_label <- paste0("(NN: ", neighbor_rank_range[1], ":", neighbor_rank_range[2], ")")

  # --- Create density plots ---
  plot_dist_spatial_for_graph_nn <- ggplot2::ggplot(coherence_df, ggplot2::aes(x = median_spatial_dist_of_graph_neighbors)) +
    ggplot2::geom_density(fill = "lightblue", alpha = 0.7, na.rm = TRUE) +
    ggplot2::labs(title = paste0(modality_name, ": Spatial Dist to Feature-Graph ", neighbor_range_label), x = "Median Spatial Distance", y = "Density") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  plot_sim_graph_for_spatial_nn <- ggplot2::ggplot(coherence_df, ggplot2::aes(x = median_graph_sim_of_spatial_neighbors)) +
    ggplot2::geom_density(fill = "lightcoral", alpha = 0.7, na.rm = TRUE) +
    ggplot2::labs(title = paste0(modality_name, ": Feature Sim to Spatial  ", neighbor_range_label), x = "Median Feature Similarity", y = "Density") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  plot_sim_graph_for_graph_nn <- ggplot2::ggplot(coherence_df, ggplot2::aes(x = median_graph_sim_of_graph_neighbors)) +
    ggplot2::geom_density(fill = "lightgreen", alpha = 0.7, na.rm = TRUE) +
    ggplot2::labs(title = paste0(modality_name, ": Feature Sim to Feature-Graph  ", neighbor_range_label), x = "Median Feature Similarity", y = "Density") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # plot_dist_graph_for_spatial_nn is REMOVED from the output list.

  # --- Spatial plot: pivot long ---
  spatial_coords_df_for_merge <- as.data.frame(spatial_coords_df)
  if (!"cell_id" %in% colnames(spatial_coords_df_for_merge)) spatial_coords_df_for_merge$cell_id <- rownames(spatial_coords_df_for_merge)

  # Default empty plot if no data is available
  p_spatial_metrics <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data for spatial plot") + ggplot2::theme_void()

  if (nrow(coherence_df) > 0) {
    plot_data_spatial <- merge(spatial_coords_df_for_merge, coherence_df, by = "cell_id", all.x = TRUE)

    # Define which metrics to include in the spatial plot
    spatial_plot_metrics_cols <- c(
      "median_spatial_dist_of_graph_neighbors", # Metric 1
      "median_graph_sim_of_spatial_neighbors",     # Metric 2
      "combined_coherence_score"                   # New CCS
    )

    # Define labels only for the included metrics
    metric_labels_for_spatial_plot <- c(
      median_spatial_dist_of_graph_neighbors = "Spatial Dist. to Feature-Graph NN",
      median_graph_sim_of_spatial_neighbors = "Feature Sim. to Spatial NN",
      combined_coherence_score = "Combined Coherence Score"
    )

    plot_data_long <- plot_data_spatial %>%
      tidyr::pivot_longer(cols = all_of(spatial_plot_metrics_cols),
                          names_to = "metric_type", values_to = "metric_value")

    # Ensure the factor levels match the desired order for facets
    plot_data_long$metric_type <- factor(plot_data_long$metric_type,
                                         levels = spatial_plot_metrics_cols,
                                         labels = metric_labels_for_spatial_plot)

    plot_data_long_normalized <- plot_data_long %>%
      dplyr::group_by(metric_type) %>%
      dplyr::mutate(
        val_min = suppressWarnings(min(metric_value, na.rm = TRUE)),
        val_max = suppressWarnings(max(metric_value, na.rm = TRUE)),
        norm_val = dplyr::case_when(
          val_max == val_min ~ rep(0.5, n()), # Handle constant values
          TRUE ~ (metric_value - val_min) / (val_max - val_min)
        ),
        # For distance metrics (Metric 1), higher normalized value means closer (1-norm_val)
        # For similarity (Metric 2) and score (CCS) metrics, higher normalized value is better (norm_val)
        plot_color_value = dplyr::case_when(
          metric_type == "Spatial Dist. to Feature-Graph NN" ~ 1 - norm_val,
          TRUE ~ norm_val
        )
      ) %>%
      dplyr::ungroup()

    if (nrow(plot_data_long_normalized) > 0 && !all(is.na(plot_data_long_normalized$plot_color_value))) {
      p_spatial_metrics <- ggplot2::ggplot(plot_data_long_normalized, ggplot2::aes(x = x, y = y, color = plot_color_value)) +
        ggplot2::geom_point(size = 0.1, alpha = 0.6) +
        ggplot2::facet_wrap(~metric_type, ncol = 3) +
        ggplot2::scale_color_viridis_c(option = "viridis", direction = 1, na.value = "grey80", name = "Coherence (Normalized)") +
        ggplot2::labs(title = paste0(modality_name, ": Spatial Coherence Metrics (k=", k_neighbors, ")"), x = "Spatial X", y = "Spatial Y") +
        ggplot2::coord_fixed() + ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "bottom")
    } else {
      message("Warning: Spatial plot data (after normalization/NA handling) is empty or all NA for coloring. Returning empty spatial plot for this modality.")
    }
  }

  # --- Diagnostics ---
  diagnostics <- list(
    n_cells = num_cells,
    k_neighbors = k_neighbors,
    neighbor_rank_range = neighbor_rank_range,
    n_workers = n_workers,
    batch_size = batch_size,
    used_pca_embedding = exists("feature_matrix_embed"),
    pca_dims = if (exists("feature_matrix_embed")) nrow(feature_matrix_embed) else NA_integer_
  )

  # --- Return ---
  out_list <- list(
    coherence_df = coherence_df,
    plot_dist_spatial_for_graph_nn = plot_dist_spatial_for_graph_nn,
    plot_sim_graph_for_spatial_nn = plot_sim_graph_for_spatial_nn,
    plot_sim_graph_for_graph_nn = plot_sim_graph_for_graph_nn,
    # plot_dist_graph_for_spatial_nn is REMOVED from the list
    plot_spatial_coherence_metrics = p_spatial_metrics,
    diagnostics = diagnostics
  )
  return(out_list)
}


#' Validate Graph Coherence Against Random Cell Connections
#'
#' \code{validateGraphCoherenceWithRandom()} assesses the spatial coherence of a feature graph
#' by comparing observed metrics for top-k neighbors against random cell sets (averaged).
#' Returns per-cell metrics and a combined plot comparing distributions.
#'
#' @param feature_graph A sparse matrix (cells x cells) representing feature-based similarities. Rownames must be cell IDs.
#' @param spatial_coords_df A data.frame (cells x 2) of spatial x, y coordinates. Rownames must be cell IDs and match `feature_graph`.
#' @param k_neighbors Integer, the number of neighbors to consider for comparison (both for graph and random).
#' @param num_random_comparisons Integer, number of random neighbor sets to generate for averaging.
#' @param modality_name Character string (e.g., "RNA", "ATAC") for labeling results and plots.
#' @param n_workers Integer, number of parallel workers. NULL defaults to detectCores() - 1.
#' @return A list containing coherence metrics dataframe and plots comparing graph vs. random.
#' @export
validateGraphCoherenceWithRandom <- function(
    feature_graph,
    spatial_coords_df,
    k_neighbors = 10,
    num_random_comparisons = 5,
    modality_name = "Modality",
    n_workers = NULL
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install 'Matrix'.")
  if (!requireNamespace("RANN", quietly = TRUE)) stop("Install 'RANN'.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install 'data.table'.")
  if (!requireNamespace("future.apply", quietly = TRUE)) stop("Install 'future.apply'.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install 'ggplot2'.")
  if (!requireNamespace("cowplot", quietly = TRUE)) stop("Install 'cowplot'.")

  if (!is(feature_graph, "Matrix")) stop("`feature_graph` must be a sparse Matrix (e.g., dgCMatrix).")
  if (!is.data.frame(spatial_coords_df) || ncol(spatial_coords_df) < 2) stop("`spatial_coords_df` must be a data.frame with at least two columns (x,y).")
  if (is.null(rownames(feature_graph)) || is.null(rownames(spatial_coords_df))) {
    stop("`feature_graph` and `spatial_coords_df` must have matching cell IDs (rownames).")
  }

  # Align cell sets
  graph_cell_ids <- rownames(feature_graph)
  spatial_cell_ids <- rownames(spatial_coords_df)
  common_cells <- intersect(graph_cell_ids, spatial_cell_ids)
  if (length(common_cells) == 0) stop("No common cell IDs found across inputs.")
  if (length(common_cells) < nrow(feature_graph) || length(common_cells) < nrow(spatial_coords_df)) {
    message("Warning: Subsetting to ", length(common_cells), " common cells.")
  }
  feature_graph <- feature_graph[common_cells, common_cells, drop = FALSE]
  spatial_coords_df <- spatial_coords_df[common_cells, , drop = FALSE]
  num_cells <- nrow(feature_graph)
  cell_ids <- rownames(feature_graph)

  if (num_cells <= 1) {
    empty_plots <- function() ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data") + ggplot2::theme_void()
    return(list(
      coherence_df = data.frame(cell_id = character(0), modality = character(0), graph_type = character(0)),
      plot_comparison = empty_plots(),
      diagnostics = list(note = "Too few cells")
    ))
  }

  if (k_neighbors >= num_cells) {
    warning("k_neighbors >= num_cells; adjusting to num_cells - 1")
    k_neighbors <- max(1, num_cells - 1)
  }
  if (k_neighbors == 0) stop("k_neighbors must be at least 1 for comparison.")

  # Parallel plan
  if (is.null(n_workers)) {
    n_workers <- max(1L, parallel::detectCores() - 1L)
  } else {
    n_workers <- as.integer(n_workers); if (n_workers < 1L) n_workers <- 1L
  }
  old_plan <- future::plan()
  options(future.globals.maxSize = +Inf)
  future::plan(future::multisession, workers = n_workers)
  on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)

  # Normalize spatial coordinates to 0..1 by max dimension range (same approach as assessSpatialGraphCoherenceOverall)
  spatial_coords_mat_raw <- as.matrix(spatial_coords_df[, c("x", "y")])
  x_range <- range(spatial_coords_mat_raw[,1])
  y_range <- range(spatial_coords_mat_raw[,2])
  max_dim_range <- max(diff(x_range), diff(y_range))
  if (max_dim_range == 0) {
    spatial_coords_mat_normalized <- spatial_coords_mat_raw
  } else {
    spatial_coords_mat_normalized <- sweep(spatial_coords_mat_raw, 2, c(min(x_range), min(y_range)), "-")
    spatial_coords_mat_normalized <- spatial_coords_mat_normalized / max_dim_range
  }

  # Symmetrize feature_graph and remove self-loops
  feature_graph_sym <- (feature_graph + Matrix::t(feature_graph)) / 2
  diag(feature_graph_sym) <- 0

  # Helper to compute per-cell median spatial distance and median graph similarity for given neighbor sets
  calculate_coherence_metrics <- function(cell_indices_to_process, graph_matrix, spatial_matrix_norm, k, random_neighbor_set = NULL) {
    nb <- length(cell_indices_to_process)
    median_spatial_dist_vec <- rep(NA_real_, nb)
    median_graph_sim_vec <- rep(NA_real_, nb)

    # If graph is empty, return NA vectors
    if (nrow(graph_matrix) == 0 || Matrix::nnzero(graph_matrix) == 0) {
      return(list(med_dist = median_spatial_dist_vec, med_sim = median_graph_sim_vec))
    }

    is_random <- !is.null(random_neighbor_set)

    if (!is_random) {
      # Use graph-based neighbors: rank neighbors by similarity and pick top-k
      gs <- data.table::as.data.table(Matrix::summary(graph_matrix))
      if (nrow(gs) == 0) return(list(med_dist = median_spatial_dist_vec, med_sim = median_graph_sim_vec))
      data.table::setnames(gs, c("i","j","x"))
      gs <- gs[i != j]
      # compute normalized spatial distances for edges
      x1 <- spatial_matrix_norm[gs$i, 1]; y1 <- spatial_matrix_norm[gs$i, 2]
      x2 <- spatial_matrix_norm[gs$j, 1]; y2 <- spatial_matrix_norm[gs$j, 2]
      gs[, spatial_distance := sqrt((x1 - x2)^2 + (y1 - y2)^2)]
      graph_neighbors_ranked <- gs[, {
        current_k <- min(k, .N)
        if (current_k > 0) head(.SD[order(-x)], current_k) else .SD[0]
      }, by = i]
      if (nrow(graph_neighbors_ranked) == 0) return(list(med_dist = median_spatial_dist_vec, med_sim = median_graph_sim_vec))
      median_spatial_results <- graph_neighbors_ranked[, .(median_dist = median(spatial_distance, na.rm = TRUE)), by = i]
      median_graph_sim_results <- graph_neighbors_ranked[, .(median_sim = median(x, na.rm = TRUE)), by = i]
      median_spatial_dist_vec[median_spatial_results$i] <- median_spatial_results$median_dist
      median_graph_sim_vec[median_graph_sim_results$i] <- median_graph_sim_results$median_sim
      return(list(med_dist = median_spatial_dist_vec, med_sim = median_graph_sim_vec))
    } else {
      # random_neighbor_set is a list indexed by cell index (1...num_cells) containing neighbor indices
      for (ii in seq_along(cell_indices_to_process)) {
        cell_idx <- cell_indices_to_process[ii]
        rand_neigh <- random_neighbor_set[[cell_idx]]
        if (length(rand_neigh) < k) next
        x1 <- spatial_matrix_norm[cell_idx, 1]; y1 <- spatial_matrix_norm[cell_idx, 2]
        x2 <- spatial_matrix_norm[rand_neigh, 1]; y2 <- spatial_matrix_norm[rand_neigh, 2]
        dists <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
        median_spatial_dist_vec[ii] <- median(dists, na.rm = TRUE)
        # similarity from graph_matrix: row cell_idx and columns rand_neigh
        sims <- as.vector(graph_matrix[cell_idx, rand_neigh])
        median_graph_sim_vec[ii] <- median(sims, na.rm = TRUE)
      }
      return(list(med_dist = median_spatial_dist_vec, med_sim = median_graph_sim_vec))
    }
  }

  # Graph-based metrics
  message("  Calculating metrics for graph neighbors...")
  graph_metrics <- calculate_coherence_metrics(
    cell_indices_to_process = seq_len(num_cells),
    graph_matrix = feature_graph_sym,
    spatial_matrix_norm = spatial_coords_mat_normalized,
    k = k_neighbors,
    random_neighbor_set = NULL
  )

  # Random-based metrics averaged over multiple random comparisons
  message("  Generating and computing random neighbor comparisons (n=", num_random_comparisons, ")...")
  random_med_dist_matrix <- matrix(NA_real_, nrow = num_cells, ncol = num_random_comparisons)
  random_med_sim_matrix  <- matrix(NA_real_, nrow = num_cells, ncol = num_random_comparisons)

  # Use future.apply for parallel repeats if available
  random_results <- future.apply::future_lapply(seq_len(num_random_comparisons), function(rc) {
    set.seed(100 + rc)
    current_random_sets <- vector("list", num_cells)
    all_idx <- seq_len(num_cells)
    for (i in seq_len(num_cells)) {
      pool <- all_idx[all_idx != i]
      if (length(pool) == 0) {
        current_random_sets[[i]] <- integer(0)
      } else if (length(pool) <= k_neighbors) {
        current_random_sets[[i]] <- sample(pool, length(pool))
      } else {
        current_random_sets[[i]] <- sample(pool, k_neighbors)
      }
    }
    calculate_coherence_metrics(
      cell_indices_to_process = seq_len(num_cells),
      graph_matrix = feature_graph_sym,
      spatial_matrix_norm = spatial_coords_mat_normalized,
      k = k_neighbors,
      random_neighbor_set = current_random_sets
    )
  })

  for (i in seq_along(random_results)) {
    random_med_dist_matrix[, i] <- random_results[[i]]$med_dist
    random_med_sim_matrix[, i]  <- random_results[[i]]$med_sim
  }
  avg_random_med_dist <- rowMeans(random_med_dist_matrix, na.rm = TRUE)
  avg_random_med_sim  <- rowMeans(random_med_sim_matrix, na.rm = TRUE)

  # Assemble data.frames
  coherence_df_graph <- data.frame(
    cell_id = cell_ids,
    modality = modality_name,
    graph_type = "Graph Neighbors",
    median_spatial_dist = graph_metrics$med_dist,
    median_graph_sim = graph_metrics$med_sim,
    stringsAsFactors = FALSE
  )
  coherence_df_random <- data.frame(
    cell_id = cell_ids,
    modality = modality_name,
    graph_type = "Random Cells",
    median_spatial_dist = avg_random_med_dist,
    median_graph_sim = avg_random_med_sim,
    stringsAsFactors = FALSE
  )

  full_coherence_df <- rbind(coherence_df_graph, coherence_df_random)
  full_coherence_df <- na.omit(full_coherence_df)
  if (nrow(full_coherence_df) == 0) {
    empty_plots <- function() ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data") + ggplot2::theme_void()
    return(list(coherence_df = data.frame(), plot_comparison = empty_plots(), diagnostics = list(note = "Empty after NA removal")))
  }

  # Create comparison plots (density overlays for spatial dist and graph sim)
  plot_med_dist <- ggplot2::ggplot(full_coherence_df, ggplot2::aes(x = median_spatial_dist, fill = graph_type)) +
    ggplot2::geom_density(alpha = 0.6, na.rm = TRUE) +
    ggplot2::labs(title = paste0(modality_name, ": Median Spatial Distance to Neighbors (k=", k_neighbors, ")"),
                  x = "Median Spatial Distance (Lower = Better Coherence)", y = "Density", fill = "Neighbor Type") +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  plot_med_sim <- ggplot2::ggplot(full_coherence_df, ggplot2::aes(x = median_graph_sim, fill = graph_type)) +
    ggplot2::geom_density(alpha = 0.6, na.rm = TRUE) +
    ggplot2::labs(title = paste0(modality_name, ": Median Graph Similarity to Neighbors (k=", k_neighbors, ")"),
                  x = "Median Graph Similarity (Higher = Better Coherence)", y = "Density", fill = "Neighbor Type") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::theme_minimal() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  plot_comparison <- cowplot::plot_grid(plot_med_dist, plot_med_sim, ncol = 1, align = "v")

  diagnostics <- list(
    n_cells = num_cells,
    k_neighbors = k_neighbors,
    num_random_comparisons = num_random_comparisons,
    n_workers = n_workers
  )

  return(list(coherence_df = full_coherence_df, plot_comparison = plot_comparison, diagnostics = diagnostics))
}


# Create labelEdges for ATAC-derived annotations (e.g., TF motif sets, ChIP peaks, bulk pREs)
# Inputs:
#  - region_table: data.frame with at least seqnames,start,end and a column 'cluster' naming annotation (label)
#  - peaks_gr: GRanges of ATAC peaks (must match columns of atac_counts_matrix_filtered)
#  - ATACMat_cells_by_peak: cell x peak sparse matrix (rows = cells, cols = peaks)
# Returns: a cell x label matrix (rows cells, cols labels)
createLabelEdgesFromRegions <- function(region_table, peaks_gr, ATACMat_cells_by_peak) {
  if (missing(region_table) || missing(peaks_gr) || missing(ATACMat_cells_by_peak)) {
    stop("Must provide region_table, peaks_gr (GRanges) and ATACMat (cells x peaks).")
  }
  if (!is(peaks_gr, "GRanges")) stop("peaks_gr must be a GRanges.")
  if (! (is(ATACMat_cells_by_peak, "matrix") || is(ATACMat_cells_by_peak, "Matrix"))) stop("ATACMat must be a matrix/Matrix with dimensions cells x peaks.")
  if (ncol(ATACMat_cells_by_peak) != length(peaks_gr)) stop("Number of peaks must match columns of ATACMat.")

  # region_table must contain seqnames,start,end and cluster columns
  if (!all(c("seqnames","start","end","cluster") %in% colnames(region_table))) {
    stop("region_table must contain columns: seqnames, start, end, cluster")
  }

  # Convert region_table to GRanges
  regs <- GRanges(seqnames = region_table$seqnames,
                  ranges = IRanges(start = as.integer(region_table$start), end = as.integer(region_table$end)))
  regs$cluster <- as.character(region_table$cluster)
  regs$region_id <- seq_len(length(regs))

  # find overlaps between peaks and regions
  ov <- findOverlaps(peaks_gr, regs)
  if (length(ov) == 0) {
    stop("No overlaps between provided regions and peaks.")
  }
  # For each region (cluster), get list of peaks indices
  df_ov <- data.frame(peak_idx = queryHits(ov), region_idx = subjectHits(ov), cluster = regs$cluster[subjectHits(ov)], stringsAsFactors = FALSE)

  # For each cluster (label), identify unique peaks
  cluster_list <- split(df_ov$peak_idx, df_ov$cluster)
  cluster_list <- lapply(cluster_list, unique)

  # Build labelEdges: for each cluster compute per-cell accessibility score (mean or sum)
  label_edges_mat <- sapply(names(cluster_list), function(cl) {
    peaks_idx <- cluster_list[[cl]]
    if (length(peaks_idx) == 0) return(rep(0, nrow(ATACMat_cells_by_peak)))
    # sum across peaks per cell, then normalize by total peak counts per cell to adjust for depth
    per_cell_sum <- Matrix::rowSums(ATACMat_cells_by_peak[, peaks_idx, drop = FALSE])
    depth <- Matrix::rowSums(ATACMat_cells_by_peak)
    depth[depth == 0] <- 1
    per_cell_sum / depth
  }, simplify = "matrix")
  # Ensure rownames / colnames
  if (!is.null(rownames(ATACMat_cells_by_peak))) rownames(label_edges_mat) <- rownames(ATACMat_cells_by_peak)
  colnames(label_edges_mat) <- names(cluster_list)
  label_edges_mat[label_edges_mat < 0] <- 0
  return(as.matrix(label_edges_mat))
}



# Enhanced: Robust ATAC label edge construction combining Options 1 & 4
buildAtacLabelEdgesFromGeneActivity <- function(
    atac_gene_activity_matrix,  # cells x genes
    n_clusters = 30,             # Increased from 15
    min_cluster_size = 50,       # Minimum cells per cluster
    n_init = 50,                 # Multiple k-means initializations
    min_weight = 0.05,           # Filter very weak edges
    pca_dims = 30,               # Number of PCs to use
    balance_clusters = TRUE      # Apply Option 2 normalization
) {
  library(Matrix)

  if (!is.matrix(atac_gene_activity_matrix) && !is(atac_gene_activity_matrix, "Matrix")) {
    stop("Provide matrix-like atac_gene_activity_matrix")
  }

  message("Building ATAC label edges with robust k-means...")
  message("  Parameters: n_clusters=", n_clusters, ", n_init=", n_init,
          ", min_weight=", min_weight)

  mat <- as.matrix(atac_gene_activity_matrix)  # cells x genes
  n_cells <- nrow(mat)
  n_genes <- ncol(mat)

  if (n_cells < n_clusters) {
    warning("Fewer cells than clusters requested. Reducing n_clusters to ", n_cells)
    n_clusters <- n_cells
  }

  # ========================================================================
  # STEP 1: PCA on cells (not genes) for cell clustering
  # ========================================================================
  message("  Step 1: Computing PCA on cell x gene matrix...")

  # Normalize and scale for PCA
  mat_scaled <- scale(mat)
  mat_scaled[is.na(mat_scaled)] <- 0

  pca_res <- tryCatch({
    stats::prcomp(mat_scaled, rank. = min(pca_dims, ncol(mat_scaled)))
  }, error = function(e) {
    message("  PCA failed, using raw matrix")
    return(list(x = mat_scaled))
  })

  pca_coords <- pca_res$x[, 1:min(pca_dims, ncol(pca_res$x)), drop = FALSE]
  message("  Using ", ncol(pca_coords), " PCs for clustering")


  # ========================================================================
  # STEP 2: Robust k-means with multiple initializations
  # ========================================================================
  message("  Step 2: Running k-means with ", n_init, " initializations...")

  best_kmeans <- NULL
  best_wss <- Inf

  set.seed(42)
  for (i in 1:n_init) {
    km <- tryCatch({
      stats::kmeans(pca_coords, centers = n_clusters, nstart = 1,
                    iter.max = 100, algorithm = "Lloyd")
    }, error = function(e) NULL)

    if (!is.null(km) && km$tot.withinss < best_wss) {
      best_wss <- km$tot.withinss
      best_kmeans <- km
    }

    if (i %% 10 == 0) message("    Completed ", i, "/", n_init, " initializations")
  }

  if (is.null(best_kmeans)) {
    stop("K-means clustering failed completely")
  }

  message("  Best within-cluster sum of squares: ", round(best_wss, 2))

  # Check cluster sizes
  cluster_sizes <- table(best_kmeans$cluster)
  message("  Cluster sizes: ", paste(cluster_sizes, collapse = ", "))

  # Merge small clusters into nearest larger cluster
  small_clusters <- which(cluster_sizes < min_cluster_size)
  if (length(small_clusters) > 0) {
    message("  Merging ", length(small_clusters), " small clusters (< ",
            min_cluster_size, " cells)")

    cluster_assignments <- best_kmeans$cluster
    centroids <- best_kmeans$centers

    for (small_clust in small_clusters) {
      # Find nearest large cluster
      small_centroid <- centroids[small_clust, ]
      distances_to_others <- apply(centroids[-small_clust, , drop = FALSE], 1, function(x) {
        sqrt(sum((x - small_centroid)^2))
      })
      nearest_large <- which.min(distances_to_others)

      # Reassign cells
      cluster_assignments[cluster_assignments == small_clust] <- nearest_large
    }

    # Update cluster assignments
    best_kmeans$cluster <- cluster_assignments

    # Renumber clusters to be sequential
    unique_clusters <- sort(unique(cluster_assignments))
    cluster_map <- setNames(1:length(unique_clusters), unique_clusters)
    best_kmeans$cluster <- cluster_map[as.character(cluster_assignments)]

    n_clusters <- length(unique_clusters)
    message("  Final number of clusters after merging: ", n_clusters)
  }


  # ========================================================================
  # STEP 3: Create soft label edges based on distance to centroids
  # ========================================================================
  message("  Step 3: Creating soft label edges...")

  # Compute distance from each cell to each centroid
  centroids <- best_kmeans$centers
  cell_to_centroid_dist <- matrix(0, nrow = n_cells, ncol = n_clusters)

  for (k in 1:n_clusters) {
    centroid <- centroids[k, ]
    cell_to_centroid_dist[, k] <- sqrt(rowSums((pca_coords -
                                                  matrix(centroid, nrow = n_cells,
                                                         ncol = length(centroid), byrow = TRUE))^2))
  }

  # Convert distances to similarities using Gaussian kernel
  sigma <- median(cell_to_centroid_dist)
  if (sigma == 0) sigma <- 1

  label_edges <- exp(-cell_to_centroid_dist^2 / (2 * sigma^2))

  # Row-normalize so each cell's edges sum to 1
  row_sums <- rowSums(label_edges)
  row_sums[row_sums == 0] <- 1  # Avoid division by zero
  label_edges <- label_edges / row_sums

  # Filter very weak connections
  label_edges[label_edges < min_weight] <- 0

  # Re-normalize after filtering
  row_sums <- rowSums(label_edges)
  row_sums[row_sums == 0] <- 1
  label_edges <- label_edges / row_sums


  # ========================================================================
  # STEP 4: Balance clusters to prevent dominance
  # ========================================================================
  if (balance_clusters) {
    message("  Step 4: Balancing cluster weights...")

    col_means <- colMeans(label_edges)
    max_mean <- max(col_means)
    max_ratio <- 3  # Don't let any cluster be >3x the average

    for (i in 1:ncol(label_edges)) {
      if (col_means[i] > max_mean / max_ratio) {
        scale_factor <- (max_mean / max_ratio) / col_means[i]
        label_edges[, i] <- label_edges[, i] * scale_factor
        message("    Scaled down cluster ", i, " by ", round(1/scale_factor, 2), "x")
      }
    }

    # Re-normalize rows after balancing
    row_sums <- rowSums(label_edges)
    row_sums[row_sums == 0] <- 1
    label_edges <- label_edges / row_sums
  }


  # ========================================================================
  # STEP 5: Finalize matrix
  # ========================================================================
  rownames(label_edges) <- rownames(mat)
  colnames(label_edges) <- paste0("ATAC_cluster_", 1:n_clusters)

  # Summary statistics
  message("\n  Summary:")
  message("    Output dimensions: ", nrow(label_edges), " cells x ",
          ncol(label_edges), " clusters")
  message("    Mean edge weight per cluster: ",
          paste(round(colMeans(label_edges), 3), collapse = ", "))
  message("    Cells per cluster (hard assignment): ",
          paste(table(apply(label_edges, 1, which.max)), collapse = ", "))

  return(as.matrix(label_edges))
}


# Create motif_regions_df from peaks using motifmatchr
# Input:
#  - peaks_gr: GRanges of ATAC peaks (must match columns in atac_counts_matrix_filtered)
#  - pfm_list: optional PFMatrixList or PWMatrixList (JASPAR, etc). If NULL, will try JASPAR2020 via TFBSTools/JASPAR2020.
#  - min_count: keep motifs present in at least this many peaks
# Output:
#  - data.frame with columns seqnames,start,end,cluster (cluster = motif name)
createMotifRegionsDf <- function(peaks_gr, pfm_list = NULL, min_count = 50, genome = "hg38") {
  if (!is(peaks_gr, "GRanges")) stop("peaks_gr must be a GRanges")
  # try load motifs if not provided
  if (is.null(pfm_list)) {
    if (!requireNamespace("JASPAR2020", quietly = TRUE) || !requireNamespace("TFBSTools", quietly = TRUE) || !requireNamespace("motifmatchr", quietly = TRUE)) {
      stop("Please install JASPAR2020, TFBSTools and motifmatchr to run createMotifRegionsDf without pfm_list.")
    }
    opts <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts=list(collection="CORE", tax_group="vertebrates", matrixtype="PWM"))
    pfm_list <- opts
  }
  if (!requireNamespace("motifmatchr", quietly = TRUE)) stop("Install motifmatchr")
  # match motifs to peaks
  mm <- motifmatchr::matchMotifs(pfm_list, peaks_gr, genome = genome, out = "matches")
  # mm is a MotifMatches object; convert to list of motif->peak indices
  mm_mat <- as.matrix(motifmatchr::motifMatches(mm)) # peaks x motifs (logical)
  motif_counts <- colSums(mm_mat, na.rm = TRUE)
  keep_motifs <- names(motif_counts)[motif_counts >= min_count]
  if (length(keep_motifs) == 0) stop("No motifs passed min_count threshold; reduce min_count or provide pfm_list.")
  # for each motif, get peaks where it occurs and convert those peaks to region rows
  motif_regions <- lapply(keep_motifs, function(m) {
    peak_idx <- which(mm_mat[, m])
    gr_sub <- peaks_gr[peak_idx]
    data.frame(seqnames = as.character(seqnames(gr_sub)), start = start(gr_sub), end = end(gr_sub), cluster = m, stringsAsFactors = FALSE)
  })
  motif_regions_df <- do.call(rbind, motif_regions)
  rownames(motif_regions_df) <- NULL
  return(motif_regions_df)
}




# Helper function for generic spatial score plotting
plot_score_spatial <- function(df, score_col, title, ...) {
  if (nrow(df) == 0 || !score_col %in% colnames(df) || all(is.na(df[[score_col]]))) {
    message("No valid scores to plot spatially.")
    return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No score data") + ggplot2::theme_void())
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = .data[[score_col]])) +
    ggplot2::coord_fixed() + ggplot2::ggtitle(title) + ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(na.value = "grey50") +
    ggplot2::xlab("Spatial X") + ggplot2::ylab("Spatial Y")
  return(p)
}

# Helper function for score distribution plotting
plot_score_distribution <- function(df, score_col, title, score_method_label) {
  if (nrow(df) == 0 || !score_col %in% colnames(df) || all(is.na(df[[score_col]]))) {
    message("No valid scores for distribution plot.")
    return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No score data") + ggplot2::theme_void())
  }

  # Filter out NA scores for density calculation
  valid_scores_df <- df[!is.na(df[[score_col]]), c(score_col, "cell_id"), drop = FALSE]

  if (nrow(valid_scores_df) > 0) {
    p <- ggplot2::ggplot(valid_scores_df, ggplot2::aes(x = .data[[score_col]])) +
      ggplot2::geom_density(fill = "lightgreen", alpha = 0.7, na.rm = TRUE) +
      ggplot2::geom_vline(xintercept = median(valid_scores_df[[score_col]], na.rm = TRUE), color = "darkgreen", linetype = "dashed") +
      ggplot2::labs(title = title,
                    x = paste0("Score (", score_method_label, ")"), y = "Density") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    return(p)
  } else {
    message("No valid scores to plot distribution.")
    return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No score data") + ggplot2::theme_void())
  }
}


#' Spatial Contextual Feature Similarity Score
#'
#' For each cell, compute:
#'  - S_spatial: median cosine similarity (in a feature space) to its k spatial nearest neighbours
#'  - S_random: median cosine similarity to k random cells
#'  - score: S_spatial - S_random (or a normalized ratio)
#'
#' This uses the same cosine machinery as assessSpatialGraphCoherenceOverall().
#'
#' @param feature_matrix_for_similarity sparse matrix (features x cells). Colnames are cell IDs.
#' @param spatial_coords_df data.frame with columns x,y; rownames are cell IDs.
#' @param k_neighbors integer, number of spatial neighbours (and random cells) to use.
#' @param modality_name character, just for messages/plots.
#' @param n_workers integer, for parallel; NULL -> detect.
#' @param batch_size integer, batch size for parallel.
#' @param n_pca_embed integer or NA. If >0, will PCA-embed features before cosine (like assessSpatialGraphCoherenceOverall).
#' @param score_method "difference" or "normalized_ratio".
#'
#' @return list with:
#'   - scores_df: data.frame with columns cell_id, x, y, S_spatial, S_random, score
#'   - spatial_plot: ggplot (score in spatial coordinates)
#'   - score_distribution_plot: ggplot (score density)
#' @export
computeSpatialContextualFeatureScore <- function(
    feature_matrix_for_similarity,
    spatial_coords_df,
    k_neighbors = 10,
    modality_name = "Modality",
    n_workers = NULL,
    batch_size = 1000L,
    n_pca_embed = 50L,
    score_method = c("difference", "normalized_ratio")
) {
  use_irlba <- requireNamespace("irlba", quietly = TRUE)

  if (!is(feature_matrix_for_similarity, "Matrix"))
    stop("feature_matrix_for_similarity must be a sparse Matrix (features x cells).")
  if (!is.data.frame(spatial_coords_df) || ncol(spatial_coords_df) < 2)
    stop("spatial_coords_df must be data.frame with at least columns x,y.")
  if (is.null(colnames(feature_matrix_for_similarity)) || is.null(rownames(spatial_coords_df)))
    stop("feature_matrix_for_similarity colnames and spatial_coords_df rownames must be cell IDs.")

  # Align cells
  feature_cell_ids <- colnames(feature_matrix_for_similarity)
  spatial_cell_ids <- rownames(spatial_coords_df)
  common_cells <- intersect(feature_cell_ids, spatial_cell_ids)
  if (length(common_cells) == 0) stop("No common cells between feature matrix and spatial coords.")
  if (length(common_cells) < length(feature_cell_ids) || length(common_cells) < nrow(spatial_coords_df)) {
    message("Subsetting to ", length(common_cells), " common cells.")
  }

  feature_matrix_for_similarity <- feature_matrix_for_similarity[, common_cells, drop = FALSE]
  spatial_coords_df <- spatial_coords_df[common_cells, , drop = FALSE]
  num_cells <- length(common_cells)
  cell_ids <- common_cells

  if (num_cells <= 1) {
    empty_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data") +
      ggplot2::theme_void()
    return(list(
      scores_df = data.frame(),
      spatial_plot = empty_plot,
      score_distribution_plot = empty_plot
    ))
  }

  if (k_neighbors >= num_cells) {
    warning("k_neighbors >= num_cells, setting to num_cells-1")
    k_neighbors <- max(1, num_cells - 1)
  }

  # Parallel plan
  if (is.null(n_workers)) {
    n_workers <- max(1L, parallel::detectCores() - 1L)
  } else {
    n_workers <- max(1L, as.integer(n_workers))
  }
  old_plan <- future::plan()
  options(future.globals.maxSize = +Inf)
  future::plan(future::multisession, workers = n_workers)
  on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)

  message("Computing spatial contextual feature scores for ", modality_name,
          " (n_cells=", num_cells, ", k=", k_neighbors, ")")

  # --- Prepare spatial matrix (normalized) ---
  spatial_coords_mat_raw <- as.matrix(spatial_coords_df[cell_ids, c("x", "y")])
  x_range <- range(spatial_coords_mat_raw[, 1])
  y_range <- range(spatial_coords_mat_raw[, 2])
  max_dim_range <- max(diff(x_range), diff(y_range))
  if (max_dim_range == 0) {
    spatial_coords_mat_normalized <- spatial_coords_mat_raw
  } else {
    spatial_coords_mat_normalized <- sweep(spatial_coords_mat_raw, 2,
                                           c(min(x_range), min(y_range)), "-")
    spatial_coords_mat_normalized <- spatial_coords_mat_normalized / max_dim_range
  }

  # --- Feature matrix + optional PCA embedding (same pattern as assessSpatialGraphCoherenceOverall) ---
  feature_matrix_for_similarity_to_use <- feature_matrix_for_similarity
  feature_matrix_norms <- sqrt(Matrix::colSums(feature_matrix_for_similarity_to_use *
                                                 feature_matrix_for_similarity_to_use))

  if (!is.null(n_pca_embed) && !is.na(n_pca_embed) && n_pca_embed > 0 && use_irlba) {
    message("PCA embedding (n_pca_embed=", n_pca_embed, ") for faster cosine...")
    dense_mat <- as.matrix(feature_matrix_for_similarity)  # F x C
    sv <- irlba::irlba(dense_mat, nv = n_pca_embed)
    feature_matrix_embed <- t(sv$v %*% diag(sv$d))  # k x C

    if (ncol(feature_matrix_embed) == length(cell_ids)) {
      colnames(feature_matrix_embed) <- cell_ids
      rownames(feature_matrix_embed) <- paste0("PC", seq_len(n_pca_embed))
      feature_matrix_for_similarity_to_use <- feature_matrix_embed
      feature_matrix_norms <- sqrt(Matrix::colSums(feature_matrix_for_similarity_to_use *
                                                     feature_matrix_for_similarity_to_use))
    } else {
      warning("PCA embedding dimension mismatch; using original feature matrix.")
    }
    rm(dense_mat, sv); gc()
  }

  # --- Spatial nearest neighbours (like Metric 2 in assessSpatialGraphCoherenceOverall) ---
  message("Finding spatial nearest neighbours (RANN::nn2)...")
  all_spatial_nn_results <- RANN::nn2(data = spatial_coords_mat_normalized,
                                      k = k_neighbors + 1)
  all_nn_idx <- all_spatial_nn_results$nn.idx

  # --- Batch processing: compute S_spatial and S_random using cosine similarity ---
  process_batch <- function(sources) {
    nb <- length(sources)
    S_spatial_vec <- rep(NA_real_, nb)
    S_random_vec  <- rep(NA_real_, nb)

    all_indices <- seq_len(num_cells)

    for (ii in seq_along(sources)) {
      i <- sources[ii]

      # spatial neighbors: drop self; take top k
      raw_targets_full <- all_nn_idx[i, ]
      raw_targets_full <- raw_targets_full[raw_targets_full != i & raw_targets_full > 0]
      if (length(raw_targets_full) == 0) {
        S_spatial_vec[ii] <- NA_real_
      } else {
        spatial_targets <- raw_targets_full[1:min(length(raw_targets_full), k_neighbors)]
        neighbor_cell_names <- cell_ids[spatial_targets]

        # map to columns in feature matrix
        matrix_colnames <- colnames(feature_matrix_for_similarity_to_use)
        neighbor_col_idx <- match(neighbor_cell_names, matrix_colnames)
        neighbor_col_idx <- neighbor_col_idx[!is.na(neighbor_col_idx)]
        if (length(neighbor_col_idx) == 0) {
          S_spatial_vec[ii] <- NA_real_
        } else {
          mat2 <- feature_matrix_for_similarity_to_use[, neighbor_col_idx, drop = FALSE]
          vec1 <- feature_matrix_for_similarity_to_use[, i, drop = FALSE]

          dot_products <- as.vector(crossprod(vec1, mat2))
          norm1 <- feature_matrix_norms[i]
          norm2 <- feature_matrix_norms[neighbor_col_idx]
          denom <- norm1 * norm2
          cos_sim <- ifelse(denom == 0, 0, dot_products / denom)
          cos_sim <- pmax(0, pmin(1, cos_sim))

          S_spatial_vec[ii] <- median(cos_sim, na.rm = TRUE)
        }
      }

      # random neighbors: sample k different cells (no self)
      pool <- all_indices[all_indices != i]
      if (length(pool) == 0) {
        S_random_vec[ii] <- NA_real_
      } else {
        rand_targets <- sample(pool, min(k_neighbors, length(pool)))
        rand_cell_names <- cell_ids[rand_targets]
        matrix_colnames <- colnames(feature_matrix_for_similarity_to_use)
        rand_col_idx <- match(rand_cell_names, matrix_colnames)
        rand_col_idx <- rand_col_idx[!is.na(rand_col_idx)]
        if (length(rand_col_idx) == 0) {
          S_random_vec[ii] <- NA_real_
        } else {
          mat_rand <- feature_matrix_for_similarity_to_use[, rand_col_idx, drop = FALSE]
          vec1 <- feature_matrix_for_similarity_to_use[, i, drop = FALSE]

          dot_products_rand <- as.vector(crossprod(vec1, mat_rand))
          norm1 <- feature_matrix_norms[i]
          norm2_rand <- feature_matrix_norms[rand_col_idx]
          denom_rand <- norm1 * norm2_rand
          cos_sim_rand <- ifelse(denom_rand == 0, 0, dot_products_rand / denom_rand)
          cos_sim_rand <- pmax(0, pmin(1, cos_sim_rand))

          S_random_vec[ii] <- median(cos_sim_rand, na.rm = TRUE)
        }
      }
    }

    list(idx = sources,
         S_spatial = S_spatial_vec,
         S_random  = S_random_vec)
  }

  seq_starts <- seq(1L, num_cells, by = batch_size)
  batch_idx_list <- lapply(seq_starts, function(s) {
    e <- min(s + batch_size - 1L, num_cells)
    seq.int(s, e)
  })

  message("Computing spatial vs random feature similarity in batches...")
  batch_results <- future.apply::future_lapply(batch_idx_list, FUN = process_batch)

  S_spatial_all <- rep(NA_real_, num_cells)
  S_random_all  <- rep(NA_real_, num_cells)
  for (br in batch_results) {
    S_spatial_all[br$idx] <- br$S_spatial
    S_random_all[br$idx]  <- br$S_random
  }

  # --- Compute score ---
  score_all <- rep(NA_real_, num_cells)
  if (score_method == "difference") {
    # treat missing as 0 in difference
    score_all <- (ifelse(is.na(S_spatial_all), 0, S_spatial_all) -
                    ifelse(is.na(S_random_all),  0, S_random_all))
  } else if (score_method == "normalized_ratio") {
    epsilon <- 1e-6
    for (i in seq_len(num_cells)) {
      s <- S_spatial_all[i]
      r <- S_random_all[i]
      if (!is.na(s) && !is.na(r)) {
        denom <- s + r
        if (abs(denom) < epsilon) {
          score_all[i] <- s - r
        } else {
          score_all[i] <- (s - r) / (denom + epsilon)
        }
      } else if (!is.na(s) && is.na(r)) {
        score_all[i] <- s
      } else if (is.na(s) && !is.na(r)) {
        score_all[i] <- -r
      } else {
        score_all[i] <- NA_real_
      }
    }
  }

  scores_df <- data.frame(
    cell_id  = cell_ids,
    x        = spatial_coords_df[cell_ids, "x"],
    y        = spatial_coords_df[cell_ids, "y"],
    S_spatial = S_spatial_all,
    S_random  = S_random_all,
    score     = score_all,
    stringsAsFactors = FALSE
  )

  # --- Plots using your helpers (already defined in your script) ---
  spatial_plot <- plot_score_spatial(
    df = scores_df,
    score_col = "score",
    title = paste0(modality_name, ": Contextual Feature Similarity Score")
  )

  score_distribution_plot <- plot_score_distribution(
    df = scores_df,
    score_col = "score",
    title = paste0(modality_name, ": Contextual Feature Similarity Score Distribution"),
    score_method_label = score_method
  )

  list(
    scores_df = scores_df,
    spatial_plot = spatial_plot,
    score_distribution_plot = score_distribution_plot
  )
}



# Align labels to a joint universe by zero-padding missing cells
# This ensures that RNA labels exist for RNA cells and 0-placeholders exist for ATAC cells
alignLabelEdgesToJointUniverse <- function(labelEdges, all_graph_ids) {
  if (is.null(labelEdges)) return(NULL)

  # Create empty matrix with all graph cells as rows
  full_le <- matrix(0,
                    nrow = length(all_graph_ids),
                    ncol = ncol(labelEdges),
                    dimnames = list(all_graph_ids, colnames(labelEdges)))

  # Fill in the labels for the cells that actually have them
  common <- intersect(rownames(labelEdges), all_graph_ids)
  if (length(common) == 0) {
    warning("No overlap between labelEdges and graph IDs.")
    return(full_le)
  }

  full_le[common, ] <- as.matrix(labelEdges[common, , drop = FALSE])
  return(full_le)
}



# Subset logic for joint graphs
subset_joint_cells_for_mapping <- function(cg, rna_le, atac_le, max_cells = 5000L) {
  all_ids <- rownames(cg)

  # Determine which IDs belong to which modality based on the input label edges
  # (Assuming labelEdges rownames are the ground truth for cell identity)
  rna_pool <- intersect(all_ids, rownames(rna_le))
  atac_pool <- if(!is.null(atac_le)) intersect(all_ids, rownames(atac_le)) else c()

  if (length(all_ids) > max_cells) {
    message("  Subsetting joint graph to ", max_cells, " cells...")
    set.seed(123)
    # We try to keep a balanced representation of both modalities if possible
    keep_ids <- sample(all_ids, max_cells)
  } else {
    keep_ids <- all_ids
  }

  # Subset the graph
  cg_sub <- cg[keep_ids, keep_ids, drop = FALSE]

  # Align/Pad the labels to the subsetted IDs
  rna_le_sub  <- alignLabelEdgesToJointUniverse(rna_le, keep_ids)
  atac_le_sub <- alignLabelEdgesToJointUniverse(atac_le, keep_ids)

  # Create modality grouping: 1 for RNA-derived cells, 2 for ATAC-derived cells
  # If a cell is in both (rare), it defaults to 1
  modality_groups <- ifelse(keep_ids %in% rownames(rna_le), 1, 2)

  list(
    graph = cg_sub,
    rna_le = rna_le_sub,
    atac_le = atac_le_sub,
    cell_ids = keep_ids,
    groups = modality_groups
  )
}


# Align labels to the Joint Graph by zero-padding the "other" modality
alignLabelEdgesToJointGraph <- function(labelEdges, all_graph_ids) {
  if (is.null(labelEdges)) return(NULL)

  # Create empty matrix (All Cells x Labels)
  full_le <- matrix(0,
                    nrow = length(all_graph_ids),
                    ncol = ncol(labelEdges),
                    dimnames = list(all_graph_ids, colnames(labelEdges)))

  # Only fill the rows for which we actually have data
  common <- intersect(rownames(labelEdges), all_graph_ids)
  full_le[common, ] <- as.matrix(labelEdges[common, , drop = FALSE])
  return(full_le)
}

# Joint subsetting logic
subset_joint_cells_for_mapping <- function(cg, rna_le, atac_le, max_cells = 5000L) {
  all_ids <- rownames(cg)

  # Weighted sampling to ensure we get a mix of RNA and ATAC cells
  rna_pool <- intersect(all_ids, rownames(rna_le))
  atac_pool <- intersect(all_ids, rownames(atac_le))

  set.seed(123)
  if (length(all_ids) > max_cells) {
    # Sample proportionally from both pools
    keep_rna <- sample(rna_pool, min(length(rna_pool), max_cells/2))
    keep_atac <- sample(atac_pool, min(length(atac_pool), max_cells/2))
    keep_ids <- c(keep_rna, keep_atac)
  } else {
    keep_ids <- all_ids
  }

  cg_sub <- cg[keep_ids, keep_ids]

  # Align/Pad labels to this specific subset
  rna_le_aligned  <- alignLabelEdgesToJointGraph(rna_le, keep_ids)
  atac_le_aligned <- alignLabelEdgesToJointGraph(atac_le, keep_ids)

  # Define modality groups: RNA=1, ATAC=2
  groups <- ifelse(keep_ids %in% rna_pool, 1, 2)

  list(graph = cg_sub, rna_le = rna_le_aligned, atac_le = atac_le_aligned,
       cell_ids = keep_ids, groups = groups)
}
