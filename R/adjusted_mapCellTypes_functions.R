# compute z score
compute_zscore <- function(info, info_rand, nround)
{
  info_mean = matrix(0, nrow(info), ncol(info))
  info_var = matrix(0, nrow(info), ncol(info))
  for(i in 1:nround)
  {
    info_mean = info_mean + info_rand[[i]]
    info_var = info_var + info_rand[[i]]^2
  }

  info_mean = info_mean/nround
  info_var = info_var/nround
  info_std = info_var - info_mean^2
  if(any(info_std < 0)) {
    message('some of variances are zero or negative when computing Z-score')
    info_std[info_std < 0] = 0
  }
  zscore = (info - info_mean)/sqrt(info_std)
  zscore[zscore <0] = 0
  return(list('zscore' = zscore, 'mean' = info_mean, 'var' = info_var))
}


compute_zscore_fast <- function(info, info_rand, nround) {
  # Handle edge cases
  if(nround < 2) {
    warning("Need at least 2 permutation rounds for Z-score calculation")
    return(list('zscore' = matrix(0, nrow(info), ncol(info)),
                'mean' = info,
                'var' = matrix(0, nrow(info), ncol(info))))
  }

  # Ensure all matrices have same dimensions
  dims_check <- sapply(info_rand, function(x) identical(dim(x), dim(info)))
  if(!all(dims_check)) {
    stop("All permutation matrices must have same dimensions as info")
  }

  # Stack into 3D array safely
  info_array <- array(0, dim = c(nrow(info), ncol(info), nround))
  for(i in 1:nround) {
    info_array[,,i] <- as.matrix(info_rand[[i]])
  }

  # Vectorized mean and variance with dimension checks
  if(nround == 1) {
    info_mean <- info_array[,,1]
    info_var <- matrix(0, nrow(info), ncol(info))
  } else {
    info_mean <- apply(info_array, c(1, 2), mean)
    info_var <- apply(info_array, c(1, 2), var)
  }

  # Handle negative/zero variance
  info_var[info_var < 1e-10] <- 1e-10
  info_var[is.na(info_var)] <- 1e-10

  # Z-score calculation
  zscore <- (info - info_mean) / sqrt(info_var)
  zscore[zscore < 0] <- 0
  zscore[is.na(zscore)] <- 0

  return(list('zscore' = zscore, 'mean' = info_mean, 'var' = info_var))
}


simulate_rand0_ <- function(labelEdges2) {
  labelEdges2 <- as.matrix(labelEdges2)
  if(nrow(labelEdges2) == 0) return(labelEdges2)

  # 1. Calculate cell margins (handling dimension drops)
  cell_margin = apply(labelEdges2, 1, function(x){
    x[x > 1] = 1
    l = sapply(seq(0, 1, by = 0.2), function(y) mean(x <= y))
    c(l[1], rep(diff(l), each = 2)/2)
  })
  if (is.vector(cell_margin)) cell_margin <- matrix(cell_margin, ncol = 1)

  # 2. Calculate label margins
  label_margin = apply(labelEdges2, 2, function(x){
    x[x > 1] = 1
    l = sapply(seq(0, 1, by = 0.1), function(y) sum(x <= y))
    c(l[1], diff(l))
  })
  if (is.vector(label_margin)) label_margin <- matrix(label_margin, ncol = ncol(labelEdges2))

  labelEdges2_rand = matrix(0, nrow(labelEdges2), ncol(labelEdges2))
  cuts = seq(0, 1, by = 0.1)

  for(j in 1:ncol(labelEdges2)) {
    probs = label_margin[, j]
    idx = 1:nrow(labelEdges2)

    for(k in length(probs):2) {
      if (probs[k] == 0) next

      if(length(idx) <= probs[k]) {
        aa = idx
      } else {
        # FIX: Use sample.int to prevent R from expanding single integers into 1:x
        # and ensure the probability vector matches length(idx)
        prob_vec <- cell_margin[k, idx] + 0.01
        aa_idx = sample.int(length(idx), probs[k], prob = prob_vec)
        aa = idx[aa_idx]
      }

      if(k > 1 && length(aa) > 0) {
        labelEdges2_rand[aa, j] = runif(length(aa), cuts[k-1], cuts[k])
      }
      idx = setdiff(idx, aa)
    }
  }

  colnames(labelEdges2_rand) = colnames(labelEdges2)
  rownames(labelEdges2_rand) = rownames(labelEdges2)
  return(labelEdges2_rand)
}


simulate_rand0_fast <- function(labelEdges2) {
  labelEdges2 <- as.matrix(labelEdges2)
  if(nrow(labelEdges2) == 0) return(labelEdges2)

  # Handle single-cell case
  if(nrow(labelEdges2) == 1) {
    return(labelEdges2)  # Can't permute a single cell
  }

  # Handle single-label case
  if(ncol(labelEdges2) == 1) {
    # Simple shuffle for single column
    labelEdges2_rand <- matrix(sample(labelEdges2[,1]), ncol = 1)
    colnames(labelEdges2_rand) <- colnames(labelEdges2)
    rownames(labelEdges2_rand) <- rownames(labelEdges2)
    return(labelEdges2_rand)
  }

  # Precompute cuts once
  cuts <- seq(0, 1, by = 0.1)
  n_bins <- length(cuts) - 1

  # Discretize all values at once (vectorized)
  bin_indices <- findInterval(pmin(labelEdges2, 1), cuts, rightmost.closed = TRUE)

  # Calculate cell margins with dimension preservation
  cell_margin <- matrix(0, nrow = n_bins + 1, ncol = nrow(labelEdges2))
  for(i in 1:(n_bins + 1)) {
    if(is.matrix(bin_indices)) {
      cell_margin[i, ] <- colMeans(bin_indices <= i)
    } else {
      cell_margin[i, ] <- mean(bin_indices <= i)
    }
  }

  # Safe diff calculation
  if(nrow(cell_margin) > 1) {
    cell_margin_diff <- apply(cell_margin, 2, diff)
    if(!is.matrix(cell_margin_diff)) {
      cell_margin_diff <- matrix(cell_margin_diff, nrow = n_bins)
    }
    cell_margin <- rbind(cell_margin[1, , drop = FALSE],
                         cell_margin_diff / 2)
  }

  # Calculate label margins with dimension preservation
  label_margin <- matrix(0, nrow = n_bins + 1, ncol = ncol(labelEdges2))
  for(i in 1:(n_bins + 1)) {
    if(is.matrix(bin_indices)) {
      label_margin[i, ] <- colSums(bin_indices <= i)
    } else {
      label_margin[i, ] <- sum(bin_indices <= i)
    }
  }

  # Safe diff calculation
  if(nrow(label_margin) > 1) {
    label_margin_diff <- apply(label_margin, 2, diff)
    if(!is.matrix(label_margin_diff)) {
      label_margin_diff <- matrix(label_margin_diff, nrow = n_bins, ncol = ncol(labelEdges2))
    }
    label_margin <- rbind(label_margin[1, , drop = FALSE], label_margin_diff)
  }

  # Vectorized assignment
  labelEdges2_rand <- matrix(0, nrow(labelEdges2), ncol(labelEdges2))

  for(j in 1:ncol(labelEdges2)) {
    probs <- label_margin[, j]
    idx <- 1:nrow(labelEdges2)

    for(k in length(probs):2) {
      if (probs[k] == 0 || length(idx) == 0) next

      n_sample <- min(probs[k], length(idx))

      # Safe probability vector extraction
      if(is.matrix(cell_margin)) {
        prob_vec <- pmax(cell_margin[k, idx], 0.01)
      } else {
        prob_vec <- rep(0.01, length(idx))
      }

      # Use sample.int with explicit length check
      if(length(idx) == 1) {
        aa <- idx
      } else if(n_sample >= length(idx)) {
        aa <- idx
      } else {
        aa_idx <- sample.int(length(idx), n_sample, prob = prob_vec)
        aa <- idx[aa_idx]
      }

      if(k > 1 && length(aa) > 0) {
        labelEdges2_rand[aa, j] <- runif(length(aa), cuts[k-1], cuts[k])
      }
      idx <- setdiff(idx, aa)
    }
  }

  colnames(labelEdges2_rand) <- colnames(labelEdges2)
  rownames(labelEdges2_rand) <- rownames(labelEdges2)
  return(labelEdges2_rand)
}


# check label Edges and make labelEdges having the same rows as cellGraph
checkLabelEdges_ <- function(labelEdges, groups, cellGraph, tips = NULL, suffix = NULL)
{
  if(missing(labelEdges) || (!is(labelEdges, "data.frame") & (!is(labelEdges, "matrix") & (!is(labelEdges, "Matrix"))))){
    stop("Must provide a dataframe or matrix of cell-to-label edges")
  }
  if(is.null(colnames(labelEdges)) || is.null(rownames(labelEdges))){
    stop("labelEdges must have cell barcodes as rownames and cell type labels as colnames")
  }
  if(!is(groups, 'numeric') || length(groups) != nrow(labelEdges))
  {
    stop("groups must be numeric and length of groups must be the same as the number of rows of labelEdges")
  }
  #groups = as.numeric(as.factor(groups)) # change groups to numeric
  names(groups) = rownames(labelEdges)
  stopifnot("groups contain NA" = all(!is.na(groups)))

  if(nrow(labelEdges) == nrow(cellGraph) && all(rownames(labelEdges) == rownames(cellGraph))) return(list(labelEdges, groups))
  rowIdx = match(rownames(cellGraph), rownames(labelEdges))
  if(all(is.na(rowIdx))) warning('no cell barcodes are common in labelEdges and cellGraph')
  groups = groups[rowIdx]
  labelEdges = labelEdges[rowIdx, ]
  groups[is.na(groups)] = 0
  labelEdges[is.na(labelEdges)] = 0
  rownames(labelEdges) = names(groups) = rownames(cellGraph)
  if(!is.null(tips)) {
    if(!all(tips %in% colnames(labelEdges)))
    {
      stop('missing cell-to-label edges for tree tips: ', tips)
    }
    labelEdges = labelEdges[, tips]
  }
  if(!is.null(suffix)) colnames(labelEdges) = paste(colnames(labelEdges), suffix, sep = '_')
  return(list(labelEdges,groups))
}


# Convergence checker for early stopping
check_convergence <- function(info_rand_list, window = 20, threshold = 0.01) {
  n <- length(info_rand_list)
  if(n < window + 10) return(FALSE)

  # Verify all matrices exist and have same dimensions
  if(any(sapply(info_rand_list, is.null))) return(FALSE)

  # Compare recent window to previous window
  recent_mats <- info_rand_list[(n - window + 1):n]
  prev_mats <- info_rand_list[(n - 2*window + 1):(n - window)]

  # Average matrices in each window (safe element-wise)
  recent_mean <- Reduce('+', recent_mats) / window
  prev_mean <- Reduce('+', prev_mats) / window

  # Check mean absolute difference
  diff <- mean(abs(recent_mean - prev_mean), na.rm = TRUE)

  if(is.na(diff) || is.infinite(diff)) return(FALSE)

  message(sprintf("  Convergence check at round %d: diff = %.6f (threshold = %.6f)",
                  n, diff, threshold))

  return(diff < threshold)
}


mapCellTypes_ <- function(cellGraph, labelEdgesList, labelEdgeWeights = NULL, wtrees = NULL, treeList = NULL, compute.Zscore = TRUE,  nround = 50,
                         groupsList = NULL, sampleDepth =2000, batch_size = 20, conv_threshold = 0.01, conv_window = 20,...) # groups, permutation within cell groups, cells same order as cellEdges
{

  args = list(...)
  if('steps' %in% names(args)){
    steps = args[['steps']]
  }else{
    steps = Inf
  }
  if('tensorflow' %in% names(args)){
    tensorflow = args[['tensorflow']]
  }else{
    tensorflow = F
  }

  if(missing(cellGraph) || (!is(cellGraph, "data.frame") & (!is(cellGraph, "matrix") & (!is(cellGraph, "Matrix"))))){
    stop("Must provide a dataframe or matrix of cell-to-cell similarity graph")
  }
  if(is.null(colnames(cellGraph)) || is.null(rownames(cellGraph))){
    stop("cellGraph must have same rownames and colnames as cell barcodes")
  }
  if(any(colnames(cellGraph) != rownames(cellGraph)))
  {
    stop('colname and rowname of cellGraph must be the same')
  }
  diag(cellGraph) = 0

  if(missing(labelEdgesList) || !is(labelEdgesList, "list")){
    stop("Must provide a list of cell-to-label edges matrices")
  }
  if(!is.null(groupsList) & (!is(groupsList, "list") || length(groupsList) != length(labelEdgesList))){
    stop('groupsList must be a list of vectors and must have the same length as labelEdgesList')
  }
  if(is.null(groupsList))
  {
    groupsList = lapply(labelEdgesList, function(labelEdges) rep(1, nrow(labelEdges)))
  }
  if(!is.null(treeList) & length(treeList)!=length(labelEdgesList))
  {
    stop('treeList must have the same length as labelEdgesList')
  }

  for (i in seq_along(labelEdgesList)){
    labelEdges = labelEdgesList[[i]]
    groups = groupsList[[i]]
    if(!is.null(treeList))
    {
      stopifnot('tree must be provided as a phylo object'=is(treeList[[i]], 'phylo'))
      stopifnot('tree must have tip labels' = !is.null(treeList[[i]]$tip.label))
      res = checkLabelEdges_(labelEdges, groups, cellGraph, tips = treeList[[i]]$tip.label, suffix = i) #reorder columns of labelEdges as the cell type matrix
    }else{
      res = checkLabelEdges_(labelEdges, groups, cellGraph, suffix = i) # adding '_x' to cell type labels to avoid duplicate
    }

    labelEdgesList[[i]] = res[[1]]
    groupsList[[i]] = res[[2]]
  }

  if(!is.null(treeList))
  {
    if(!is.null(wtrees))
    {
      if(!(is(wtrees, "data.frame") | is(wtrees, "matrix")) | nrow(wtrees) != length(labelEdgesList) | ncol(wtrees) !=2 )
      {
        stop('tree edge weights must be a dataframe/matrix with two columns and same number of rows as labelEdgesList')
      }
    }else{
      wtrees = matrix(1, length(labelEdgesList), 2)
    }

    res = lapply(seq_along(treeList), function(i){
      tr = treeList[[i]]
      w = wtrees[i, ]
      tree2Mat(tr,w[1], w[2], i) # adding '_x' to tips to avoid duplicates
    })
    allCellTypes = unlist(sapply(res, function(x) x[[2]]))
    ncellTypes = sapply(res, function(x) length(x[[2]]))
    cellTypesM = Matrix::bdiag(sapply(res, function(x) x[1])) # becomes a sparse matrix
    ##cellTypesM = as.matrix(cellTypesM)
    colnames(cellTypesM) = rownames(cellTypesM) = allCellTypes

  }else{
    allCellTypes = unlist(sapply(labelEdgesList, colnames))
    ncellTypes = sapply(labelEdgesList, ncol)
    cellTypesM = matrix(0, length(allCellTypes), length(allCellTypes))
    colnames(cellTypesM) = rownames(cellTypesM) = allCellTypes
  }
  l_all = length(allCellTypes)

  if(is.null(labelEdgeWeights))
  {
    message('tunning labelEdgeWeights... ')
    edgeWeights <- tuneEdgeWeights(cellGraph,
                                   labelEdgesList,
                                   sampleDepth = sampleDepth,
                                   ...)
    opt_idx = which.max(edgeWeights$cellHomogeneity)
    message('cellHomogeneity at optimal edgeWeight:')
    print(edgeWeights[opt_idx,])
    labelEdgeWeights = unlist(edgeWeights[opt_idx, ])
  }else if(!is(labelEdgeWeights, 'numeric') || length(labelEdgeWeights)!=length(labelEdgesList)){
    stop("Must provide a numeric value weight for each set of cell-to-label edges")
  }

  message('Pre-computing static graph components...')

  # These don't change across permutations
  static_cell_graph <- cellGraph
  static_cellTypesM <- cellTypesM

  # Pre-allocate combined graph structure
  n_labels <- nrow(cellTypesM)
  n_cells <- nrow(cellGraph)
  total_size <- n_labels + n_cells

  # Create index map for fast assembly
  idx_labels <- 1:n_labels
  idx_cells <- (n_labels + 1):total_size

  # construct the whole graph permuted or original
  if(!compute.Zscore) nround = 0
  message('run CellWalker:')
  library(foreach)
  # construct the whole graph permuted or original
  if(!compute.Zscore) nround = 0
  message('run CellWalker with adaptive convergence:')

  # Initialize
  info1_rand <- list()
  batch_size <- 20  # Run 20 permutations per batch
  max_rounds <- nround
  converged <- FALSE
  rounds_completed <- 0

  # Run permutation round 0 (observed data) first
  message('Computing observed influence (round 0)...')
  expandLabelEdges = lapply(seq_along(labelEdgesList), function(i) {
    labelEdges = labelEdgesList[[i]]
    if(!is.null(treeList)) {
      labelEdgeWeights[i] * cbind(labelEdges, matrix(0, dim(labelEdges)[1], treeList[[i]]$Nnode))
    } else {
      labelEdges * labelEdgeWeights[i]
    }
  })

  cell2label = do.call('cbind', expandLabelEdges)
  combinedGraph = rbind(cbind(cellTypesM, t(cell2label)),
                        cbind(cell2label, cellGraph))
  infMat <- randomWalk(combinedGraph, tensorflow = tensorflow, steps = steps)
  aa = infMat[1:l_all, 1:l_all]
  colnames(aa) = rownames(aa) = allCellTypes
  info1_rand[[1]] <- as.matrix(aa)

  # Validate input dimensions before starting permutations
  message("Validating input data...")
  for(i in seq_along(labelEdgesList)) {
    le <- labelEdgesList[[i]]
    gr <- groupsList[[i]]

    message(sprintf("  Label set %d: %d cells x %d labels", i, nrow(le), ncol(le)))

    # Check for groups with valid data
    for(g in unique(gr[gr != 0])) {
      cells_in_group <- which(gr == g)
      group_data <- le[cells_in_group, , drop = FALSE]
      message(sprintf("    Group %d: %d cells, rowSums range [%.3f, %.3f]",
                      g, length(cells_in_group),
                      min(rowSums(group_data)), max(rowSums(group_data))))

      if(nrow(group_data) == 0 || ncol(group_data) == 0) {
        warning(sprintf("Group %d in label set %d has invalid dimensions", g, i))
      }
    }
  }

  # Now run permutations in batches with convergence checking
  while(rounds_completed < max_rounds && !converged) {
    # Determine batch size for this iteration
    remaining <- max_rounds - rounds_completed
    current_batch <- min(batch_size, remaining)
    batch_start <- rounds_completed + 1
    batch_end <- rounds_completed + current_batch

    message(sprintf('Running permutation rounds %d-%d...', batch_start, batch_end))

    # Run batch in parallel
    batch_results <- foreach(r = batch_start:batch_end,
                             .packages = c('Matrix'),
                             .errorhandling = 'pass') %dopar% {
                               tryCatch({
                                 ## resample labelEdges
                                 expandLabelEdges = lapply(seq_along(labelEdgesList), function(i) {
                                   labelEdges_rand = labelEdges = labelEdgesList[[i]]
                                   groups = groupsList[[i]]
                                   if(any(groups != 0)) {
                                     for(g in unique(groups)) {
                                       if(g == 0) next
                                       cells = which(groups == g)

                                       # Extra safety: skip if no cells or invalid subset
                                       if(length(cells) == 0) next
                                       cell_subset <- labelEdges[cells, , drop = FALSE]
                                       if(nrow(cell_subset) == 0 || ncol(cell_subset) == 0) next

                                       labelEdges_rand[cells, ] = simulate_rand0_fast(cell_subset)
                                     }
                                   }
                                   if(!is.null(treeList)) {
                                     labelEdgeWeights[i] * cbind(labelEdges_rand, matrix(0, dim(labelEdges)[1], treeList[[i]]$Nnode))
                                   } else {
                                     labelEdges_rand * labelEdgeWeights[i]
                                   }
                                 })

                                 cell2label = do.call('cbind', expandLabelEdges)
                                 combinedGraph = rbind(cbind(cellTypesM, t(cell2label)),
                                                       cbind(cell2label, cellGraph))

                                 infMat <- randomWalk(combinedGraph, tensorflow = tensorflow, steps = steps)
                                 aa = infMat[1:l_all, 1:l_all]
                                 colnames(aa) = rownames(aa) = allCellTypes
                                 return(as.matrix(aa))
                               }, error = function(e) {
                                 message("Error in permutation round ", r, ": ", e$message)
                                 return(NULL)  # Return NULL on error
                               })
                             }

    # Filter out NULL results (failed rounds)
    batch_results <- batch_results[!sapply(batch_results, is.null)]

    # Check if we got any valid results
    if(length(batch_results) == 0) {
      stop("All permutation rounds in batch failed. Check data dimensions and structure.")
    }

    # Append batch results
    info1_rand <- c(info1_rand, batch_results)
    rounds_completed <- batch_end

    # Check convergence (only after sufficient rounds)
    if(rounds_completed >= 40) {
      converged <- check_convergence(info1_rand[-1], window = conv_window, threshold = conv_threshold)
      if(converged) {
        message(sprintf('*** Converged after %d rounds (requested %d) ***',
                        rounds_completed, max_rounds))
      }
    }
  }

  message(sprintf('Completed %d permutation rounds', rounds_completed))

  info1 = info1_rand[[1]]
  actual_nround <- length(info1_rand) - 1  # Exclude round 0

  if(compute.Zscore && actual_nround > 0){
    # Pass only the permutation results (exclude round 0)
    perm_results <- info1_rand[2:length(info1_rand)]

    # DEBUG: Check dimensions
    message("Debug info:")
    message("  info1 dimensions: ", paste(dim(info1), collapse = " x "))
    message("  Number of permutation results: ", length(perm_results))
    message("  First perm result dimensions: ", paste(dim(perm_results[[1]]), collapse = " x "))

    # Verify we have valid data
    if(length(perm_results) == 0) {
      warning("No permutation results available for Z-score calculation")
      zscore = NULL
      infomean = NULL
      infovar = NULL
    } else {
      reslist = compute_zscore_fast(info1, perm_results, actual_nround)
      zscore = reslist$zscore
      infomean = reslist$mean
      infovar = reslist$var
    }
  } else {
    zscore = NULL
    infomean = NULL
    infovar = NULL
  }
  idx = cumsum(ncellTypes)
  sts = c(1, idx[1:(length(idx)-1)]+1)
  params = expand.grid(1:length(idx), 1:length(idx))
  params = params[params[,1]!=params[,2,], ]
  info = Map(function(u,v) info1[sts[u]:idx[u],sts[v]:idx[v]], params[,2],  params[,1])
  if(compute.Zscore){
    zscore = Map(function(u,v) zscore[sts[u]:idx[u],sts[v]:idx[v]], params[,2],  params[,1])
    infomean = Map(function(u,v) infomean[sts[u]:idx[u],sts[v]:idx[v]], params[,2],  params[,1])
    infovar = Map(function(u,v) infovar[sts[u]:idx[u],sts[v]:idx[v]], params[,2],  params[,1])
  }
  cellWalk = list(infMat=info, zscore= zscore, infomean = infomean, infovar = infovar, labelEdgeWeights = labelEdgeWeights)
  class(cellWalk) = "cellWalk2"
  return(cellWalk)
}
