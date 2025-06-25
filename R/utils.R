#' Combine potential outcomes matrices according to a treatment program
#'
#' @param M0 matrix DxN, untreated potential outcome
#' @param M1 matrix DxN, treated potential outcomes
#' @param Tr vector length N, binary treatment program assigning each
#' element 1,...,N to treated (Tr[i] = 1) or untreated (Tr[i] = 0)
#'
#' @returns matrix DxN, M combines columns of M0 and M1 according to treatment
#' program Tr, such that M[,i] = M0[,i] \* (1 - Tr[i]) + M1[,i] \* Tr[i]
#' @export
combine_mat <- function(M0, M1, Tr) {
  G <- ncol(M1)
  do.call(cbind, lapply(1:G, function(g) {
    if (Tr[g] == 1) {
      return(M1[,g])
    } else {
      return(M0[,g])
    }
  }))
}

#' Column sums of a matrix or returning a 1-d vector
#'
#' @param mat matrix or vector
#'
#' @returns colSums(mat) if matrix, mat if vector
#' @noRd
colSums_wrapper <- function(mat) {
  # allows a vector to be passed in (returns itself)
  # otherwise row means
  if (!('matrix' %in% class(mat))) {
    return(mat)
  }
  return(colSums(mat))
}

#' Row-wise means of a matrix or returning a 1-d vector
#'
#' @param mat matrix or vector
#'
#' @returns rowMeans(mat) if matrix, mat if vector
#' @noRd
rowMeans_wrapper <- function(mat) {
  # allows a vector to be passed in (returns itself)
  # otherwise row means
  if (!('matrix' %in% class(mat))) {
    return(mat)
  }
  return(rowMeans(mat))
}

#' Pairwise cosine similarity between rows or columns of matrices
#'
#' @param mat1 matrix, first matrix for comparison
#' @param mat2 matrix, second matrix for comparison
#' @param name1 string, to name rows or cols of similarity matrix
#' @param name2 string, to name rows or cols of similarity matrix
#' @param which string, one of c("rows","cols")
#'
#' @return matrix
#' @export
pairwise_sim <- function(
    mat1, mat2,
    name1 = NULL,
    name2 = NULL,
    which = 'cols'
) {
  row_names = colnames(mat1)
  col_names = colnames(mat2)

  if (which == 'cols') {
    mat1 = t(mat1)
    mat2 = t(mat2)
  }

  if (ncol(mat1) != ncol(mat2)) {
    overlap_dim = ifelse(which == "rows","cols","rows")
    stop(paste0(
      "Different number of ", overlap_dim, ": ",
      ncol(mat1), " != ", ncol(mat2)
    ))
  }

  sim_mat = do.call(rbind, lapply(1:nrow(mat1), function(row_mat1) {
    sapply(1:nrow(mat2), function(row_mat2) {
      lsa::cosine(mat1[row_mat1,], mat2[row_mat2,])
    })
  }))

  if (!is.null(name1)) {
    rownames(sim_mat) = paste0(name1, 1:nrow(sim_mat))
  } else {
    rownames(sim_mat) = row_names
  }

  if (!is.null(name2)) {
    colnames(sim_mat) = paste0(name2, 1:ncol(sim_mat))
  } else {
    colnames(sim_mat) = col_names
  }

  return(sim_mat)
}

#' Assign factors to a reference based on cosine similarities
#'
#' @param sim_mat similarity matrix between learned and reference factors
#'
#' @return matrix
#' @export
assign_factors <- function(sim_mat, keep_all = FALSE) {
  reassignment <- RcppHungarian::HungarianSolver(-1 * sim_mat)
  rows = reassignment$pairs[,1]
  cols = reassignment$pairs[,2]
  if (keep_all) {
    for (row in setdiff(1:nrow(sim_mat), rows)) {
      rows <- c(rows, row)
    }
    for (col in setdiff(1:ncol(sim_mat), cols)) {
      cols <- c(cols, col)
    }
  }
  reassigned_sim_mat <- sim_mat[rows, cols]

  if (nrow(sim_mat) == 1 | ncol(sim_mat) == 1) {
    reassigned_sim_mat = matrix(reassigned_sim_mat)
    colnames(reassigned_sim_mat) = colnames(sim_mat)[cols]
    rownames(reassigned_sim_mat) = rownames(sim_mat)[rows]
  }

  reassigned_sim_mat
}

#' Returns indices to split dataset into two parts
#'
#' @param N integer, sample size
#' @param prop numeric, proportion of data in first part of data
#' @param seed integer, random seed or NULL
#' @param force_second vector, indices to force in the second part of data
#'
#' @returns vector, indices for first part of data
#' @noRd
split_dat <- function(N, prop, seed = NULL, force_second = c()) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n = round(prop*N) # number in first part
  idx <- 1:N

  # random split
  idx_first = sample(which(!(idx %in% force_second)), size = n)

  return(idx_first)
}

#' NMF wrapper
#'
#' @param M matrix DxN, NMF input matrix
#' @param rank integer K << D, N, latent rank
#' @param nrun integer, number of NMF runs (passed to NMF::nmf)
#' @param seed integer, random seed (passed to NMF::nmf)
#' @param method string, NMF method (passed to NMF::nmf)
#'
#' @returns NMF result, value from NMF::nmf
#' @export
#'
#' @import NMF
nmf_wrapper <- function(M, rank, nrun = 5, seed = NULL, method = 'brunet', ...) {
  problem_rows <- rowSums(M) == 0
  problem_cols <- colSums(M) == 0
  if (sum(problem_rows) == 0 & sum(problem_cols) == 0){
    nmf_res <- NMF::nmf(M, rank = rank, nrun = nrun, seed = seed, method = method, ...)
  } else {
    # NMF fails when one row is all 0
    # run NMF excluding those rows, then add back into factors with 0s
    M_sub <- M[!problem_rows,!problem_cols]
    nmf_res_sub <- NMF::nmf(M_sub, rank = rank, nrun = nrun, seed = seed, method = method, ...)
    nmf_res <- expand_problems(nmf_res_sub, problem_rows, problem_cols)
  }
  return(nmf_res)
}

#' Expand NMF problem rows and columns
#'
#' @param nmf_res_sub NMF output on subset of data without problem rows
#' @param problem_rows indices of problem rows
#' @param problem_cols incides of problem cols
#'
#' @returns Updated NMF output
#' @noRd
expand_problems <- function(nmf_res_sub, problem_rows, problem_cols) {
  What_sub <- nmf_res_sub@fit@W
  Hhat_sub <- nmf_res_sub@fit@H

  # put 0 back in factors for problem rows
  What <- matrix(0, nrow = nrow(What_sub) + sum(problem_rows), ncol = ncol(What_sub))
  What[!problem_rows,] <- What_sub

  stopifnot(all((rowSums(What) == 0) == problem_rows))
  nmf_res_sub@fit@W = What

  # put 0 back in weights for problem cols
  Hhat <- matrix(0, nrow = nrow(Hhat_sub), ncol = ncol(Hhat_sub) + sum(problem_cols))
  Hhat[,!problem_cols] <- Hhat_sub

  stopifnot(all((colSums(Hhat) == 0) == problem_cols))
  nmf_res_sub@fit@H = Hhat

  return(nmf_res_sub)
}

#' Extract and process NMF estimates
#'
#' @param nmf_res result from NMF::nmf
#' @param reference_P matrix DxK, reference factor matrix, if available
#'
#' @returns list with estimated factors (P), factor loadings (C), and alignment
#' to reference_P if provided
#' @export
extract_nmf_info <- function(nmf_res, reference_P = NULL) {
  Phat <- nmf_res@fit@W
  Chat <- nmf_res@fit@H

  N <- ncol(Phat)

  # rescale so E is on the scale of mutation counts
  # and columns of P sum to 1
  Chat <- sweep(Chat, 1, colSums(Phat), `*`)
  Phat <- sweep(Phat, 2, colSums(Phat), `/`)

  # reorder estimated signatures to match reference
  minsim = NULL
  reassigned = NULL
  if (!is.null(reference_P)) {
    sim <- pairwise_sim(reference_P, Phat, name2 = 'est')
    colnames(Phat) <- paste0('est', 1:N)
    rownames(Chat) <- paste0('est', 1:N)
    reassigned <- assign_factors(sim)
    minsim = min(diag(reassigned))
    Phat <- Phat[, colnames(reassigned)]
    Chat <- Chat[colnames(reassigned), ]

    colnames(Phat) <- rownames(reassigned)
    rownames(Chat) <- rownames(reassigned)
  } else {
    print("no reference")
  }

  return(list(
    P = Phat,
    C = Chat,
    minsim = minsim,
    reassigned = reassigned
  ))
}

bootstrap_data <- function(M, Tr) {
  G = ncol(M)
  idx = sample(1:G, G, replace = TRUE)

  return(list(
    M = M[,idx],
    Tr = Tr[idx]
  ))
}
