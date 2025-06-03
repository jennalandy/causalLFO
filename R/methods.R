#' Poisson-likelihood non-negative linear model with fixed factor matrix using
#' Lee and Seung 1999-style multiplicative updates
#'
#' @param M matrix DxN, data matrix
#' @param fixed_P matrix DxK, fixed factor matrix
#' @param maxiter int, max number of gradient descent updates
#' @param tol float, tolerance for convergence
#'
#' @returns C, matrix KxN, factor loadings matrix
#' @export
poisson_nnlm <- function(M, fixed_P, maxiter = 10000, tol = 1e-6) {
  K <- nrow(M); G <- ncol(M); N <- ncol(fixed_P)

  # deterministic initialization of equal distribution of mutations across signatures
  C <- matrix(colSums_wrapper(M)/N, nrow = N, ncol = G, byrow = TRUE)

  eps <- 1e-10
  C_prev <- C
  for (iter in 1:maxiter) {
    numerator <- t(fixed_P) %*% (M / (fixed_P%*%C + eps))
    denominator <- matrix(colSums_wrapper(fixed_P), nrow = N, ncol = G)
    C <- C * numerator / (denominator + eps)
    delta <- sum(abs(C - C_prev)) / (sum(abs(C_prev)) + eps)
    if (delta < tol) {
      break
    }

    C_prev <- C
  }

  return(C)
}
