#' All data Algorithm
#'
#' @param M
#' @param Tr
#' @param rank
#' @param reference_P
#' @param seed
#' @param nrun
#' @param method
#' @param force_second
#'
#' @returns
#' @export
#'
#' @examples
all_data <- function(
  M, Tr, rank,
  reference_P = NULL,
  seed = NULL, nrun = 5, method = 'brunet',
  force_second = NULL
) {
  # NMF on all data
   nmf_res <- nmf_wrapper(
    M, rank = rank, nrun = nrun, seed = seed, method = method
  )
  # Rescale and align to reference
  nmf_res <- extract_nmf_info(nmf_res, reference_P = reference_P)

  # Difference of means estimate of ATE on all data
  ATE <- rowMeans_wrapper(nmf_res$C[,Tr==1]) -
         rowMeans_wrapper(nmf_res$C[,Tr==0])

  return(list(
    ATE = ATE,
    sim_mat = nmf_res$reassigned,
    Chat = nmf_res$C,
    Phat = nmf_res$P
  ))
}

#' Random split algorithm
#'
#' @param M
#' @param Tr
#' @param rank
#' @param reference_P
#' @param prop
#' @param seed
#' @param nrun
#' @param force_second
#'
#' @returns
#' @export
#'
#' @examples
random_split <- function(
  M, Tr, rank,
  reference_P = NULL,
  prop = 0.5,
  seed = NULL, nrun = 5, method = "brunet",
  force_second = c()
) {
  split_idx = split_dat(ncol(M), prop, seed = seed, force_second = force_second)

  # Need at least one treated and one untreated in the second half
  while (mean(Tr[-split_idx]) == 1 | mean(Tr[-split_idx]) == 0) {
    print("trying split again")
    split_idx = split_dat(ncol(M), prop, seed = seed, force_second = force_second)
  }

  # NMF on first part of data
  nmf_res <- nmf_wrapper(
    M[,split_idx], rank = rank, nrun = nrun, seed = seed, method = method
  )
  # Rescale and align to reference
  nmf_res <- extract_nmf_info(nmf_res, reference_P = reference_P)

  # NNLM with fixed P to estimate C on second part of data
  Chat <- poisson_nnlm(M[,-split_idx], fixed_P = nmf_res$P)

  # Difference of means estimate of ATE on second part of data
  ATE <- rowMeans_wrapper(Chat[,Tr[-split_idx]==1]) -
         rowMeans_wrapper(Chat[,Tr[-split_idx]==0])

  Chat_expanded <- matrix(NA, nrow = rank, ncol = ncol(M))
  Chat_expanded[,-split_idx] <- Chat

  return(list(
    ATE = ATE,
    sim_mat = nmf_res$reassigned,
    Chat = Chat_expanded,
    split_idx = split_idx,
    Phat = nmf_res$P
  ))
}

#' Impute algorithm
#'
#' @param M
#' @param Tr
#' @param N
#' @param reference_P
#' @param seed
#' @param nrun
#' @param method
#' @param force_second
#'
#' @returns
#' @export
#'
#' @examples
impute <- function(
    M, Tr, rank, reference_P = NULL,
    seed = NULL, nrun = 5, method = "brunet",
    force_second = NULL
) {
  G <- ncol(M)

  # Perform variance stabilizing
  M_star <- sqrt(M)

  # Impute potential outcomes on the M_star (VST of observable) scale
  ATE_M_star <- rowMeans_wrapper(M_star[,Tr==1]) -
    rowMeans_wrapper(M_star[,Tr==0])

  # Impute potential outcomes
  M1_star <- do.call(cbind, lapply(1:G, function(g) {
    M_star[,g] + ATE_M_star*(1-Tr[g])
  }))
  M0_star <- do.call(cbind, lapply(1:G, function(g) {
    M_star[,g] - ATE_M_star*Tr[g]
  }))

  # Back-transform potential outcomes
  M1 <- M1_star**2
  M0 <- M0_star**2
  theoretical_correction = 0.25*(1 + 1/sum(Tr == 0) + 1/sum(Tr == 1))
  M1[,Tr == 0] <- M1_star[,Tr == 0]**2 + theoretical_correction
  M0[,Tr == 1] <- M0_star[,Tr == 1]**2 + theoretical_correction
  Mimp = combine_mat(M0, M1, 1-Tr)

  # NMF on M(T) observed
  nmf_res <- nmf_wrapper(
    M, rank = rank, nrun = nrun, seed = seed, method = method
  )
  # Rescale and align to reference
  nmf_res <- extract_nmf_info(nmf_res, reference_P = reference_P)
  Chat <- nmf_res$C

  # NNLM on M(1-T) imputed
  Chatimp <- poisson_nnlm(Mimp, fixed_P = nmf_res$P)

  # Construct C0, C1
  Chat0 <- matrix(nrow = nrow(nmf_res$C), ncol = ncol(nmf_res$C))
  Chat1 <- matrix(nrow = nrow(nmf_res$C), ncol = ncol(nmf_res$C))
  Chat0[,Tr == 0] <- nmf_res$C[,Tr == 0]
  Chat1[,Tr == 1] <- nmf_res$C[,Tr == 1]
  Chat0[,Tr == 1] <- Chatimp[,Tr == 1]
  Chat1[,Tr == 0] <- Chatimp[,Tr == 0]

  # mean ITE estimator on Chat
  ATE_C <- rowMeans_wrapper(Chat1 - Chat0)

  return(list(
    ATE = ATE_C,
    sim_mat = nmf_res$reassigned,
    Chat = Chat,
    Chat_imputed = Chatimp,
    Phat = nmf_res$P
  ))
}

#' Stabilize algorithm
#'
#' @param M
#' @param Tr
#' @param N
#' @param reference_P
#' @param seed
#' @param nrun
#' @param method
#' @param force_second
#'
#' @returns
#' @export
#'
#' @examples
stabilize <- function(
    M, Tr, rank, reference_P = NULL,
    seed = NULL, nrun = 5, method = "brunet",
    force_second = NULL
) {
  G <- ncol(M)
  Chat <- matrix(nrow = nrow(M), ncol = rank)

  # NMF on M0
  nmf_res <- nmf_wrapper(
    M[,Tr == 0], rank = rank, nrun = nrun, seed = seed, method = method
  )
  # Rescale and align to reference
  nmf_res <- extract_nmf_info(nmf_res, reference_P = reference_P)
  Chat[,Tr == 0] <- nmf_res$C

  # NNLM on M1
  Chat[,Tr == 1] <- poisson_nnlm(M[,Tr==1], fixed_P = nmf_res$P)

  # DM estimator on Chat
  ATE <- rowMeans_wrapper(Chat[,Tr==1]) -
         rowMeans_wrapper(Chat[,Tr==0])

  return(list(
    ATE = ATE,
    sim_mat = nmf_res$reassigned,
    Chat = Chat,
    Phat = nmf_res$P
  ))
}


#' Impute and Stabilize algorithm
#'
#' @param M
#' @param Tr
#' @param N
#' @param reference_P
#' @param seed
#' @param nrun
#' @param method
#' @param force_second
#'
#' @returns
#' @export
#'
#' @examples
impute_and_stabilize <- function(
    M, Tr, rank, reference_P = NULL,
    seed = NULL, nrun = 5, method = "brunet",
    force_second = NULL
) {
  G <- ncol(M)

  # Perform variance stabilizing
  M_star <- sqrt(M)

  # Impute potential outcomes on the M_star (VST of observable) scale
  ATE_M_star <- rowMeans_wrapper(M_star[,Tr==1]) -
                rowMeans_wrapper(M_star[,Tr==0])

  # Impute potential outcomes
  M1_star <- do.call(cbind, lapply(1:G, function(g) {
    M_star[,g] + ATE_M_star*(1-Tr[g])
  }))
  M0_star <- do.call(cbind, lapply(1:G, function(g) {
    M_star[,g] - ATE_M_star*Tr[g]
  }))

  # Back-transform potential outcomes
  M1 <- M1_star**2
  M0 <- M0_star**2
  theoretical_correction = 0.25*(1 + 1/sum(Tr == 0) + 1/sum(Tr == 1))
  M1[,Tr == 0] <- M1_star[,Tr == 0]**2 + theoretical_correction
  M0[,Tr == 1] <- M0_star[,Tr == 1]**2 + theoretical_correction

  # NMF on M0
  nmf_res <- nmf_wrapper(
    M0, rank = rank, nrun = nrun, seed = seed, method = method
  )
  # Rescale and align to reference
  nmf_res <- extract_nmf_info(nmf_res, reference_P = reference_P)
  Chat0 <- nmf_res$C

  # NNLM on M1
  Chat1 <- poisson_nnlm(M1, fixed_P = nmf_res$P)

  # mean ITE estimator on Chat
  ATE_C <- rowMeans_wrapper(Chat1 - Chat0)

  return(list(
    ATE = ATE_C,
    sim_mat = nmf_res$reassigned,
    Chat = combine_mat(Chat0, Chat1, Tr),
    Chat_imputed = combine_mat(Chat0, Chat1, 1-Tr),
    Phat = nmf_res$P
  ))
}
