#' All Data algorithm
#' @description All Data algorithm to estimate ATE on latent factor-modeled outcomes.
#' Fits NMF on all data, then estimates ATE from estimated latent outcomes.
#' \strong{Subject to measurement interference and not recommended by the authors.}
#'
#' @param M \code{DxN} matrix, observed outcomes for \code{N} samples
#' @param Tr boolean vector of length \code{N}, treatment assignment for \code{N} samples
#' @param rank integer, latent rank \code{K}
#' @param reference_P \code{DxK} matrix, optional latent factor reference matrix.
#' Results will be aligned to these factors if provided.
#' @param seed integer, random seed passed to \link[NMF]{nmf}
#' @param nrun integer, number of NMF runs to perform, passed to \link[NMF]{nmf}
#' @param method string, specification of the NMF algorithm, passed to \link[NMF]{nmf}.
#' Default is \code{"brunet"} to optimize Poisson likelihood.
#' @param bootstrap boolean, whether to employ bootstrap resampling for
#' uncertainty quantification; default is \code{FALSE}
#' @param bootstrap_reps integer, number of bootstrap samples to use; default is \code{500}
#' @param bootstrap_filename string, file name to save bootstrap estimates
#' and aligned factor matrices
#'
#' @return If \code{bootstrap = TRUE}, returns the output of \link{bootstrap_wrapper}. Otherwise:
#' \itemize{
#'   \item \code{ATE}: estimated average treatment effect, vector of length \code{K = rank}
#'   \item \code{sim_mat}: cosine similarity matrix \code{K x K} between estimated (\code{Phat}) and reference (\code{reference_P}) factors if provided
#'   \item \code{Chat}: estimated latent causal outcome matrix \code{K x N} corresponding to assigned treatment levels, one column per sample
#'   \item \code{Phat}: estimated latent factor matrix \code{D x K}
#' }
#' @export
all_data <- function(
  M, Tr, rank,
  reference_P = NULL,
  seed = NULL, nrun = 5, method = 'brunet',
  bootstrap = FALSE, bootstrap_reps = 500,
  bootstrap_filename = "all_data"
) {
  if (mean(Tr) == 1 | mean(Tr) == 0) {
    stop(glue("No diversity in treatment, all Tr = {mean(Tr)}"))
  }
  if (bootstrap) {
    return(bootstrap_wrapper(
      all_data, bootstrap_filename,
      M, Tr, rank,
      seed, nrun, method,
      reference_P = reference_P,
      bootstrap_reps = bootstrap_reps
    ))
  }
  # NMF on all data
   nmf_res <- nmf_wrapper(
    M, rank = rank, nrun = nrun, seed = seed, method = method
  )
  # Rescale and align to reference
  nmf_res <- extract_nmf_info(nmf_res, reference_P = reference_P)

  # Difference of means estimate of ATE on all data
  ATE <- rowMeans_wrapper(nmf_res$C[,Tr==1]) -
         rowMeans_wrapper(nmf_res$C[,Tr==0])

  out <- list(
    ATE = ATE,
    sim_mat = nmf_res$reassigned,
    Chat = nmf_res$C,
    Phat = nmf_res$P
  )
  class(out) <- c("causalLFO_result", "list")
  return(out)
}

#' Random Split algorithm
#' @description Random Split algorithm to estimate ATE on latent factor-modeled outcomes.
#' Fits NMF on a subset of data, a Poisson non-negative linear model on the rest with
#' fixed factors, then estimates ATE from estimated latent outcomes in the second subset.
#' \strong{Subject to measurement interference and not recommended by the authors.}
#'
#' @param M \code{DxN} matrix, observed outcomes for \code{N} samples
#' @param Tr boolean vector of length \code{N}, treatment assignment for \code{N} samples
#' @param rank integer, latent rank \code{K}
#' @param reference_P \code{DxK} matrix, optional latent factor reference matrix.
#' Results will be aligned to these factors if provided.
#' @param prop float, proportion to include in the first (NMF) part of data
#' @param force_second vector, indices to force into the second (NNLM) part of data
#' @param seed integer, random seed passed to \link[NMF]{nmf}
#' @param nrun integer, number of NMF runs to perform, passed to \link[NMF]{nmf}
#' @param method string, specification of the NMF algorithm, passed to \link[NMF]{nmf}.
#' Default is \code{"brunet"} to optimize Poisson likelihood.
#' @param bootstrap boolean, whether to employ bootstrap resampling for
#' uncertainty quantification; default is \code{FALSE}
#' @param bootstrap_reps integer, number of bootstrap samples to use; default is \code{500}
#' @param bootstrap_filename string, file name to save bootstrap estimates
#' and aligned factor matrices
#'
#' @return If \code{bootstrap = TRUE}, returns the output of \link{bootstrap_wrapper}. Otherwise:
#' \itemize{
#'   \item \code{ATE}: estimated average treatment effect, vector of length \code{K = rank}
#'   \item \code{sim_mat}: cosine similarity matrix \code{K x K} between estimated (\code{Phat}) and reference (\code{reference_P}) factors if provided
#'   \item \code{Chat}: estimated latent causal outcome matrix \code{K x N} corresponding to assigned treatment levels, one column per sample
#'   \item \code{Phat}: estimated latent factor matrix \code{D x K}
#' }
#'
#' @export
random_split <- function(
  M, Tr, rank,
  reference_P = NULL,
  prop = 0.5, force_second = c(),
  seed = NULL, nrun = 5, method = "brunet",
  bootstrap = FALSE, bootstrap_reps = 500,
  bootstrap_filename = "random_split"
) {
  if (mean(Tr) == 1 | mean(Tr) == 0) {
    stop(glue("No diversity in treatment, all Tr = {mean(Tr)}"))
  }
  if (bootstrap) {
    return(bootstrap_wrapper(
      random_split, bootstrap_filename,
      M, Tr, rank,
      seed, nrun, method,
      reference_P = reference_P,
      bootstrap_reps = bootstrap_reps,
      prop = prop,
      force_second = force_second
    ))
  }
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

  out <- list(
    ATE = ATE,
    sim_mat = nmf_res$reassigned,
    Chat = Chat_expanded,
    split_idx = split_idx,
    Phat = nmf_res$P
  )
  class(out) <- c("causalLFO_result", "list")
  return(out)
}

#' Impute algorithm
#' @description Impute algorithm to estimate ATE on latent factor-modeled outcomes.
#' Imputes counterfactual outcomes under Poisson distributional assumptions,
#' fits NMF on observed data, a Poisson non-negative linear model on imputed data,
#' then estimates ATE as the mean difference in estimated latent outcomes between treated and untreated.
#' \strong{Intended as an ablation of \link{impute_and_stabilize} and not recommended by the authors.}
#'
#' @param M \code{DxN} matrix, observed outcomes for \code{N} samples
#' @param Tr boolean vector of length \code{N}, treatment assignment for \code{N} samples
#' @param rank integer, latent rank \code{K}
#' @param reference_P \code{DxK} matrix, optional latent factor reference matrix.
#' Results will be aligned to these factors if provided.
#' @param seed integer, random seed passed to \link[NMF]{nmf}
#' @param nrun integer, number of NMF runs to perform, passed to \link[NMF]{nmf}
#' @param method string, specification of the NMF algorithm, passed to \link[NMF]{nmf}.
#' Default is \code{"brunet"} to optimize Poisson likelihood.
#' @param bootstrap boolean, whether to employ bootstrap resampling for
#' uncertainty quantification; default is \code{FALSE}
#' @param bootstrap_reps integer, number of bootstrap samples to use; default is \code{500}
#' @param bootstrap_filename string, file name to save bootstrap estimates
#' and aligned factor matrices
#'
#' @return If \code{bootstrap = TRUE}, returns the output of \link{bootstrap_wrapper}. Otherwise:
#' \itemize{
#'   \item \code{ATE}: estimated average treatment effect, vector of length \code{K = rank}
#'   \item \code{sim_mat}: cosine similarity matrix \code{K x K} between estimated (\code{Phat}) and reference (\code{reference_P}) factors if provided
#'   \item \code{Chat}: estimated latent causal outcome matrix \code{K x N} corresponding to assigned treatment levels, one column per sample
#'   \item \code{Chat_imputed}: estimated latent causal outcome matrix \code{K x N} corresponding to counterfactual treatment levels, one column per sample
#'   \item \code{Phat}: estimated latent factor matrix \code{D x K}
#' }
#'
#' @export
impute <- function(
    M, Tr, rank, reference_P = NULL,
    seed = NULL, nrun = 5, method = "brunet",
    bootstrap = FALSE, bootstrap_reps = 500,
    bootstrap_filename = "impute"
) {
  if (mean(Tr) == 1 | mean(Tr) == 0) {
    stop(glue("No diversity in treatment, all Tr = {mean(Tr)}"))
  }
  if (bootstrap) {
    return(bootstrap_wrapper(
      impute, bootstrap_filename,
      M, Tr, rank,
      seed, nrun, method,
      reference_P = reference_P,
      bootstrap_reps = bootstrap_reps
    ))
  }
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

  out <- list(
    ATE = ATE_C,
    sim_mat = nmf_res$reassigned,
    Chat = Chat,
    Chat_imputed = Chatimp,
    Phat = nmf_res$P
  )
  class(out) <- c("causalLFO_result", "list")
  return(out)
}

#' Stabilize algorithm
#' @description Stabilize algorithm to estimate ATE on latent factor-modeled outcomes.
#' Fits NMF on untreated samples, a Poisson non-negative linear model on treated samples,
#' then estimates ATE using estimated latent outcomes.
#' \strong{Intended as an ablation of \link{impute_and_stabilize} and not recommended by the authors.}
#'
#' @param M \code{DxN} matrix, observed outcomes for \code{N} samples
#' @param Tr boolean vector of length \code{N}, treatment assignment for \code{N} samples
#' @param rank integer, latent rank \code{K}
#' @param reference_P \code{DxK} matrix, optional latent factor reference matrix.
#' Results will be aligned to these factors if provided.
#' @param seed integer, random seed passed to \link[NMF]{nmf}
#' @param nrun integer, number of NMF runs to perform, passed to \link[NMF]{nmf}
#' @param method string, specification of the NMF algorithm, passed to \link[NMF]{nmf}.
#' Default is \code{"brunet"} to optimize Poisson likelihood.
#' @param bootstrap boolean, whether to employ bootstrap resampling for
#' uncertainty quantification; default is \code{FALSE}
#' @param bootstrap_reps integer, number of bootstrap samples to use; default is \code{500}
#' @param bootstrap_filename string, file name to save bootstrap estimates
#' and aligned factor matrices
#'
#' @return If \code{bootstrap = TRUE}, returns the output of \link{bootstrap_wrapper}. Otherwise:
#' \itemize{
#'   \item \code{ATE}: estimated average treatment effect, vector of length \code{K = rank}
#'   \item \code{sim_mat}: cosine similarity matrix \code{K x K} between estimated (\code{Phat}) and reference (\code{reference_P}) factors if provided
#'   \item \code{Chat}: estimated latent causal outcome matrix \code{K x N} corresponding to assigned treatment levels, one column per sample
#'   \item \code{Phat}: estimated latent factor matrix \code{D x K}
#' }
#'
#' @export
stabilize <- function(
    M, Tr, rank, reference_P = NULL,
    seed = NULL, nrun = 5, method = "brunet",
    bootstrap = FALSE, bootstrap_reps = 500,
    bootstrap_filename = "stabilize"
) {
  if (mean(Tr) == 1 | mean(Tr) == 0) {
    stop(glue("No diversity in treatment, all Tr = {mean(Tr)}"))
  }
  if (bootstrap) {
    return(bootstrap_wrapper(
      stabilize, bootstrap_filename,
      M, Tr, rank,
      seed, nrun, method,
      reference_P = reference_P,
      bootstrap_reps = bootstrap_reps
    ))
  }
  G <- ncol(M)
  Chat <- matrix(nrow = rank, ncol = ncol(M))

  # NMF on M0
  nmf_res <- nmf_wrapper(
    M[,Tr==0], rank = rank, nrun = nrun, seed = seed, method = method
  )
  # Rescale and align to reference
  nmf_res <- extract_nmf_info(nmf_res, reference_P = reference_P)
  Chat[,Tr==0] <- nmf_res$C

  # NNLM on M1
  Chat[,Tr==1] <- poisson_nnlm(M[,Tr==1], fixed_P = nmf_res$P)

  # DM estimator on Chat
  ATE <- rowMeans_wrapper(Chat[,Tr==1]) -
         rowMeans_wrapper(Chat[,Tr==0])

  out <- list(
    ATE = ATE,
    sim_mat = nmf_res$reassigned,
    Chat = Chat,
    Phat = nmf_res$P
  )
  class(out) <- c("causalLFO_result", "list")
  return(out)
}


#' Impute and Stabilize algorithm
#' @description Impute and Stabilize algorithm to estimate ATE on latent factor-modeled outcomes.
#' Imputes counterfactual outcomes under Poisson distributional assumptions,
#' fits NMF on untreated data (mix of observed and imputed), a Poisson non-negative linear model on treated data,
#' then estimates ATE as the mean difference in estimated latent outcomes between treated and untreated.
#'
#' @param M \code{DxN} matrix, observed outcomes for \code{N} samples
#' @param Tr boolean vector of length \code{N}, treatment assignment for \code{N} samples
#' @param rank integer, latent rank \code{K}
#' @param reference_P \code{DxK} matrix, optional latent factor reference matrix.
#' Results will be aligned to these factors if provided.
#' @param seed integer, random seed passed to \link[NMF]{nmf}
#' @param nrun integer, number of NMF runs to perform, passed to \link[NMF]{nmf}
#' @param method string, specification of the NMF algorithm, passed to \link[NMF]{nmf}.
#' Default is \code{"brunet"} to optimize Poisson likelihood.
#' @param bootstrap boolean, whether to employ bootstrap resampling for
#' uncertainty quantification; default is \code{FALSE}
#' @param bootstrap_reps integer, number of bootstrap samples to use; default is \code{500}
#' @param bootstrap_filename string, file name to save bootstrap estimates
#' and aligned factor matrices
#'
#' @return If \code{bootstrap = TRUE}, returns the output of \link{bootstrap_wrapper}. Otherwise:
#' \itemize{
#'   \item \code{ATE}: estimated average treatment effect, vector of length \code{K = rank}
#'   \item \code{sim_mat}: cosine similarity matrix \code{K x K} between estimated (\code{Phat}) and reference (\code{reference_P}) factors if provided
#'   \item \code{Chat}: estimated latent causal outcome matrix \code{K x N} corresponding to assigned treatment levels, one column per sample
#'   \item \code{Chat_imputed}: estimated latent causal outcome matrix \code{K x N} corresponding to counterfactual treatment levels, one column per sample
#'   \item \code{Phat}: estimated latent factor matrix \code{D x K}
#' }
#'
#' @export
impute_and_stabilize <- function(
    M, Tr, rank, reference_P = NULL,
    seed = NULL, nrun = 5, method = "brunet",
    bootstrap = FALSE, bootstrap_reps = 500,
    bootstrap_filename = "impute_and_stabilize"
) {
  if (mean(Tr) == 1 | mean(Tr) == 0) {
    stop(glue("No diversity in treatment, all Tr = {mean(Tr)}"))
  }
  if (bootstrap) {
    return(bootstrap_wrapper(
      impute_and_stabilize, bootstrap_filename,
      M, Tr, rank,
      seed, nrun, method,
      reference_P = reference_P,
      bootstrap_reps = bootstrap_reps
    ))
  }
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

  out <- list(
    ATE = ATE_C,
    sim_mat = nmf_res$reassigned,
    Chat = combine_mat(Chat0, Chat1, Tr),
    Chat_imputed = combine_mat(Chat0, Chat1, 1-Tr),
    Phat = nmf_res$P
  )
  out <- structure(out, class = c("causalLFO_result"))
  return(out)
}

#' Bootstrap wrapper for causalLFO
#' @description Performs bootstrap resampling and applies the provided algorithm to each
#' sample. Aligns estimated factor matrices and returns bootsrapped mean, se, and 95% confidence
#' interval for the ATE. Called internally by causalLFO algorithms (\link{all_data}, \link{random_split},
#' \link{impute}, \link{stabilize}, and \link{impute_and_stabilize}) when \code{bootstrap = TRUE}.
#'
#' @param algorithm function, a \code{causalLFO} algorithm (e.g., \link{random_split})
#' @param filename string, file name to save bootstrap estimates and aligned
#' factor matrices
#' @param M \code{DxN} matrix, observed outcomes for \code{N} samples
#' @param Tr boolean vector of length \code{N}, treatment assignment for \code{N} samples
#' @param rank integer, latent rank \code{K}
#' @param reference_P \code{DxK} matrix, optional latent factor reference matrix.
#' Results will be aligned to these factors if provided.
#' @param bootstrap_reps integer, number of bootstrap samples to use; default is \code{500}
#' @param ... other parameters to be passed to \code{algorithm}, e.g., \code{"prop"}
#' for the \link{random_split} algorithm
#'
#' @return A \code{list} with the following elements:
#' \itemize{
#'   \item \code{mean}: bootstrap mean ATE, vector of length \code{K = rank}
#'   \item \code{se}: bootstrap standard error, vector of length \code{K = rank}
#'   \item \code{lower95}: lower end of 95\% confidence interval via bootstrap quantiles, vector of length \code{K = rank}
#'   \item \code{upper95}: upper end of 95\% confidence interval via bootstrap quantiles, vector of length \code{K = rank}
#'   \item \code{Phat}: consensus factor matrix, elementwise mean across aligned bootstrap replicates
#'   \item \code{Chat}: consensus factor loadings matrix, generated via NNLM on \code{Phat} and \code{M}
#'   \item \code{est_file}: string, path to file containing all ATE estimates
#'   \item \code{all_Ps_file}: string, path to file containing all aligned factor matrices
#' }
#'
#' @import glue
#' @export
bootstrap_wrapper <- function(
    algorithm, filename,
    M, Tr, rank,
    seed = NULL, nrun = 5, method = 'brunet',
    reference_P = NULL,
    bootstrap_reps = 500,
    ...
) {
  external_reference <- !is.null(reference_P)
  estimates <- matrix(nrow = bootstrap_reps, ncol = rank)
  aligned_Ps <- list()
  split_idxs <- list()
  for (rep in 1:bootstrap_reps) {
    boot_data <- bootstrap_data(M, Tr)
    while(mean(boot_data$Tr) == 0 | mean(boot_data$Tr) == 1) {
      print("trying bootstrap again")
      boot_data <- bootstrap_data(M, Tr)
    }

    # Call the algorithm extra arguments `...`
    est <- algorithm(
      M = boot_data$M,
      Tr = boot_data$Tr,
      rank = rank,
      seed = seed,
      nrun = nrun,
      method = method,
      reference_P = reference_P,
      ... # e.g., prop, force_second, etc.
    )
    estimates[rep,] <- est$ATE
    aligned_Ps[[rep]] <- est$Phat
    if ("split_idx" %in% names(est)) {
      split_idxs[[rep]] <- est$split_idx
    }
    if (!external_reference) {
      Pbar <- Reduce(`+`, aligned_Ps)/length(aligned_Ps)
      reference_P = Pbar
    }

    if (rep == 1) {
      # create new file
      colnames(estimates) <- colnames(reference_P)
      write.csv(
        estimates[rep, , drop = FALSE],
        file = glue::glue("{filename}.csv"),
        row.names = FALSE
      )
    } else {
      # add new line, don't overwrite file
      write.table(
        estimates[rep, , drop = FALSE],
        file = glue::glue("{filename}.csv"),
        sep = ",", append = TRUE,
        col.names = FALSE, row.names = FALSE
      )
    }
  }
  saveRDS(aligned_Ps, file = glue::glue("{filename}_aligned_Ps.rds"))

  # bootstrap samples have different individuals/orders, can't average Chats
  # instead, back-calculate from bootstrap Phat after the fact
  bootstrap_Phat = Reduce(`+`, aligned_Ps)/length(aligned_Ps)
  bootstrap_Chat = poisson_nnlm(M, fixed_P = bootstrap_Phat)

  out = list(
    mean = apply(estimates, 2, mean, na.rm = TRUE),
    se = apply(estimates, 2, sd, na.rm = TRUE),
    lower95 = apply(estimates, 2, function(col){quantile(col, 0.025)}),
    upper95 = apply(estimates, 2, function(col){quantile(col, 0.975)}),
    Phat = bootstrap_Phat,
    Chat = bootstrap_Chat,
    est_file = glue::glue("{filename}.csv"),
    all_Ps_file = glue::glue("{filename}_aligned_Ps.rds")
  )
  if ("split_idx" %in% names(est)) {
    out[["split_idxs"]] <- split_idxs
  }
  out <- structure(out, class = c("causalLFO_bootstrap_result"))
  saveRDS(out, file = glue::glue("{filename}_res.rds"))
  return(out)
}
