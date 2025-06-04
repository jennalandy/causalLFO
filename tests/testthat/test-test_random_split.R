test_that("random_split works", {
  true_P = matrix(rexp(96*3, rate = 1), nrow = 96)
  true_P = sweep(true_P, 2, colSums(true_P), '/')
  true_C = matrix(rexp(3*10, rate = 0.01), nrow = 3)
  Tr = sample(c(0, 1), 10, replace = TRUE)
  true_C[1, Tr == 1] <- true_C[1, Tr == 1] + 30
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- rpois(96, lambda = true_P %*% true_C[,i])
  }

  res <- random_split(M, Tr, rank = 3, reference_P = true_P)
  expect_equal(length(res$ATE), 3)
  expect_equal(dim(res$Phat), c(96, 3))
  expect_equal(dim(res$Chat), c(3, 10))
  expect_true(sum(is.na(res$Chat)) > 0)
})

test_that("random_split works with force_second", {
  true_P = matrix(rexp(96*3, rate = 1), nrow = 96)
  true_P = sweep(true_P, 2, colSums(true_P), '/')
  true_C = matrix(rexp(3*10, rate = 0.01), nrow = 3)
  Tr = sample(c(0, 1), 10, replace = TRUE)
  true_C[1, Tr == 1] <- true_C[1, Tr == 1] + 30
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- rpois(96, lambda = true_P %*% true_C[,i])
  }

  res <- random_split(M, Tr, rank = 3, reference_P = true_P, force_second = c(1,3))
  expect_equal(length(res$ATE), 3)
  expect_equal(dim(res$Phat), c(96, 3))
  expect_equal(dim(res$Chat), c(3, 10))
  expect_equal(sum(is.null(res$Chat[,3])), 0)
  expect_equal(sum(is.na(res$Chat[,3])), 0)
})

test_that("bootstrap_wrapper of random_split works called internally", {
  true_P = matrix(rexp(96*3, rate = 1), nrow = 96)
  true_P = sweep(true_P, 2, colSums(true_P), '/')
  true_C = matrix(rexp(3*10, rate = 0.01), nrow = 3)
  Tr = sample(c(0, 1), 10, replace = TRUE)
  true_C[1, Tr == 1] <- true_C[1, Tr == 1] + 30
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- rpois(96, lambda = true_P %*% true_C[,i])
  }

  res <- random_split(
    M, Tr, rank = 3, reference_P = true_P,
    bootstrap = TRUE, bootstrap_file = "random_split",
    bootstrap_reps = 5
  )
  expect_equal(length(res$mean), 3)
  expect_equal(length(res$se), 3)
  expect_equal(length(res$lower95), 3)
  expect_equal(length(res$upper95), 3)
  expect_true(all(res$lower95 < res$upper95))
  expect_equal(dim(res$Phat), c(96, 3))
  expect_equal(dim(res$Chat), c(3, 10))
  expect_equal(res$est_file, "random_split.csv")
  expect_equal(res$all_Ps_file, "random_split_aligned_Ps.rds")

  est <- read.csv(res$est_file)
  expect_equal(nrow(est), 5)
})

test_that("bootstrap_wrapper of random_split works with force_second called internally", {
  true_P = matrix(rexp(96*3, rate = 1), nrow = 96)
  true_P = sweep(true_P, 2, colSums(true_P), '/')
  true_C = matrix(rexp(3*10, rate = 0.01), nrow = 3)
  Tr = sample(c(0, 1), 10, replace = TRUE)
  true_C[1, Tr == 1] <- true_C[1, Tr == 1] + 30
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- rpois(96, lambda = true_P %*% true_C[,i])
  }

  res <- random_split(
    M, Tr, rank = 3, reference_P = true_P,
    force_second = c(1,3),
    bootstrap = TRUE, bootstrap_file = "random_split_force_second",
    bootstrap_reps = 5
  )
  expect_true("split_idxs" %in% names(res))
  # check indices in force_second are never in split_idx
  expect_equal(sum(sapply(res$split_idxs, function(split_idx) { 1 %in% split_idx })), 0)
  expect_equal(sum(sapply(res$split_idxs, function(split_idx) { 3 %in% split_idx })), 0)

  expect_equal(length(res$mean), 3)
  expect_equal(length(res$se), 3)
  expect_equal(length(res$lower95), 3)
  expect_equal(length(res$upper95), 3)
  expect_true(all(res$lower95 < res$upper95))
  expect_equal(dim(res$Phat), c(96, 3))
  expect_equal(dim(res$Chat), c(3, 10))
  expect_equal(res$est_file, "random_split_force_second.csv")
  expect_equal(res$all_Ps_file, "random_split_force_second_aligned_Ps.rds")

  est <- read.csv(res$est_file)
  expect_equal(nrow(est), 5)
})
