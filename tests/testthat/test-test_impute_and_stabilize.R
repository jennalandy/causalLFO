test_that("impute_and_stabilize works", {
  true_P = matrix(rexp(96*3, rate = 1), nrow = 96)
  true_P = sweep(true_P, 2, colSums(true_P), '/')
  true_C = matrix(rexp(3*10, rate = 0.01), nrow = 3)
  Tr = sample(c(0, 1), 10, replace = TRUE)
  true_C[1, Tr == 1] <- true_C[1, Tr == 1] + 30
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- rpois(96, lambda = true_P %*% true_C[,i])
  }

  res <- impute_and_stabilize(M, Tr, rank = 3, reference_P = true_P)
  expect_equal(length(res$ATE), 3)
  expect_equal(dim(res$Phat), c(96, 3))
  expect_equal(dim(res$Chat), c(3, 10))
})

test_that("bootstrap_wrapper for impute_and_stabilize works called internally", {
  true_P = matrix(rexp(96*3, rate = 1), nrow = 96)
  true_P = sweep(true_P, 2, colSums(true_P), '/')
  true_C = matrix(rexp(3*10, rate = 0.01), nrow = 3)
  Tr = sample(c(0, 1), 10, replace = TRUE)
  true_C[1, Tr == 1] <- true_C[1, Tr == 1] + 30
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- rpois(96, lambda = true_P %*% true_C[,i])
  }

  res <- impute_and_stabilize(
    M, Tr, rank = 3, reference_P = true_P,
    bootstrap = TRUE, bootstrap_file = "impute_and_stabilize",
    bootstrap_reps = 10
  )
  expect_equal(length(res$mean), 3)
  expect_equal(length(res$se), 3)
  expect_equal(length(res$lower95), 3)
  expect_equal(length(res$upper95), 3)
  expect_true(res$lower95 < res$upper95)
  expect_equal(dim(res$Phat), c(96, 3))
  expect_equal(dim(res$Chat), c(3, 10))
  expect_equal(res$est_file, "impute_and_stabilize.csv")
  expect_equal(res$all_Ps_file, "impute_and_stabilize_aligned_Ps.rds")

  est <- read.csv(res$est_file)
  expect_equal(nrow(est), 10)
})
