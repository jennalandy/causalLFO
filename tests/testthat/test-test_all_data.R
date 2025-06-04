test_that("all_data works", {
  true_P = matrix(rexp(96*3, rate = 1), nrow = 96)
  true_P = sweep(true_P, 2, colSums(true_P), '/')
  true_C = matrix(rexp(3*10, rate = 0.01), nrow = 3)
  Tr = sample(c(0, 1), 10, replace = TRUE)
  true_C[1, Tr == 1] <- true_C[1, Tr == 1] + 30
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- rpois(96, lambda = true_P %*% true_C[,i])
  }

  res <- all_data(M, Tr, rank = 3, reference_P = true_P)
  expect_equal(length(res$ATE), 3)
  expect_equal(dim(res$Phat), c(96, 3))
  expect_equal(dim(res$Chat), c(3, 10))
})

test_that("bootstrap_wrapper of all_data works called externally", {
  true_P = matrix(rexp(96*3, rate = 1), nrow = 96)
  true_P = sweep(true_P, 2, colSums(true_P), '/')
  true_C = matrix(rexp(3*10, rate = 0.01), nrow = 3)
  Tr = sample(c(0, 1), 10, replace = TRUE)
  true_C[1, Tr == 1] <- true_C[1, Tr == 1] + 30
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- rpois(96, lambda = true_P %*% true_C[,i])
  }

  res <- bootstrap_wrapper(
    all_data, "all_data_ext",
    M, Tr, rank = 3, reference_P = true_P,
    bootstrap_reps = 5
  )
  expect_equal(length(res$mean), 3)
  expect_equal(length(res$se), 3)
  expect_equal(length(res$lower95), 3)
  expect_equal(length(res$upper95), 3)
  expect_true(all(res$lower95 < res$upper95))
  expect_equal(dim(res$Phat), c(96, 3))
  expect_equal(dim(res$Chat), c(3, 10))
  expect_equal(res$est_file, "all_data_ext.csv")
  expect_equal(res$all_Ps_file, "all_data_ext.rds")

  est <- read.csv(res$est_file)
  expect_equal(nrow(est), 5)
})

test_that("bootstrap_wrapper of all_data works called internally", {
  true_P = matrix(rexp(96*3, rate = 1), nrow = 96)
  true_P = sweep(true_P, 2, colSums(true_P), '/')
  true_C = matrix(rexp(3*10, rate = 0.01), nrow = 3)
  Tr = sample(c(0, 1), 10, replace = TRUE)
  true_C[1, Tr == 1] <- true_C[1, Tr == 1] + 30
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- rpois(96, lambda = true_P %*% true_C[,i])
  }

  res <- all_data(
    M, Tr, rank = 3, reference_P = true_P,
    bootstrap = TRUE, bootstrap_file = "all_data_int",
    bootstrap_reps = 5
  )
  expect_equal(length(res$mean), 3)
  expect_equal(length(res$se), 3)
  expect_equal(length(res$lower95), 3)
  expect_equal(length(res$upper95), 3)
  expect_true(all(res$lower95 < res$upper95))
  expect_equal(dim(res$Phat), c(96, 3))
  expect_equal(dim(res$Chat), c(3, 10))
  expect_equal(res$est_file, "all_data_int.csv")
  expect_equal(res$all_Ps_file, "all_data_int_aligned_Ps.rds")

  est <- read.csv(res$est_file)
  expect_equal(nrow(est), 5)
})
