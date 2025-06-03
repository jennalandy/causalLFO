test_that("nnlm works", {
  true_P = matrix(rpois(96*3, lambda = 19), nrow = 96)
  true_C = matrix(rpois(3*10, lambda = 1), nrow = 3)
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- true_P %*% true_C[,i]
  }

  Chat = poisson_nnlm(M, fixed_P = true_P)
  rmse = sqrt(sum((Chat - true_C)**2))

  expect_true(rmse < 0.1)
})

test_that("nnlm works on sparse", {
  true_P = matrix(rpois(96*3, lambda = 1), nrow = 96)
  true_C = matrix(rpois(3*10, lambda = 1), nrow = 3)
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- true_P %*% true_C[,i]
  }

  Chat = poisson_nnlm(M, fixed_P = true_P)
  rmse = sqrt(sum((Chat - true_C)**2))

  expect_true(rmse < 0.1)
})
