test_that("combine_mat works", {
  M0 = matrix(1, nrow = 2, ncol = 2)
  M1 = matrix(2, nrow = 2, ncol = 2)
  Tr = c(0, 1)
  M = combine_mat(M0, M1, Tr)
  expect_equal(M[,1], c(1, 1))
  expect_equal(M[,2], c(2, 2))

  Tr = c(1, 0)
  M = combine_mat(M0, M1, Tr)
  expect_equal(M[,1], c(2, 2))
  expect_equal(M[,2], c(1, 1))
})

test_that("colSums_wrapper works", {
  mat = matrix(1, nrow = 2, ncol = 2)
  expect_equal(colSums_wrapper(mat), c(2, 2))

  vec = c(1,1)
  expect_equal(colSums_wrapper(vec), vec)
})

test_that("rowMeans_wrapper works", {
  mat = matrix(1, nrow = 2, ncol = 2)
  expect_equal(rowMeans_wrapper(mat), c(1, 1))

  vec = c(1,1)
  expect_equal(rowMeans_wrapper(vec), vec)
})

test_that("split_data works", {
  split_idx <- split_dat(10, prop = 0.5)
  expect_equal(length(split_idx), 5)

  split_idx <- split_dat(10, prop = 0.2)
  expect_equal(length(split_idx), 2)

  split_idx <- split_dat(10, prop = 0.5, force_second = c(1, 5))
  expect_false(1 %in% split_idx)
  expect_false(5 %in% split_idx)
})

test_that("nmf_wrapper works", {
  M <- matrix(rpois(96*10, lambda = 10), nrow = 96)
  nmf_res <- nmf_wrapper(M, rank = 3, nrun = 1)

  expect_equal(nrow(nmf_res@fit@W), 96)
  expect_equal(ncol(nmf_res@fit@W), 3)
  expect_equal(nrow(nmf_res@fit@H), 3)
  expect_equal(ncol(nmf_res@fit@H), 10)
})

test_that("nmf_wrapper works with problem rows", {
  M <- matrix(rpois(96*10, lambda = 10), nrow = 96)
  M[10,] <- 0
  M[15,] <- 0
  nmf_res <- nmf_wrapper(M, rank = 3, nrun = 1)

  expect_equal(nrow(nmf_res@fit@W), 96)
  expect_equal(ncol(nmf_res@fit@W), 3)
  expect_equal(nrow(nmf_res@fit@H), 3)
  expect_equal(ncol(nmf_res@fit@H), 10)
})

test_that("extract_nmf_info works", {
  M <- matrix(rpois(96*10, lambda = 10), nrow = 96)
  nmf_res <- nmf_wrapper(M, rank = 3, nrun = 1)
  nmf_res <- extract_nmf_info(nmf_res)

  expect_equal(nrow(nmf_res$P), 96)
  expect_equal(ncol(nmf_res$P), 3)
  expect_equal(nrow(nmf_res$C), 3)
  expect_equal(ncol(nmf_res$C), 10)
  expect_null(nmf_res$minsim)
  expect_null(nmf_res$reassigned)
})

test_that("extract_nmf_info works with reference_P", {
  true_P = matrix(rpois(96*3, lambda = 19), nrow = 96)
  true_C = matrix(rpois(3*10, lambda = 1), nrow = 3)
  M = matrix(nrow = 96, ncol = 10)
  for (i in 1:10){
    M[,i] <- true_P %*% true_C[,i]
  }
  nmf_res <- nmf_wrapper(M, rank = 3, nrun = 1)
  nmf_res <- extract_nmf_info(nmf_res, reference_P = true_P)

  expect_equal(nrow(nmf_res$P), 96)
  expect_equal(ncol(nmf_res$P), 3)
  expect_equal(nrow(nmf_res$C), 3)
  expect_equal(ncol(nmf_res$C), 10)
  expect_true(nmf_res$minsim > 0.9)
  expect_equal(nrow(nmf_res$reassigned), 3)
  expect_equal(ncol(nmf_res$reassigned), 3)
  expect_equal(min(diag(nmf_res$reassigned)), nmf_res$minsim, tolerance = 0.001)
})
