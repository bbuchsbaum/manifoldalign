test_that("createSimFun preserves sparsity and symmetry", {
  S <- matrix(c(1, 0.2, 0.3,
                0.2, 1, 0.4,
                0.3, 0.4, 1), nrow = 3, byrow = TRUE,
              dimnames = list(letters[1:3], letters[1:3]))

  expect_warning(simfun_warn <- createSimFun(S, na_value = 5),
                 "ignores non-zero na_value", fixed = TRUE)
  expect_type(simfun_warn, "closure")

  simfun <- createSimFun(S)
  labels <- c("a", "b", NA, "c", "ghost")
  M <- simfun(labels)

  expect_s4_class(M, "dgCMatrix")
  expect_true(Matrix::isSymmetric(M))
  expect_equal(as.numeric(Matrix::diag(M)), rep(0, length(labels)))
  expect_gt(Matrix::nnzero(M[1:3, 1:3]), 0)
  expect_equal(Matrix::nnzero(M[5, ]), 0)
  expect_equal(Matrix::nnzero(M[, 5]), 0)
})

test_that("lowrank_align agrees across solver backends", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("tibble")

  set.seed(42)
  X1 <- matrix(rnorm(5 * 3), 5, 3)
  X2 <- matrix(rnorm(4 * 3), 4, 3)

  labels1 <- c("L1", "L2", NA, "L1", "L3")
  labels2 <- c("L1", NA, "L3", "L2")

  if (!("package:tibble" %in% search())) {
    suppressWarnings(library(tibble))
  }

  d1 <- multidesign::multidesign(X1, tibble::tibble(y = labels1))
  d2 <- multidesign::multidesign(X2, tibble::tibble(y = labels2))
  hd <- multidesign::hyperdesign(list(block1 = d1, block2 = d2))

  label_pool <- sort(unique(stats::na.omit(c(labels1, labels2))))
  S <- diag(length(label_pool))
  dimnames(S) <- list(label_pool, label_pool)
  simfun <- createSimFun(S)

  total_samples <- length(labels1) + length(labels2)
  total_features <- ncol(X1) + ncol(X2)

  res_exp <- lowrank_align(hd, y,
                           simfun = simfun,
                           ncomp = 2,
                           solver = "explicit",
                           sv_thresh = 0.5,
                           lambda = 0.01)

  res_op <- lowrank_align(hd, y,
                          simfun = simfun,
                          ncomp = 2,
                          solver = "operator",
                          sv_thresh = 0.5,
                          lambda = 0.01)

  expect_s3_class(res_exp, c("lowrank_align", "multiblock_biprojector"))
  expect_s3_class(res_op, c("lowrank_align", "multiblock_biprojector"))

  expect_equal(dim(res_exp$s), c(total_samples, 2))
  expect_equal(dim(res_op$s), c(total_samples, 2))
  expect_equal(nrow(res_exp$v), total_features)
  expect_equal(nrow(res_op$v), total_features)

  expect_equal(sum(is.na(res_exp$labels)), sum(is.na(c(labels1, labels2))))
  expect_equal(sum(is.na(res_op$labels)), sum(is.na(c(labels1, labels2))))

  overlap <- svd(crossprod(res_exp$s, res_op$s))
  expect_gt(min(overlap$d), 0.8)
})
