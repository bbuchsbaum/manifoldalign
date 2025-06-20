test_that("cone_align_multiple recovers exact permutations on isomorphic graphs", {

  skip_on_cran()               # 150-200 ms locally, but leave CI time free

  set.seed(42)
  n  <- 80                     # nodes
  m  <- 4                      # graphs
  k  <- 8                      # embedding dim

  ## 1. base matrix with smooth structure ------------------------------
  X0 <- matrix(rnorm(n * 5), n)   # 5 "features" – only topology matters

  ## 2. make m permuted copies -----------------------------------------
  perms <- replicate(m, sample(n), simplify = FALSE)
  mats  <- lapply(perms, function(p) X0[p, , drop = FALSE])

  ## 3. run alignment ---------------------------------------------------
  res <- cone_align_multiple(mats, ncomp = k, max_iter = 40, tol = 1e-3)

  ## 4. assess accuracy -------------------------------------------------
  recovered <- res$assignment             # list of length m
  ref_inv <- match(seq_len(n), perms[[1]])
  for (g in seq_len(m)) {
    expect_equal(recovered[[g]][ref_inv], match(seq_len(n), perms[[g]]))
  }
})
