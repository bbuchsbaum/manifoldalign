test_that("cone_align_multiple runs in reasonable time and memory", {

  skip_on_cran()
  skip_if_not(Sys.info()[["sysname"]] != "Windows")  # peak RAM check uses 'ps'

  # The ps library is not strictly required, so skip if not available
  if (!requireNamespace("ps", quietly = TRUE)) {
    skip("ps package not available for performance test")
  }

  set.seed(7)
  n  <- 500                 # 500 nodes
  m  <- 4
  mats <- lapply(1:m, \(i) matrix(rnorm(n * 4), n, 4))

  p0   <- ps::ps_handle()
  # Defensive check: ensure p0 is a valid handle before trying to use it
  if (!inherits(p0, "ps_handle")) {
      skip("Could not get valid process handle for memory test")
  }
  
  mem0 <- ps::ps_memory_info(p0)$rss

  t0 <- proc.time()
  res <- cone_align_multiple(mats, ncomp = 4, max_iter = 15, tol = 1e-2)
  t1 <- proc.time() - t0

  mem_peak <- ps::ps_memory_info(p0)$rss_peak - mem0

  # Check: reasonable time and memory usage
  expect_lt(t1[["elapsed"]], 20)           # Should be faster than 20s
  expect_lt(mem_peak / 1e9, 2)             # Should use < 2 GB peak memory

  # Check: valid output
  expect_s3_class(res, "cone_align_multiple")
  expect_type(res$assignment, "list")
  expect_length(res$assignment, m)
}) 