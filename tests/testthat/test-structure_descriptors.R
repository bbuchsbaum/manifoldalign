test_that("hks_time_grid fixed spacing behaves as expected", {
  vals <- c(0, 0.5, 2)

  t_lin <- hks_time_grid(vals, q = 4, mode = "fixed", time_range = c(0.1, 1.1),
                         spacing = "linear", time_scale = 2)
  expect_equal(t_lin, seq(0.1, 1.1, length.out = 4) * 2)

  t_log <- hks_time_grid(vals, q = 3, mode = "fixed", time_range = c(0.1, 100),
                         spacing = "log", time_scale = 10)
  expect_equal(t_log, exp(seq(log(0.1), log(100), length.out = 3)) * 10)
})

test_that("hks_time_grid auto derives a reasonable range from eigenvalues", {
  vals <- c(0, 0.5, 2)
  t <- hks_time_grid(vals, q = 2, mode = "auto", spacing = "log", time_scale = 1)
  expect_equal(t, c(0.5, 2), tolerance = 1e-10)
})

test_that("compute_hks_descriptors matches a simple identity basis case", {
  values <- c(0, 1, 2)
  vectors <- diag(3)
  times <- c(1, 2)

  desc <- compute_hks_descriptors(
    values = values,
    vectors = vectors,
    times = times,
    eigen_tol = 1e-12,
    normalize = "none"
  )

  expected <- rbind(
    c(0, 0),
    exp(-1 * times),
    exp(-2 * times)
  )
  expect_equal(desc, expected, tolerance = 1e-12)

  # Sparse eigenvectors should behave the same.
  vectors_sp <- Matrix::Matrix(vectors, sparse = TRUE)
  desc_sp <- compute_hks_descriptors(
    values = values,
    vectors = vectors_sp,
    times = times,
    eigen_tol = 1e-12,
    normalize = "none"
  )
  expect_equal(desc_sp, expected, tolerance = 1e-12)
})

test_that("compute_hks_descriptors normalization produces unit-norm rows/cols", {
  values <- c(0, 1, 2)
  vectors <- diag(3)
  times <- c(1, 2, 3)

  desc_both <- compute_hks_descriptors(
    values = values,
    vectors = vectors,
    times = times,
    eigen_tol = 1e-12,
    normalize = "both"
  )

  # Row 1 is identically zero; other rows should be unit norm.
  rn <- sqrt(rowSums(desc_both * desc_both))
  expect_equal(rn[2:3], c(1, 1), tolerance = 1e-12)

  # Column norms are not preserved after row-normalization; just sanity check.
  cn <- sqrt(colSums(desc_both * desc_both))
  expect_true(all(is.finite(cn)))

  desc_col <- compute_hks_descriptors(
    values = values,
    vectors = vectors,
    times = times,
    eigen_tol = 1e-12,
    normalize = "column"
  )
  cn2 <- sqrt(colSums(desc_col * desc_col))
  expect_true(all(abs(cn2 - 1) < 1e-12))
})

test_that("compute_diffusion_coordinates matches identity basis case", {
  values <- c(0, 1, 2)
  vectors <- diag(3)

  coords <- compute_diffusion_coordinates(
    values = values,
    vectors = vectors,
    time = 1,
    eigen_tol = 1e-12,
    normalize = "none"
  )

  expected <- cbind(
    c(0, exp(-1), 0),
    c(0, 0, exp(-2))
  )
  expect_equal(coords, expected, tolerance = 1e-12)

  coords_list <- compute_diffusion_coordinates(
    values = values,
    vectors = vectors,
    time = c(1, 2),
    eigen_tol = 1e-12
  )
  expect_true(is.list(coords_list))
  expect_equal(length(coords_list), 2L)
  expect_equal(coords_list[[1]], expected, tolerance = 1e-12)
})

test_that("adjacency-based helpers produce consistent shapes", {
  # Path graph on 5 nodes.
  n <- 5L
  ii <- c(1:4, 2:5)
  jj <- c(2:5, 1:4)
  A <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, n))

  L <- compute_graph_laplacian(A, normalized = TRUE)
  expect_true(inherits(L, "CsparseMatrix"))
  expect_equal(dim(L), c(n, n))

  basis <- compute_laplacian_basis(A, k = 2, normalized = TRUE)
  expect_equal(length(basis$values), 2L)
  expect_equal(dim(basis$vectors), c(n, 2L))
  expect_true(all(is.finite(basis$values)))
  expect_true(all(basis$values > 0))

  desc <- compute_hks_from_adjacency(A, k_embed = 2, q = 6, time_mode = "fixed",
                                    time_range = c(0.1, 2), spacing = "log",
                                    time_scale = 1, normalize = "both")
  expect_equal(dim(desc), c(n, 6L))
  rn <- sqrt(rowSums(desc * desc))
  expect_true(all(abs(rn - 1) < 1e-8))

  coords <- compute_diffusion_coordinates_from_adjacency(A, k_embed = 2, time = 1)
  expect_equal(dim(coords), c(n, 2L))
})
