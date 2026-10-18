test_that("diagnostics: orthogonal sync residuals are near zero on noiseless data", {
  skip_on_cran()
  N <- 3; k <- 2
  rot <- function(theta) {
    c <- cos(theta); s <- sin(theta)
    matrix(c(c,-s,s,c), 2, 2)
  }
  R1 <- diag(k)
  R2 <- rot(pi/6)     # 30 degrees
  R3 <- rot(-pi/10)   # -18 degrees
  T_list <- list(
    manifoldalign::new_align_transform("O", R1, from = 1, to = "consensus", k = k),
    manifoldalign::new_align_transform("O", R2, from = 2, to = "consensus", k = k),
    manifoldalign::new_align_transform("O", R3, from = 3, to = "consensus", k = k)
  )
  E <- data.frame(i = c(1,1,2), j = c(2,3,3))
  # Row-action contract (?new_align_transform): X_i g_ij = X_j and
  # X_i T_i = consensus, so consistency is T_i = g_ij T_j, g_ij = R_i R_j'.
  G <- list(
    manifoldalign::new_align_transform("O", R1 %*% t(R2), from = 1, to = 2, k = k),
    manifoldalign::new_align_transform("O", R1 %*% t(R3), from = 1, to = 3, k = k),
    manifoldalign::new_align_transform("O", R2 %*% t(R3), from = 2, to = 3, k = k)
  )
  W <- rep(1, nrow(E))
  diags <- manifoldalign:::compute_diagnostics_("O", T_list, G, E, W)
  expect_true(is.numeric(diags$mean_edge_residual))
  expect_lt(diags$mean_edge_residual, 1e-10)
  expect_false(is.null(diags$cycle_residuals))
  expect_lt(diags$mean_cycle_residual, 1e-10)

  # The certified consensus really aligns the observations.
  set.seed(1)
  X <- matrix(rnorm(20), 10, 2)
  X2 <- X %*% G[[1]]$op
  X3 <- X %*% G[[2]]$op
  Z1 <- manifoldalign::apply_transform(T_list[[1]], X)
  expect_equal(manifoldalign::apply_transform(T_list[[2]], X2), Z1, tolerance = 1e-12)
  expect_equal(manifoldalign::apply_transform(T_list[[3]], X3), Z1, tolerance = 1e-12)
})

test_that("diagnostics: permutation sync residuals and entropies on perfect maps", {
  skip_on_cran()
  N <- 3; L <- 3; n <- c(3,3,3)
  # Define per-domain permutations
  p1 <- c(1,2,3)
  p2 <- c(2,3,1)
  p3 <- c(3,1,2)
  # Build T_list: each is L x n_i with one-hot per column indicating consensus row
  make_T <- function(p) {
    M <- matrix(0, L, length(p)); M[cbind(p, seq_along(p))] <- 1; M
  }
  T1 <- make_T(p1); T2 <- make_T(p2); T3 <- make_T(p3)
  T_list <- list(
    manifoldalign::new_align_transform("perm", T1, from = 1, to = "consensus"),
    manifoldalign::new_align_transform("perm", T2, from = 2, to = "consensus"),
    manifoldalign::new_align_transform("perm", T3, from = 3, to = "consensus")
  )
  # Build relative maps G_ij = t(T_j) %*% T_i  (n_j x n_i)
  G12 <- t(T2) %*% T1; G13 <- t(T3) %*% T1; G23 <- t(T3) %*% T2
  G <- list(
    manifoldalign::new_align_transform("perm", G12, from = 1, to = 2),
    manifoldalign::new_align_transform("perm", G13, from = 1, to = 3),
    manifoldalign::new_align_transform("perm", G23, from = 2, to = 3)
  )
  E <- data.frame(i = c(1,1,2), j = c(2,3,3))
  W <- rep(1, nrow(E))
  diags <- manifoldalign:::compute_diagnostics_("perm", T_list, G, E, W)
  expect_true(is.numeric(diags$mean_edge_residual))
  expect_lt(diags$mean_edge_residual, 1e-10)
  if (!is.null(diags$cycle_residuals)) {
    expect_lt(diags$mean_cycle_residual, 1e-10)
  }
  # Entropy should be ~0 for perfect permutations
  expect_true(is.data.frame(diags$entropy))
  expect_true(all(diags$entropy$mean_col_entropy <= 1e-12))
  expect_true(all(diags$entropy$mean_row_entropy <= 1e-12))

  # Row action: X_j = G_ij X_i and every domain lands on the same consensus.
  X1 <- matrix(c(1, 2, 3, 10, 20, 30), 3, 2)
  Z1 <- manifoldalign::apply_transform(T_list[[1]], X1)
  expect_equal(manifoldalign::apply_transform(T_list[[2]], G12 %*% X1), Z1)
  expect_equal(manifoldalign::apply_transform(T_list[[3]], G13 %*% X1), Z1)
})

