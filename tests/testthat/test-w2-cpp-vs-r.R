# Independent R oracle: squared 1-D Wasserstein-2 distance between weighted
# empirical measures, as the integral over t in (0, 1) of
# (F^-1(t) - G^-1(t))^2 evaluated exactly on the merged breakpoints of the two
# cumulative weight functions.
w2_1d_oracle <- function(x, wx, y, wy) {
  wx <- wx / sum(wx)
  wy <- wy / sum(wy)
  ox <- order(x)
  oy <- order(y)
  x <- x[ox]; wx <- wx[ox]
  y <- y[oy]; wy <- wy[oy]
  cx <- cumsum(wx); cy <- cumsum(wy)
  br <- sort(unique(c(0, cx, cy, 1)))
  br <- br[br >= 0 & br <= 1]
  len <- diff(br)
  mid <- (head(br, -1) + tail(br, -1)) / 2
  qx <- x[pmin(findInterval(mid, cx) + 1L, length(x))]
  qy <- y[pmin(findInterval(mid, cy) + 1L, length(y))]
  sum(len * (qx - qy)^2)
}

w2_matrix_oracle <- function(Uref, Ug, wr, wg) {
  K <- ncol(Uref)
  C <- matrix(0, K, K)
  S <- matrix(0L, K, K)
  for (k in seq_len(K)) {
    for (l in seq_len(K)) {
      cpos <- w2_1d_oracle(Uref[, k], wr, Ug[, l], wg)
      cneg <- w2_1d_oracle(Uref[, k], wr, -Ug[, l], wg)
      if (cpos <= cneg) {
        C[k, l] <- cpos; S[k, l] <- 1L
      } else {
        C[k, l] <- cneg; S[k, l] <- -1L
      }
    }
  }
  list(C = C, sign = S)
}

test_that("w2_cost_matrix_1d returns valid costs and signs", {
  skip_if_missing_cd_deps()
  set.seed(123)
  n <- 50; K <- 6
  Uref <- matrix(rnorm(n * K), n, K)
  Ug   <- matrix(rnorm(n * K), n, K)
  wr <- runif(n); wg <- runif(n)

  # Call C++ kernel
  out_cpp <- w2_cost_matrix_1d(Uref, Ug, wr, wg)
  C_cpp <- as.matrix(out_cpp$C)
  S_cpp <- as.matrix(out_cpp$sign)
  expect_equal(dim(C_cpp), c(K, K))
  expect_equal(dim(S_cpp), c(K, K))
  expect_true(all(is.finite(C_cpp)))
  expect_true(all(C_cpp >= -1e-12))
  expect_true(all(S_cpp %in% c(-1L, 1L)))
  # Self-consistency: identical inputs should give near-zero diagonal and +1 signs
  out_id <- w2_cost_matrix_1d(Uref, Uref, wr, wr)
  C_id <- as.matrix(out_id$C); S_id <- as.matrix(out_id$sign)
  expect_true(max(diag(C_id)) < 1e-12)
  expect_true(all(diag(S_id) == 1L))
})

test_that("w2_cost_matrix_1d matches an independent quantile-function oracle", {
  set.seed(2024)
  for (case in 1:3) {
    n1 <- c(7, 31, 12)[case]
    n2 <- c(11, 17, 12)[case]
    K <- 4
    Uref <- matrix(rnorm(n1 * K), n1, K)
    Ug <- matrix(rnorm(n2 * K, sd = 2) + 0.5, n2, K)
    wr <- runif(n1)
    wg <- runif(n2)
    if (case == 3) {
      wr <- rep(1, n1)   # unnormalised uniform weights
      wg <- rep(3, n2)
    }
    out <- w2_cost_matrix_1d(Uref, Ug, wr, wg)
    ref <- w2_matrix_oracle(Uref, Ug, wr, wg)
    expect_equal(as.matrix(out$C), ref$C, tolerance = 1e-10)
    expect_equal(as.matrix(out$sign), ref$sign)
  }

  # Hand-checkable case: equal weights, two points each -> W2^2 = mean squared
  # difference of the sorted values.
  out <- w2_cost_matrix_1d(matrix(c(0, 1), 2, 1), matrix(c(3, 5), 2, 1),
                           c(1, 1), c(1, 1))
  expect_equal(out$C[1, 1], min(((3 - 0)^2 + (5 - 1)^2) / 2,
                                ((-5 - 0)^2 + (-3 - 1)^2) / 2))
})

test_that("w2_cost_matrix_1d rejects non-finite values and invalid weights", {
  U <- matrix(rnorm(10), 5, 2)
  w <- rep(1, 5)
  bad <- U
  bad[2, 1] <- NA
  expect_error(w2_cost_matrix_1d(bad, U, w, w), "non-finite")
  bad[2, 1] <- Inf
  expect_error(w2_cost_matrix_1d(U, bad, w, w), "non-finite")
  expect_error(w2_cost_matrix_1d(U, U, c(1, NA, 1, 1, 1), w), "non-finite")
  expect_error(w2_cost_matrix_1d(U, U, c(1, -1, 1, 1, 1), w), "non-negative")
  expect_error(w2_cost_matrix_1d(U, U, w, rep(0, 5)), "positive sum")
})
