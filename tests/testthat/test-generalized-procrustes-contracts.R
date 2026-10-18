# Contracts of generalized_procrustes(): preprocessing lists, translations for
# partial task designs, final_diff semantics and the optimality certificate.

gpa_rand_O <- function(d) qr.Q(qr(matrix(rnorm(d * d), d)))

gpa_global_rotation <- function(O_hat, O_true) {
  M <- Reduce(`+`, Map(function(Oh, Ot) t(Ot) %*% Oh, O_hat, O_true))
  sv <- svd(M)
  sv$u %*% t(sv$v)
}

# Hyperdesign with rows = tasks: X_i[j, ] = (R_i a_j + t_i)^T
gpa_hyperdesign <- function(template, R_list, shifts, task_sets, task_names) {
  doms <- lapply(seq_along(R_list), function(i) {
    A_i <- R_list[[i]] %*% template[, task_sets[[i]], drop = FALSE] + shifts[[i]]
    multidesign::multidesign(t(A_i), data.frame(task = task_names[task_sets[[i]]]))
  })
  names(doms) <- paste0("d", seq_along(doms))
  multidesign::hyperdesign(doms)
}

test_that("generalized_procrustes accepts preproc = NULL and NULL list entries", {
  skip_if_not_installed("multidesign")
  set.seed(21)
  d <- 3; L <- 6
  template <- matrix(rnorm(d * L), d, L)
  template <- template - rowMeans(template)
  R_list <- replicate(3, gpa_rand_O(d), simplify = FALSE)
  shifts <- replicate(3, rep(0, d), simplify = FALSE)
  tasks <- replicate(3, seq_len(L), simplify = FALSE)
  hd <- gpa_hyperdesign(template, R_list, shifts, tasks, paste0("t", seq_len(L)))

  # Regression: `proclist[[i]] <- NULL` shrank the list (two code paths).
  set.seed(1)
  gp_null <- generalized_procrustes(hd, task, preproc = NULL, max_iter = 200, tol = 1e-10)
  set.seed(1)
  gp_mixed <- generalized_procrustes(hd, task, preproc = list(multivarious::center(), NULL, NULL),
                                     max_iter = 200, tol = 1e-10)
  set.seed(1)
  gp_center <- generalized_procrustes(hd, task, max_iter = 200, tol = 1e-10)

  expect_length(gp_null$preproc, 3L)
  expect_true(all(vapply(gp_null$preproc, is.null, logical(1))))
  expect_named(gp_null$preproc, names(hd))
  expect_length(gp_mixed$preproc, 3L)
  expect_false(is.null(gp_mixed$preproc[[1]]))
  expect_null(gp_mixed$preproc[[2]])

  # The data are already centred, so every variant recovers the same consensus.
  for (gp in list(gp_null, gp_mixed, gp_center)) {
    Q <- gpa_global_rotation(gp$O_mats, R_list)
    expect_equal(Q %*% gp$A_est, template, tolerance = 1e-6, ignore_attr = TRUE)
  }
})

test_that("partial task designs are centred on shared tasks (no translation bias)", {
  skip_if_not_installed("multidesign")
  set.seed(5)
  d <- 3; L <- 8
  template <- matrix(rnorm(d * L, sd = 2), d, L)
  R_list <- replicate(4, gpa_rand_O(d), simplify = FALSE)
  shifts <- replicate(4, rnorm(d, sd = 3), simplify = FALSE)
  # Tasks 1:4 are observed by every domain; the rest are partial.
  task_sets <- list(c(1:4, 5, 6), c(1:4, 7), c(1:4, 5, 8), c(1:4, 6, 7, 8))
  hd <- gpa_hyperdesign(template, R_list, shifts, task_sets, paste0("t", seq_len(L)))

  set.seed(1)
  gp <- generalized_procrustes(hd, task, max_iter = 300, tol = 1e-10)
  expect_true(gp$converged)
  # The rotations are a certified global optimum of the Procrustes objective.
  expect_true(gp$certificate$certified)

  # Oracle: the template expressed relative to the mean of the shared tasks.
  target <- template - rowMeans(template[, 1:4])
  ord <- match(paste0("t", seq_len(L)), gp$task_labels)
  Q <- gpa_global_rotation(gp$O_mats, R_list)
  expect_equal(Q %*% gp$A_est[, ord], target, tolerance = 1e-6, ignore_attr = TRUE)

  # Translations: the shared-task mean of each preprocessed domain.
  expect_named(gp$translations, names(hd))
  for (i in seq_along(R_list)) {
    X_i <- t(R_list[[i]] %*% template[, task_sets[[i]]] + shifts[[i]])
    X_i <- sweep(X_i, 2, colMeans(X_i))
    expect_equal(unname(gp$translations[[i]]), unname(colMeans(X_i[1:4, ])), tolerance = 1e-10)
  }

  # Per-domain centring over different task subsets (translation = "none")
  # leaves domain-specific offsets that bias the consensus.
  set.seed(1)
  gp_none <- generalized_procrustes(hd, task, max_iter = 300, tol = 1e-10,
                                    translation = "none")
  Q_none <- gpa_global_rotation(gp_none$O_mats, R_list)
  err_none <- max(abs(Q_none %*% gp_none$A_est[, ord] -
                        (template - rowMeans(template))))
  expect_gt(err_none, 0.05)
  expect_null(gp_none$translations)
})

test_that("a design with no task shared by all domains warns", {
  skip_if_not_installed("multidesign")
  set.seed(8)
  d <- 2; L <- 4
  template <- matrix(rnorm(d * L), d, L)
  R_list <- replicate(3, gpa_rand_O(d), simplify = FALSE)
  shifts <- replicate(3, rep(0, d), simplify = FALSE)
  task_sets <- list(c(1, 2, 3), c(2, 3, 4), c(1, 4))
  hd <- gpa_hyperdesign(template, R_list, shifts, task_sets, paste0("t", seq_len(L)))
  expect_warning(
    gp <- generalized_procrustes(hd, task, max_iter = 50),
    "No task is observed by every domain"
  )
  expect_null(gp$translations)
})

test_that("final_diff is the change at the last iteration", {
  set.seed(3)
  d <- 3; L <- 5
  template <- matrix(rnorm(d * L), d, L)
  A_list <- lapply(seq_len(3), function(i) gpa_rand_O(d) %*% template +
                     0.3 * matrix(rnorm(d * L), d, L))
  labels <- replicate(3, seq_len(L), simplify = FALSE)

  set.seed(10)
  gp4 <- generalized_procrustes(A_list, labels, L, max_iter = 4, tol = 0)
  set.seed(10)
  gp3 <- generalized_procrustes(A_list, labels, L, max_iter = 3, tol = 0)
  expect_false(gp4$converged)
  step4 <- sqrt(sum(mapply(function(a, b) sum((a - b)^2), gp4$O_mats, gp3$O_mats)))
  expect_equal(gp4$final_diff, step4, tolerance = 1e-10)
})

test_that("certificate detects non-stationary rotations", {
  set.seed(4)
  d <- 3; L <- 5; n <- 3
  template <- matrix(rnorm(d * L), d, L)
  O_true <- replicate(n, gpa_rand_O(d), simplify = FALSE)
  A_list <- lapply(O_true, function(O) O %*% template)
  D <- do.call(rbind, A_list)

  O_arr <- array(unlist(O_true), c(d, d, n))
  cert_true <- manifoldalign:::.gpa_certificate(D, O_arr)
  expect_identical(cert_true$status, "computed")
  expect_true(cert_true$certified)

  O_bad <- O_arr
  O_bad[, , 2] <- gpa_rand_O(d)
  cert_bad <- manifoldalign:::.gpa_certificate(D, O_bad)
  expect_gt(cert_bad$stationarity_residual, 1e-3)
  expect_false(cert_bad$certified)

  # Large problems report stationarity only.
  cert_big <- manifoldalign:::.gpa_certificate(D, O_arr, max_dense = 2L)
  expect_identical(cert_big$status, "stationarity_only")
  expect_true(is.na(cert_big$certified))
})
