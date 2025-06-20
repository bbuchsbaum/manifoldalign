test_that("linear_sim_embed recovers known latent subspace", {
  set.seed(1)
  n     <- 400        # moderate n so O(n²) still OK
  d_lat <- 2
  d_obs <- 10

  # latent X*: two concentric Gaussian rings
  r     <- sample(c(1, 5), n, replace = TRUE) + rnorm(n, 0, .1)
  theta <- runif(n, 0, 2*pi)
  Xlat  <- cbind(r*cos(theta), r*sin(theta))

  # random orthonormal projector  P : R^2 -> R^10
  Q     <- qr.Q(qr(matrix(rnorm(d_obs*d_lat), d_obs, d_lat)))
  Xobs  <- Xlat %*% t(Q) + matrix(rnorm(n*d_obs, 0, .05), n)   # add small noise

  # target similarity from *latent* distances – gold standard
  D2    <- as.matrix(dist(Xlat))^2
  Tstar <- exp(-D2 / median(D2))

  fit   <- linear_sim_embed(Xobs, T = Tstar, ncomp = 2,
                            sigma_P = 1, alpha_p = 0.05,
                            maxit = 800, tol = 1e-7, verbose = FALSE)

  # Align to ground truth with Procrustes
  Yhat  <- fit$scores
  Yhat  <- scale(Yhat, center = TRUE, scale = FALSE)
  proc  <- svd(t(Xlat) %*% Yhat)
  Yrot  <- Yhat %*% proc$v %*% t(proc$u)

  err   <- sqrt(mean((Xlat - Yrot)^2))          # RMSE in latent space
  pcaRMSE <- sqrt(mean((Xlat - prcomp(Xobs, rank.=2)$x)^2))

  expect_lt(err, 1.05 * pcaRMSE)               # allow 5 % slack
}) 