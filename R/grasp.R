#' Graph Alignment by Spectral Corresponding Functions (GRASP) - Corrected Version
#'
#' Performs GRASP alignment on hyperdesign data structures. 
#' Projects graph data from multiple domains into aligned spectral spaces.
#'
#' @param data A hyperdesign object containing multiple graph domains
#' @param preproc Preprocessing function to apply to the data (default: center)
#' @param ncomp Number of spectral components to extract (default: 30)
#' @param q_descriptors Number of diffusion-time descriptors (default: 100)
#' @param sigma Diffusion parameter for descriptor computation (default: 0.73)
#' @param lambda Regularization parameter for base alignment (default: 0.1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param solver Assignment method: "linear" for exact assignment (default)
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the GRASP alignment
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign graph data
#' library(multidesign)
#' library(tibble)
#' 
#' # Create synthetic graph domains (node coordinates for graph construction)
#' set.seed(123)
#' X1 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 1
#' X2 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 2
#' 
#' # Create design data frames with node indices (GRASP doesn't use features)
#' design1 <- tibble(node_id = 1:50)
#' design2 <- tibble(node_id = 1:50)
#' 
#' # Create multidesign objects
#' md1 <- multidesign(X1, design1)
#' md2 <- multidesign(X2, design2)
#' 
#' # Create hyperdesign from multidesign objects
#' hd <- hyperdesign(list(domain1 = md1, domain2 = md2))
#' 
#' # Run GRASP alignment (purely spectral, uses graph structure from X)
#' result <- grasp(hd, ncomp = 20, q_descriptors = 50)
#' }
#'
#' @rdname grasp
#' @method grasp hyperdesign
#' @export
#' @importFrom multivarious init_transform prep concat_pre_processors
#' @importFrom manifoldalign block_indices
#' @importFrom stats dist
grasp.hyperdesign <- function(data, 
                             preproc = multivarious::center(), 
                             ncomp = 30,
                             q_descriptors = 100,
                             sigma = 0.73,
                             lambda = 0.1,
                             use_laplacian = TRUE,
                             solver = "linear",
                             ...) {
  
  # Input validation (following package patterns + your suggestions)
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(q_descriptors)
  chk::chk_true(q_descriptors > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  chk::chk_number(lambda)
  chk::chk_true(lambda > 0)
  chk::chk_logical(use_laplacian)
  
  # Validate solver parameter
  if (!solver %in% c("linear", "hungarian", "auction")) {
    stop("solver must be 'linear', 'hungarian', or 'auction'", call. = FALSE)
  }
  
  # Validate input data - should be a hyperdesign object containing 2 domains
  chk::chk_s3_class(data, "hyperdesign")
  
  if (!is.list(data) || length(data) == 0) {
    stop("hyperdesign data must contain domains", call. = FALSE)
  }
  
  if (length(data) != 2) {
    stop("GRASP currently supports exactly 2 domains. For multi-graph ", 
         "alignment with 3+ domains, use grasp_multiset()", call. = FALSE)
  }
  
  # Apply preprocessing per domain.
  # Avoid deprecated `prep()` / `init_transform()` paths in multivarious >= 0.3.0.
  pdata <- unclass(data)
  nb <- length(pdata)

  if (is.null(preproc)) {
    proclist <- vector("list", nb)
    for (i in seq_len(nb)) {
      pdata[[i]]$x <- as.matrix(pdata[[i]]$x)
      proclist[[i]] <- NULL
    }
  } else {
    is_preproc_obj <- inherits(preproc, "prepper") || inherits(preproc, "pre_processor")

    if (!is.list(preproc) || is_preproc_obj) {
      preproc <- rep(list(preproc), nb)
    } else if (length(preproc) == 1L) {
      preproc <- rep(preproc, nb)
    } else if (length(preproc) != nb) {
      stop(
        "Length of `preproc` list (",
        length(preproc),
        ") must match number of domains (",
        nb,
        ").",
        call. = FALSE
      )
    }

    proclist <- vector("list", nb)
    for (i in seq_len(nb)) {
      Xi <- as.matrix(pdata[[i]]$x)
      pre_i <- preproc[[i]]

      if (is.null(pre_i)) {
        pdata[[i]]$x <- Xi
        proclist[[i]] <- NULL
        next
      }

      # Functions are applied directly (avoid recipes-based preprocessing).
      if (is.function(pre_i)) {
        Xi_tf <- pre_i(Xi)
        pdata[[i]]$x <- Xi_tf
        proclist[[i]] <- pre_i
        next
      }

      # Clone stateful preprocessors before fitting.
      if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
        pre_i <- unserialize(serialize(pre_i, connection = NULL))
      }

      if (exists("fit_transform", envir = asNamespace("multivarious"), mode = "function")) {
        ft <- multivarious::fit_transform(pre_i, Xi)
        pdata[[i]]$x <- ft$transformed
        proclist[[i]] <- ft$preproc
      } else {
        # Back-compat for older multivarious versions
        templ <- multivarious::prep(pre_i)
        Xi_tf <- multivarious::init_transform(templ, Xi)
        pdata[[i]]$x <- Xi_tf
        proc_attr <- attr(Xi_tf, "preproc")
        proclist[[i]] <- if (is.null(proc_attr)) templ else proc_attr
      }
    }
  }

  names(pdata) <- names(data)
  class(pdata) <- class(data)
  names(proclist) <- names(pdata)

  sample_block_idx <- block_indices(pdata)
  sample_blocks_list <- split(sample_block_idx, row(sample_block_idx))
  names(sample_blocks_list) <- names(pdata)

  # Concatenate preprocessors (following package pattern from coupled_diagonalization.R)
  proc <- multivarious::concat_pre_processors(proclist, sample_blocks_list)

  grasp_fit(pdata, proc, ncomp, q_descriptors, sigma, lambda,
            use_laplacian, solver, feature_block_indices(pdata))
}

#' @keywords internal
grasp_fit <- function(strata, proc, ncomp, q_descriptors, sigma, lambda, 
                     use_laplacian, solver, feature_blocks) {
  
  # Block A: Spectral basis construction
  bases <- compute_grasp_basis(strata, ncomp, use_laplacian)
  
  # Block B: Multi-scale descriptors  
  descriptors <- compute_grasp_descriptors(bases, q_descriptors, sigma)
  
  # Block C: Base alignment
  alignment_result <- align_grasp_bases(bases[[1]], bases[[2]], 
                                       descriptors[[1]], descriptors[[2]], lambda)
  
  # Blocks D & E: Functional map and assignment
  assignment_result <- compute_grasp_assignment(bases[[1]], bases[[2]], 
                                               descriptors[[1]], descriptors[[2]],
                                               alignment_result$rotation,
                                               distance_method = "cosine",
                                               solver_method = solver)
  
  # Compute scores (following package pattern)
  embed1 <- as.matrix(bases[[1]]$vectors)
  embed2 <- as.matrix(bases[[2]]$vectors %*% alignment_result$rotation %*% assignment_result$mapping_matrix)
  perm <- as.integer(assignment_result$assignment)
  embed2_aligned <- embed2
  if (length(perm) == nrow(embed1)) {
    # assignment maps rows in domain 1 to rows in domain 2, so reorder directly.
    embed2_aligned <- matrix(0, nrow = nrow(embed1), ncol = ncol(embed2))
    valid <- perm > 0L & perm <= nrow(embed2)
    if (any(valid)) {
      embed2_aligned[valid, ] <- embed2[perm[valid], , drop = FALSE]
    }
  } else if (length(perm) == nrow(embed2)) {
    valid <- perm > 0L & perm <= nrow(embed2)
    if (all(valid)) {
      embed2_aligned <- embed2[perm, , drop = FALSE]
    }
  }
  scores <- rbind(embed1, embed2_aligned)
  
  # Primal vectors computation (following KEMA pattern)
  loadings <- do.call(rbind, lapply(seq_along(strata), function(i) {
    xi <- strata[[i]]$x
    alpha_i <- bases[[i]]$vectors
    Matrix::crossprod(xi, alpha_i)
  }))

  new_alignment_result(
    scores = scores,
    loadings = loadings,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "grasp",
    extras = list(
      assignment = assignment_result$assignment,
      rotation = alignment_result$rotation,
      mapping_matrix = assignment_result$mapping_matrix
    )
  )
}

#' Block A: Spectral Basis Construction - Corrected
#' 
#' @param strata List of data domains 
#' @param ncomp Number of spectral components
#' @param use_laplacian Whether to use Laplacian normalization
#' @return List of spectral bases for each domain
#' @keywords internal
compute_grasp_basis <- function(strata, ncomp, use_laplacian = TRUE) {
  # Validate input
  if (!is.list(strata) || length(strata) == 0) {
    stop("strata must be a non-empty list", call. = FALSE)
  }
  
  bases <- lapply(strata, function(stratum) {
    # Ensure stratum$x is a matrix for neighborweights
    if (!is.matrix(stratum$x)) {
      stratum$x <- as.matrix(stratum$x)
      if (!is.matrix(stratum$x)) {
        stop("Failed to convert domain data to matrix format", call. = FALSE)
      }
    }
    
    # Validate dimensions
    if (nrow(stratum$x) < 3) {
      stop("GRASP requires at least 3 samples per domain", call. = FALSE)
    }
    
    # Construct graph using package patterns (Patch 2)
    # Auto-tune sigma for graph construction via median distance heuristic
    sigma_g <- tryCatch(choose_sigma(stratum$x), error = function(e) 0.5)
    graph_weights <- tryCatch(
      neighborweights::graph_weights(
        stratum$x,
        weight_mode = "normalized",
        neighbor_mode = "knn",
        k = min(max(5, nrow(stratum$x) %/% 3), nrow(stratum$x) - 1),
        type = "normal",
        sigma = sigma_g
      ),
      error = function(e) {
        stop("Graph construction failed: ", e$message, call. = FALSE)
      }
    )

    adj_matrix <- neighborweights::adjacency(graph_weights)
    adj_matrix <- (adj_matrix + Matrix::t(adj_matrix)) / 2

    if (use_laplacian) {
      degrees <- Matrix::rowSums(adj_matrix)
      isolated <- as.vector(degrees == 0)
      if (any(isolated)) {
        warning("Found ", sum(isolated), " isolated nodes; Laplacian rows will be zeroed.")
      }
      deg_safe <- pmax(degrees, 1)
      D_inv_sqrt <- Matrix::Diagonal(x = 1 / sqrt(deg_safe))
      if (any(isolated)) {
        idx <- which(isolated)
        D_inv_sqrt[idx, idx] <- 0
      }
      L <- Matrix::Diagonal(nrow(adj_matrix)) - D_inv_sqrt %*% adj_matrix %*% D_inv_sqrt
      L <- Matrix::forceSymmetric(L)
    } else {
      L <- Matrix::forceSymmetric(adj_matrix)
    }

    total_dim <- nrow(L)
    if (total_dim <= 1) {
      stop("Graph must contain at least two nodes", call. = FALSE)
    }

    k_to_request <- if (use_laplacian) {
      min(ncomp + 1L, total_dim - 1L)
    } else {
      min(ncomp, total_dim)
    }

    if (k_to_request < 1L) {
      stop("Not enough non-trivial eigenvectors available for requested components", call. = FALSE)
    }

    decomp <- tryCatch({
      if (requireNamespace("RSpectra", quietly = TRUE) && k_to_request <= ceiling(0.8 * total_dim)) {
        RSpectra::eigs_sym(L, k = k_to_request, which = "SA")
      } else if (requireNamespace("PRIMME", quietly = TRUE) && k_to_request < total_dim) {
        PRIMME::eigs_sym(L, NEig = k_to_request, which = "SA")
      } else {
        eig_full <- eigen(as.matrix(L), symmetric = TRUE)
        list(
          values = eig_full$values[seq_len(k_to_request)],
          vectors = eig_full$vectors[, seq_len(k_to_request), drop = FALSE]
        )
      }
    }, error = function(e) {
      eig_full <- eigen(as.matrix(L), symmetric = TRUE)
      list(
        values = eig_full$values[seq_len(k_to_request)],
        vectors = eig_full$vectors[, seq_len(k_to_request), drop = FALSE]
      )
    })

    if (use_laplacian) {
      if (length(decomp$values) < 2L) {
        stop("Unable to extract non-trivial Laplacian eigenvectors", call. = FALSE)
      }
      idx <- 2:min(ncomp + 1L, length(decomp$values))
      if (length(idx) == 0L) {
        stop("No non-trivial Laplacian eigenvectors remained after skipping the constant mode", call. = FALSE)
      }
      vectors <- Matrix::Matrix(decomp$vectors[, idx, drop = FALSE], sparse = TRUE)
      values <- decomp$values[idx]
    } else {
      non_trivial_idx <- which(abs(decomp$values) > 1e-10)
      if (length(non_trivial_idx) == 0L) {
        stop("No non-trivial adjacency eigenvalues available", call. = FALSE)
      }
      ncomp_actual <- min(ncomp, length(non_trivial_idx))
      idx <- non_trivial_idx[seq_len(ncomp_actual)]
      vectors <- Matrix::Matrix(decomp$vectors[, idx, drop = FALSE], sparse = TRUE)
      values <- decomp$values[idx]
    }

    col_norms <- sqrt(Matrix::colSums(vectors * vectors))
    vectors <- vectors %*% Matrix::Diagonal(x = 1 / pmax(col_norms, 1e-12))

    list(vectors = vectors, values = values)
  })
  
  bases
}

#' Block B: Multi-scale Descriptors - Corrected
#' 
#' @param bases List of spectral bases
#' @param q_descriptors Number of descriptor functions
#' @param sigma Diffusion parameter
#' @return List of descriptor matrices
#' @keywords internal
compute_grasp_descriptors <- function(bases, q_descriptors, sigma = 0.73) {
  # Validation
  chk::chk_number(q_descriptors)
  chk::chk_true(q_descriptors > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)

  descriptors <- lapply(bases, function(basis) {
    desc_matrix <- compute_hks_descriptors(
      basis = basis,
      q = as.integer(q_descriptors),
      time_mode = "fixed",
      time_range = c(0.1, 50),
      spacing = "linear",
      time_scale = sigma,
      eigen_tol = 1e-12,
      normalize = "column",
      clamp_exponent = -745
    )
    Matrix::Matrix(desc_matrix, sparse = TRUE)
  })

  descriptors
}

#' Block C: Base Alignment - Corrected
#' 
#' @param basis1 First domain's spectral basis
#' @param basis2 Second domain's spectral basis  
#' @param desc1 First domain's descriptors
#' @param desc2 Second domain's descriptors
#' @param lambda Regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return List with rotation matrix and convergence info
#' @keywords internal
align_grasp_bases <- function(basis1, basis2, desc1, desc2, lambda = 0.1, 
                             max_iter = 200, tol = 1e-8) {
  # Validation
  chk::chk_number(lambda)
  chk::chk_true(lambda > 0)
  chk::chk_number(max_iter)
  chk::chk_number(tol)
  
  # Extract matrices
  phi1 <- basis1$vectors
  phi2 <- basis2$vectors
  lambda2_vals <- basis2$values
  
  k <- ncol(phi1)
  M <- Matrix::Diagonal(k)  # Initialize as identity
  
  # Pre-compute constant parts
  F_T_Phi1 <- Matrix::crossprod(desc1, phi1)
  G_T_Phi2 <- Matrix::crossprod(desc2, phi2)
  Lambda2 <- Matrix::Diagonal(n = length(lambda2_vals), x = lambda2_vals)
  
  # Corrected objective function (your suggestion - remove λ₁ subtraction)
  objective_fn <- function(M) {
    aligned_lambda2_matrix <- Matrix::crossprod(M, Lambda2 %*% M)
    # off() is sum of squares of off-diagonal elements (independent of λ₁)
    diag_vals <- Matrix::diag(aligned_lambda2_matrix)
    off_mat <- aligned_lambda2_matrix - Matrix::Diagonal(n = length(diag_vals), x = diag_vals)
    off_penalty <- Matrix::norm(off_mat, "F")^2
    
    desc_penalty <- Matrix::norm(F_T_Phi1 - G_T_Phi2 %*% M, "F")^2
    return(off_penalty + lambda * desc_penalty)
  }

  # Corrected gradient function
  gradient_fn <- function(M) {
    A <- Lambda2 %*% M
    S <- Matrix::crossprod(M, A)
    diag_S <- Matrix::diag(S)
    grad_off <- 4 * (A %*% S - A %*% Matrix::Diagonal(x = diag_S))

    residual <- F_T_Phi1 - G_T_Phi2 %*% M
    grad_desc <- -2 * Matrix::crossprod(G_T_Phi2, residual)

    euc_grad <- grad_off + lambda * grad_desc

    MtE <- Matrix::crossprod(M, euc_grad)
    riemannian_grad <- euc_grad - M %*% (0.5 * (MtE + Matrix::t(MtE)))
    riemannian_grad
  }

  # Improved optimization with relative tolerance (your suggestion)
  prev_obj <- objective_fn(M)
  
  for (iter in seq_len(max_iter)) {
    grad <- gradient_fn(M)
    
    # Backtracking line search with smaller initial step
    step_size <- 0.5
    for (ls_iter in 1:8) {
      M_new_euc <- M - step_size * grad
      
      # SVD retraction (your suggestion - more stable than QR for small k)
      M_new_mat <- as.matrix(M_new_euc)
      if (any(!is.finite(M_new_mat))) {
        stop("Non-finite values in matrix during SVD retraction", call. = FALSE)
      }
      svd_res <- svd(M_new_mat)
      M_new <- as(svd_res$u %*% t(svd_res$v), "CsparseMatrix")
      
      new_obj <- objective_fn(M_new)
      if (new_obj < prev_obj) break
      step_size <- step_size * 0.5
    }
    
    M <- M_new
    
    # Improved convergence check - require minimum iterations
    rel_change <- abs(prev_obj - new_obj) / max(1, prev_obj)
    if (iter > 5 && rel_change < tol) {  # Minimum 5 iterations
      break
    }
    prev_obj <- new_obj
  }
  
  list(rotation = M, iterations = iter, converged = iter < max_iter)
}

#' Blocks D & E: Functional Map and Assignment - Corrected
#' 
#' @param basis1 First domain's spectral basis
#' @param basis2 Second domain's spectral basis
#' @param desc1 First domain's descriptors  
#' @param desc2 Second domain's descriptors
#' @param M Rotation matrix from base alignment
#' @param distance_method Distance method: "cosine" or "euclidean"
#' @param solver_method Assignment solver: "linear", "hungarian", or "auction"
#' @return List with assignment and related matrices
#' @keywords internal
compute_grasp_assignment <- function(basis1, basis2, desc1, desc2, M, 
                                   distance_method = "cosine", 
                                   solver_method = "linear",
                                   alpha = 0.85) {
  # Validation
  if (!distance_method %in% c("cosine", "euclidean")) {
    stop("distance_method must be 'cosine' or 'euclidean'", call. = FALSE)
  }
  
  if (!solver_method %in% c("linear", "hungarian", "auction")) {
    stop("solver_method must be 'linear', 'hungarian', or 'auction'", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a finite scalar in [0, 1]", call. = FALSE)
  }
  
  # Block D: Functional correspondence with ridge regularization (your suggestion)
  phi1 <- basis1$vectors
  phi2_aligned <- basis2$vectors %*% M
  
  # Solve F^T*Phi1 ≈ G^T*Phi2_aligned * C for diagonal C
  A_prime <- Matrix::crossprod(desc2, phi2_aligned)  # (q, k)
  B_prime <- Matrix::crossprod(desc1, phi1)         # (q, k)
  
  # Ridge regularization for numerical stability (your suggestion)
  c_diag <- Matrix::colSums(A_prime * B_prime) / (Matrix::colSums(A_prime^2) + 1e-9)
  C <- Matrix::Diagonal(n = length(c_diag), x = c_diag)
  
  # Block E: Final embeddings and assignment
  embed1 <- as.matrix(phi1)
  embed2 <- as.matrix(phi2_aligned %*% C)

  # Augment with descriptor-driven embeddings for robustness (mirrors multiset features)
  B_proj <- as.matrix(desc1 %*% B_prime)           # n1 x k
  A_proj <- as.matrix(desc2 %*% (A_prime %*% C))   # n2 x k
  # Balance scales to avoid one dominating
  fn1 <- sqrt(sum(embed1^2)); fnB <- sqrt(sum(B_proj^2)); s1 <- if (fnB > 0) fn1 / fnB else 1
  fn2 <- sqrt(sum(embed2^2)); fnA <- sqrt(sum(A_proj^2)); s2 <- if (fnA > 0) fn2 / fnA else 1
  embed1 <- cbind(embed1, B_proj * s1)
  embed2 <- cbind(embed2, A_proj * s2)
  
  # Improved distance computation (cosine + small Euclidean blend)
  if (distance_method == "cosine") {
    # Robust cosine computation without external dependencies
    embed1_norms <- sqrt(rowSums(embed1^2) + 1e-12)
    embed2_norms <- sqrt(rowSums(embed2^2) + 1e-12)
    embed1_norm <- embed1 / embed1_norms
    embed2_norm <- embed2 / embed2_norms

    sim_matrix <- embed1_norm %*% t(embed2_norm)
    cost_cos <- 1 - sim_matrix
    cost_cos[cost_cos < 0] <- 0

    # Small Euclidean component to break cosine ties and stabilize
    e1_sq <- rowSums(embed1^2)
    e2_sq <- rowSums(embed2^2)
    dm <- sqrt(pmax(outer(e1_sq, e2_sq, "+") - 2 * (embed1 %*% t(embed2)), 0))
    if (max(dm) > 0) dm <- dm / max(dm)

    cost_matrix <- alpha * cost_cos + (1 - alpha) * dm
  } else {
    # Efficient Euclidean distance computation
    embed1_sq_norms <- rowSums(embed1^2)
    embed2_sq_norms <- rowSums(embed2^2)
    cross_product <- embed1 %*% t(embed2)

    cost_sq <- outer(embed1_sq_norms, embed2_sq_norms, "+") - 2 * cross_product
    cost_matrix <- sqrt(pmax(cost_sq, 0))
  }
  
  # Assignment solver with auction option (your suggestion)
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("clue package is required for assignment computation. ",
         "Please install it with: install.packages('clue')", call. = FALSE)
  }
  
  method_arg <- switch(solver_method,
                       auction = "auction",
                       hungarian = "dense",
                       linear = NULL)
  assignment <- solve_lsap_with_padding(cost_matrix, method = method_arg)

  list(assignment = assignment,
       cost_matrix = cost_matrix,
       mapping_matrix = C)
}

solve_lsap_with_padding <- function(cost_matrix, method = NULL) {
  cost_matrix <- as.matrix(cost_matrix)
  nr <- nrow(cost_matrix)
  nc <- ncol(cost_matrix)
  if (nr == 0 || nc == 0) {
    stop("Cost matrix must have positive dimension", call. = FALSE)
  }
  solve_core <- function(mat) {
    if (is.null(method)) {
      return(clue::solve_LSAP(mat))
    }
    tryCatch(
      clue::solve_LSAP(mat, method = method),
      error = function(e) {
        if (grepl("unused argument", e$message, fixed = TRUE)) {
          clue::solve_LSAP(mat)
        } else {
          stop(e)
        }
      }
    )
  }

  if (nr == nc) {
    sol <- solve_core(cost_matrix)
    return(as.integer(sol))
  }

  dim_pad <- max(nr, nc)
  pad_value <- max(cost_matrix)
  if (!is.finite(pad_value)) {
    pad_value <- 0
  }
  padded <- matrix(pad_value, nrow = dim_pad, ncol = dim_pad)
  padded[seq_len(nr), seq_len(nc)] <- cost_matrix

  sol <- solve_core(padded)

  as.integer(sol[seq_len(nr)])
}

#' @export
grasp <- function(data, ...) {
  UseMethod("grasp")
}

#' @rdname grasp
#' @method grasp default
#' @export
grasp.default <- function(data, ...) {
  stop("grasp() requires a hyperdesign object. ",
       "Got: ", paste(class(data), collapse = ", "), 
       ". See ?grasp for usage examples.", call. = FALSE)
}
