lsejd <- function(hd, k_neighbors=5, intrinsic_dims, C_mats, B_mats,
                  t_c=2, t_b=2, alpha_val=1, beta_val=1,
                  delta=0.9, gamma=0.9, mu_c=1, mu_b=1,
                  k_dim=2, k_prime=NULL,
                  kernel_fn=function(x,y) exp(-sum((x-y)^2)),
                  epsilon=1e-5, max_iter=100, lr=1e-3, verbose=FALSE,
                  useStiefel=FALSE, # NEW: whether to use Stiefel manifold optimization
                  stiefelMethod="bb", # "bb" or "curvilinear"
                  stiefelParams=list(rho=0.1, # used if bb
                                     eta=1,
                                     rho1=0.1, # used if curvilinear
                                     rho2=0.9,
                                     tau=1)) 
{
  library(Matrix)
  library(RSpectra)
  
  X_blocks <- xdata(hd)
  m <- length(X_blocks)
  if (is.null(k_prime)) k_prime <- max(k_dim, 10)
  n_i <- sapply(X_blocks, nrow)
  
  # Compute Laplacians and eigen-decompositions
  if (verbose) cat("Computing Laplacians and eigenbases...\n")
  Lap_list <- lapply(X_blocks, function(X) {
    D2 <- as.matrix(dist(X))
    sigmas <- apply(D2,1,median) + .Machine$double.eps
    W <- exp(-D2/(outer(sigmas,sigmas,"*")))
    diag(W) <- 0
    Deg <- rowSums(W)
    L <- diag(Deg) - W
    Matrix(L, sparse=TRUE)
  })
  
  U_list <- vector("list", m)
  Lambda_list <- vector("list", m)
  for (i in seq_len(m)) {
    L <- Lap_list[[i]]
    ev <- suppressWarnings(RSpectra::eigs(L, k=k_prime, which="SM"))
    U_list[[i]] <- ev$vectors
    Lambda_list[[i]] <- diag(ev$values)
  }
  
  # Compute tangent spaces
  if (verbose) cat("Computing tangent spaces...\n")
  T_spaces <- vector("list", m)
  for (i in seq_len(m)) {
    X <- X_blocks[[i]]
    nn_idx <- t(apply(X,1,function(xx) {
      distx <- rowSums((X - matrix(xx,nrow=nrow(X),ncol=ncol(X),byrow=TRUE))^2)
      o <- order(distx)
      o[2:(k_neighbors+1)]
    }))
    xi_i <- intrinsic_dims[i]
    T_i <- vector("list", n_i[i])
    for (r in 1:n_i[i]) {
      nbrs <- nn_idx[r,]
      Y <- scale(X[c(r,nbrs), ], center=TRUE, scale=FALSE)
      CovY <- crossprod(Y)
      eigY <- eigen(CovY, symmetric=TRUE)
      T_i[[r]] <- eigY$vectors[,1:xi_i,drop=FALSE]
    }
    T_spaces[[i]] <- T_i
  }
  
  if (verbose) cat("Extending matching/mismatching info...\n")
  extend_matching <- function(Cij, T_i, X_i, t_c, alpha_val) {
    if (ncol(Cij)==0) return(Cij)
    q <- ncol(Cij)
    n_i <- nrow(Cij)
    Cij_ext <- Matrix(0, nrow=n_i, ncol=q*t_c, sparse=TRUE)
    for (l in 1:q) {
      e <- which(Cij[,l]!=0)
      if (length(e)==0) next
      e <- e[1]
      Te <- T_i[[e]]
      xe <- X_i[e,]
      Z <- (X_i - matrix(xe,n_i,ncol(X_i),byrow=TRUE)) %*% Te
      dist_tangent <- rowSums(Z^2)
      dist_tangent[e] <- Inf
      nn <- order(dist_tangent)[1:t_c]
      for (r_ in 1:t_c) {
        Cij_ext[nn[r_], l+(r_-1)*q] <- alpha_val
      }
    }
    Cij_ext
  }
  
  extend_mismatching <- function(Cij, T_i, X_i, t_b, beta_val) {
    if (ncol(Cij)==0 || t_b==0) return(Matrix(0,nrow=nrow(Cij),ncol=0,sparse=TRUE))
    q <- ncol(Cij)
    n_i <- nrow(Cij)
    Bij_ext <- Matrix(0, nrow=n_i, ncol=q*t_b, sparse=TRUE)
    for (l in 1:q) {
      e <- which(Cij[,l]!=0)
      if (length(e)==0) next
      e <- e[1]
      Te <- T_i[[e]]
      xe <- X_i[e,]
      Z <- (X_i - matrix(xe,n_i,ncol(X_i),byrow=TRUE)) %*% Te
      dist_tangent <- rowSums(Z^2)
      dist_tangent[e] <- -Inf
      ff <- order(dist_tangent,decreasing=TRUE)[1:t_b]
      for (p_ in 1:t_b) {
        Bij_ext[ff[p_], l+(p_-1)*q] <- beta_val
      }
    }
    Bij_ext
  }
  
  Cprime_mats <- vector("list", m)
  Bprime_mats <- vector("list", m)
  for (i in seq_len(m)) {
    Cprime_mats[[i]] <- vector("list", m)
    Bprime_mats[[i]] <- vector("list", m)
  }
  
  for (i in seq_len(m)) {
    X_i <- X_blocks[[i]]
    T_i <- T_spaces[[i]]
    for (j in seq_len(m)) {
      if (i==j) {
        Cprime_mats[[i]][[j]] <- Matrix(0, nrow=n_i[i], ncol=0, sparse=TRUE)
        Bprime_mats[[i]][[j]] <- Matrix(0, nrow=n_i[i], ncol=0, sparse=TRUE)
      } else {
        Cij <- C_mats[[i]][[j]]
        if (is.null(Cij)) Cij <- Matrix(0, n_i[i],0,sparse=TRUE)
        Cij_ext <- extend_matching(Cij, T_i, X_i, t_c, alpha_val)
        Bij <- B_mats[[i]][[j]]
        if (is.null(Bij)) Bij <- Matrix(0,n_i[i],0,sparse=TRUE)
        Bij_ext <- extend_mismatching(Cij, T_i, X_i, t_b, beta_val)
        
        Cprime_mats[[i]][[j]] <- Cij_ext
        Bprime_mats[[i]][[j]] <- Bij_ext
      }
    }
  }
  
  # Precompute M, Mprime, N, Nprime as before
  if (verbose) cat("Precomputing helper matrices...\n")
  M <- vector("list", m)
  Mprime <- vector("list", m)
  N <- vector("list", m)
  Nprime <- vector("list", m)
  for (i in seq_len(m)) {
    M[[i]] <- vector("list", m)
    Mprime[[i]] <- vector("list", m)
    N[[i]] <- vector("list", m)
    Nprime[[i]] <- vector("list", m)
    for (j in seq_len(m)) {
      Cij <- C_mats[[i]][[j]]; if (is.null(Cij)) Cij <- Matrix(0, n_i[i],0,sparse=TRUE)
      M[[i]][[j]] <- t(Cij) %*% U_list[[i]]
      
      Cijp <- Cprime_mats[[i]][[j]]
      Mprime[[i]][[j]] <- t(Cijp) %*% U_list[[i]]
      
      Bij <- B_mats[[i]][[j]]; if (is.null(Bij)) Bij <- Matrix(0,n_i[i],0,sparse=TRUE)
      N[[i]][[j]] <- t(Bij) %*% U_list[[i]]
      
      Bijp <- Bprime_mats[[i]][[j]]
      Nprime[[i]][[j]] <- t(Bijp) %*% U_list[[i]]
    }
  }
  
  # Combine A_i into a single V:
  # V is (m*k_prime) x k_dim, block by block vertically stacked
  stackA <- function(A_list) {
    do.call(rbind, A_list)
  }
  
  # Extract A_i from V
  getA <- function(V, i) {
    start <- (i-1)*k_prime+1
    end <- i*k_prime
    V[start:end, , drop=FALSE]
  }
  
  # Define objective and gradient in terms of V
  F_stiefel <- function(V) {
    # Extract A_i from V
    A_list_cur <- lapply(seq_len(m), function(i) getA(V,i))
    # Compute objective with current A_i's
    # This replicates the logic from the previous objective_and_gradient function, but only objective
    # We'll reuse code from that logic:
    # 1) Eigen term
    obj_eig <- 0
    for (i in seq_len(m)) {
      Ai <- A_list_cur[[i]]
      Ri <- t(Ai) %*% Lambda_list[[i]] %*% Ai - Lambda_list[[i]][1:k_prime,1:k_prime]
      obj_eig <- obj_eig + sum(Ri^2)
    }
    
    # 2) Matching and mismatching terms
    obj_c <- 0
    obj_cp <- 0
    obj_b <- 0
    obj_bp <- 0
    
    for (i in seq_len(m)) {
      Ai <- A_list_cur[[i]]
      for (j in seq_len(m)) {
        if (i==j) next
        Aj <- A_list_cur[[j]]
        Xc <- M[[i]][[j]] %*% Ai - M[[j]][[i]] %*% Aj
        obj_c <- obj_c + sum(Xc^2)
        Xcp <- Mprime[[i]][[j]] %*% Ai - Mprime[[j]][[i]] %*% Aj
        obj_cp <- obj_cp + sum(Xcp^2)
        
        Xb <- N[[i]][[j]] %*% Ai - N[[j]][[i]] %*% Aj
        obj_b <- obj_b + sum(Xb^2)
        Xbp <- Nprime[[i]][[j]] %*% Ai - Nprime[[j]][[i]] %*% Aj
        obj_bp <- obj_bp + sum(Xbp^2)
      }
    }
    
    obj <- obj_eig + mu_c*(obj_c + delta*obj_cp) - mu_b*(obj_b + gamma*obj_bp)
    return(obj)
  }
  
  dF_stiefel <- function(V) {
    # Compute gradient w.r.t V
    # Similar approach:
    A_list_cur <- lapply(seq_len(m), function(i) getA(V,i))
    # Grad arrays w.r.t each A_i
    grad_eig <- lapply(seq_len(m), function(i) matrix(0,k_prime,k_dim))
    grad_c <- lapply(seq_len(m), function(i) matrix(0,k_prime,k_dim))
    grad_cp <- lapply(seq_len(m), function(i) matrix(0,k_prime,k_dim))
    grad_b <- lapply(seq_len(m), function(i) matrix(0,k_prime,k_dim))
    grad_bp <- lapply(seq_len(m), function(i) matrix(0,k_prime,k_dim))
    
    # Eigen grad
    for (i in seq_len(m)) {
      Ai <- A_list_cur[[i]]
      Ri <- Ai^T %*% Lambda_list[[i]] %*% Ai - Lambda_list[[i]][1:k_prime,1:k_prime]
      grad_eig[[i]] <- 4 * Lambda_list[[i]] %*% Ai %*% Ri
    }
    
    # C, Cprime, B, Bprime grads
    for (i in seq_len(m)) {
      Ai <- A_list_cur[[i]]
      for (j in seq_len(m)) {
        if (i==j) next
        Aj <- A_list_cur[[j]]
        # C
        Xc <- M[[i]][[j]] %*% Ai - M[[j]][[i]] %*% Aj
        grad_c[[i]] <- grad_c[[i]] + 2 * t(M[[i]][[j]]) %*% Xc
        grad_c[[j]] <- grad_c[[j]] + 2 * t(M[[j]][[i]]) %*% (M[[j]][[i]] %*% Aj - M[[i]][[j]] %*% Ai)
        
        # Cprime
        Xcp <- Mprime[[i]][[j]] %*% Ai - Mprime[[j]][[i]] %*% Aj
        grad_cp[[i]] <- grad_cp[[i]] + 2 * t(Mprime[[i]][[j]]) %*% Xcp
        grad_cp[[j]] <- grad_cp[[j]] + 2 * t(Mprime[[j]][[i]]) %*% (Mprime[[j]][[i]] %*% Aj - Mprime[[i]][[j]] %*% Ai)
        
        # B
        Xb <- N[[i]][[j]] %*% Ai - N[[j]][[i]] %*% Aj
        grad_b[[i]] <- grad_b[[i]] + 2 * t(N[[i]][[j]]) %*% Xb
        grad_b[[j]] <- grad_b[[j]] + 2 * t(N[[j]][[i]]) %*% (N[[j]][[i]] %*% Aj - N[[i]][[j]] %*% Ai)
        
        # Bprime
        Xbp <- Nprime[[i]][[j]] %*% Ai - Nprime[[j]][[i]] %*% Aj
        grad_bp[[i]] <- grad_bp[[i]] + 2 * t(Nprime[[i]][[j]]) %*% Xbp
        grad_bp[[j]] <- grad_bp[[j]] + 2 * t(Nprime[[j]][[i]]) %*% (Nprime[[j]][[i]] %*% Aj - Nprime[[i]][[j]] %*% Ai)
      }
    }
    
    # Combine
    grad_final_blocks <- lapply(seq_len(m), function(i) {
      grad_eig[[i]] +
        mu_c*(grad_c[[i]] + delta*grad_cp[[i]]) -
        mu_b*(grad_b[[i]] + gamma*grad_bp[[i]])
    })
    
    # Stack vertically to get gradient w.r.t V
    grad_final <- do.call(rbind, grad_final_blocks)
    grad_final
  }
  
  # Initialize A_i (and thus V)
  A_init <- lapply(seq_len(m), function(i) {
    mat <- diag(k_prime)
    mat[,1:k_dim]
  })
  Vinit <- do.call(rbind, A_init)
  
  if (!useStiefel) {
    # If not using Stiefel manifold optimization:
    # Perform simple gradient descent in Euclidean space
    if (verbose) cat("Running simple gradient descent in Euclidean space...\n")
    prev_obj <- Inf
    for (iter_ in 1:max_iter) {
      g <- dF_stiefel(Vinit)
      current_obj <- F_stiefel(Vinit)
      if (verbose && iter_ %% 10 == 0) cat("Iter:", iter_, "Obj:", current_obj, "\n")
      if (abs(prev_obj - current_obj) < epsilon) {
        if (verbose) cat("Converged.\n")
        break
      }
      prev_obj <- current_obj
      Vinit <- Vinit - lr*g
      # Re-orthonormalize columns of each block if we want to maintain orthonormality
      # This step is optional if we didn't strictly require each A_i be on stiefel manifold
      # If needed: re-orthonormalize each A_i block
      for (i in seq_len(m)) {
        Ai <- getA(Vinit,i)
        # QR decomposition to ensure orthonormality
        q <- qr.Q(qr(Ai))
        Vinit[( (i-1)*k_prime+1 ):(i*k_prime), ] <- q
      }
    }
    
    V_final <- Vinit
  } else {
    # Use Stiefel manifold optimization
    if (verbose) cat("Running Stiefel manifold optimization...\n")
    
    # optStiefel expects F and dF with signature F(V) and dF(V).
    # We already have F_stiefel and dF_stiefel.
    # set up search parameters
    if (stiefelMethod == "bb") {
      # For lineSearchBB we need rho and maybe C, we can guess C = F(V)
      # We'll just pass these as searchParams
      if (!("rho" %in% names(stiefelParams))) stiefelParams$rho <- 0.1
      if (!("C" %in% names(stiefelParams))) stiefelParams$C <- F_stiefel(Vinit)
      V_final <- optStiefel(F_stiefel, dF_stiefel, Vinit, method="bb",
                            searchParams=stiefelParams,
                            tol=epsilon, maxIters=max_iter, verbose=verbose)
    } else {
      # Curvilinear search
      # expects rho1, rho2, tau in searchParams
      if (!("rho1" %in% names(stiefelParams))) stiefelParams$rho1 <- 0.1
      if (!("rho2" %in% names(stiefelParams))) stiefelParams$rho2 <- 0.9
      if (!("tau" %in% names(stiefelParams))) stiefelParams$tau <- 1
      V_final <- optStiefel(F_stiefel, dF_stiefel, Vinit, method="curvilinear",
                            searchParams=stiefelParams,
                            tol=epsilon, maxIters=max_iter, verbose=verbose)
    }
  }
  
  # Extract final A_i and compute V_list
  # The final eigenbasis per modality:
  A_final_list <- lapply(seq_len(m), function(i) getA(V_final,i))
  V_list <- lapply(seq_len(m), function(i) U_list[[i]] %*% A_final_list[[i]])
  
  if (verbose) cat("Done.\n")
  list(V_list=V_list)
}