#' Kernel Manifold Alignment (KEMA)
#'
#' Performs Kernel Manifold Alignment for supervised/semi-supervised domain 
#' adaptation. Projects data from multiple domains into a shared latent space.
#'
#' KEMA is designed for multi-domain data where you want to find a common 
#' representation that preserves both the intrinsic geometry of each domain and 
#' the class structure across domains. It supports semi-supervised learning 
#' with missing labels (NA values).
#'
#' The algorithm offers two solver methods:
#' - "regression": Fast approximation using spectral regression (default). This 
#'   method first solves the eigenvalue problem on graph Laplacians, then uses 
#'   ridge regression to find kernel coefficients. It's much faster but may be 
#'   less accurate for non-linear kernels.
#' - "exact": Precise solution using the correct generalized eigenvalue 
#'   formulation. This method solves the mathematically correct KEMA 
#'   optimization problem but is more computationally intensive, especially for 
#'   large datasets.
#'
#' @param data Input data object. Can be a hyperdesign object (for 
#'   \code{kema.hyperdesign}) or a multidesign object (for 
#'   \code{kema.multidesign})
#' @param y Name of the label variable to use for alignment. Can contain NA 
#'   values for unlabeled samples in semi-supervised learning scenarios
#' @param ... Additional arguments passed to specific methods. See 
#'   \code{\link{kema.hyperdesign}} and \code{\link{kema.multidesign}} for 
#'   details on method-specific parameters such as \code{preproc}, 
#'   \code{ncomp}, \code{knn}, \code{sigma}, \code{u}, \code{kernel}, 
#'   \code{sample_frac}, \code{use_laplacian}, \code{solver}, \code{dweight}, 
#'   \code{rweight}, \code{simfun}, \code{disfun}, \code{lambda}, and for 
#'   \code{kema.multidesign} the additional \code{subject} parameter
#'
#' @details
#' KEMA solves the following optimization problem:
#' \deqn{\max_{\alpha} \frac{\alpha^T K A K^T \alpha}{\alpha^T K B K^T \alpha}}
#' 
#' where:
#' \itemize{
#'   \item \code{A = u*L + (1-u)*Ls} captures manifold structure (L) and 
#'     same-class alignment (Ls)
#'   \item \code{B = rweight*Lr + dweight*Ld + lambda*I} captures class 
#'     separation and regularization
#'   \item \code{K} is the block-diagonal kernel matrix across domains
#' }
#'
#' The trade-off parameter \code{u} controls the balance between preserving 
#' manifold geometry and enforcing class alignment. The solver parameter 
#' determines the computational approach:
#' \itemize{
#'   \item "regression": Fast two-step approximation (eigendecomposition + 
#'     ridge regression)
#'   \item "exact": Direct solution of the generalized eigenvalue problem
#' }
#'
#' For large datasets, use \code{sample_frac < 1} to enable REKEMA, which uses 
#' landmark points to reduce computational complexity from O(n^2) to O(r^2) where 
#' r is the number of landmarks.
#'
#' @return A \code{multiblock_biprojector} object containing:
#' \itemize{
#'   \item \code{s}: Scores (embedded coordinates) for all samples
#'   \item \code{v}: Primal vectors (feature weights) for out-of-sample 
#'     projection
#'   \item \code{sdev}: Standard deviations of the components
#'   \item \code{alpha}: Dual coefficients in kernel space
#'   \item Additional metadata for reconstruction and validation
#' }
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign data
#' library(multivarious)
#' 
#' # Create synthetic multi-domain data
#' domain1 <- list(x = matrix(rnorm(100), 50, 2), 
#'                labels = sample(c("A", "B"), 50, TRUE))
#' domain2 <- list(x = matrix(rnorm(100), 50, 2), 
#'                labels = sample(c("A", "B"), 50, TRUE))
#' hd <- list(domain1 = domain1, domain2 = domain2)
#' 
#' # Run KEMA with default settings
#' result <- kema(hd, y = labels, ncomp = 2)
#' 
#' # Semi-supervised learning with missing labels
#' domain1$labels[1:10] <- NA  # Mark some samples as unlabeled
#' result_semi <- kema(hd, y = labels, ncomp = 2)
#' 
#' # Use exact solver for highest accuracy
#' result_exact <- kema(hd, y = labels, solver = "exact", ncomp = 2)
#' 
#' # Use REKEMA for large datasets
#' result_rekema <- kema(hd, y = labels, sample_frac = 0.5, ncomp = 2)
#' }
#'
#' @references
#' Tuia, D., & Camps-Valls, G. (2016). Kernel manifold alignment for domain 
#' adaptation. PLoS ONE, 11(2), e0148655.
#'
#' @seealso \code{\link{kema.hyperdesign}}, \code{\link{kema.multidesign}}
#' @export
kema <- function(data, y, ...) {
  UseMethod("kema")
}

#' Generalized Orthogonal Procrustes Alignment
#'
#' Performs Generalized Orthogonal Procrustes alignment to find orthogonal 
#' transformations. Aligns data from multiple domains by minimizing squared 
#' differences between corresponding task observations.
#'
#' This method extends classical Procrustes analysis to handle partial task 
#' observations, where each domain may observe only a subset of a global set of 
#' tasks. The algorithm uses an efficient Generalized Power Method (GPM) with 
#' sparse matrix operations and robust initialization to find optimal orthogonal 
#' transformations.
#'
#' @param data Input data object. For the hyperdesign method, this should be a 
#'   hyperdesign object containing multiple data domains
#' @param y Name of the task/label variable to use for alignment across domains
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' The Generalized Procrustes problem seeks to find orthogonal matrices 
#' \eqn{O_i} for each domain \eqn{i} that minimize:
#' \deqn{\sum_{i,j} \sum_{k \in T_{ij}} ||O_i^T A_i(:,k) - O_j^T A_j(:,k)||^2}
#' 
#' where \eqn{T_{ij}} represents the set of tasks observed by both domains 
#' \eqn{i} and \eqn{j}.
#'
#' The algorithm handles several key challenges:
#' \itemize{
#'   \item \strong{Partial observations}: Each domain may observe different 
#'     subsets of tasks
#'   \item \strong{Sparse structure}: Uses efficient sparse matrix operations 
#'     for scalability
#'   \item \strong{Robust initialization}: SVD-based initialization with 
#'     fallback to random orthogonal matrices
#'   \item \strong{Convergence guarantees}: Monotonic improvement with 
#'     configurable tolerance
#' }
#'
#' Key features:
#' \itemize{
#'   \item Vectorized operations using sparse matrix algebra
#'   \item Efficient projection onto the orthogonal group O(d)
#'   \item Optional tightness certificate for global optimality validation
#'   \item Flexible preprocessing and multiple input formats
#' }
#'
#' @return The return value depends on the specific method:
#' \itemize{
#'   \item For hyperdesign objects: A list containing orthogonal transformation 
#'     matrices, consensus matrix, convergence information, and domain metadata
#'   \item For direct matrix input: A list with transformation matrices and 
#'     alignment results
#' }
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign data
#' library(multidesign)
#' 
#' # Create example domains with partial task overlap
#' d1_data <- matrix(rnorm(50), 5, 10)  # 5 tasks x 10 features
#' d1_design <- data.frame(task = factor(c("A", "B", "C", "D", "E")))
#' d1 <- multidesign(d1_data, d1_design)
#' 
#' d2_data <- matrix(rnorm(40), 4, 10)  # 4 tasks x 10 features
#' d2_design <- data.frame(task = factor(c("A", "C", "D", "F")))
#' d2 <- multidesign(d2_data, d2_design)
#' 
#' # Create hyperdesign
#' hd <- hyperdesign(list(domain1 = d1, domain2 = d2))
#' 
#' # Perform alignment
#' result <- generalized_procrustes(hd, task)
#' 
#' # Check convergence and results
#' print(result$converged)
#' print(dim(result$A_est))  # Features x total tasks
#' }
#'
#' @references
#' Gower, J. C. (1975). Generalized procrustes analysis. Psychometrika, 40(1), 
#' 33-51.
#' 
#' Ten Berge, J. M. F. (1977). Orthogonal procrustes rotation for two or more 
#' matrices. Psychometrika, 42(2), 267-276.
#'
#' @seealso \code{\link{generalized_procrustes.hyperdesign}}
#' @export
generalized_procrustes <- function(data, ...) {
  UseMethod("generalized_procrustes")
}
