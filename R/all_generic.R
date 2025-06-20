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

#' Graph Alignment by Spectral Corresponding Functions (GRASP)
#'
#' Performs GRASP alignment on graph data structures. GRASP tackles the "unrestricted" 
#' graph alignment problem where only graph structures are given, finding optimal 
#' node-to-node correspondences between graphs by matching functions defined on nodes.
#'
#' The core idea of GRASP is to transfer functional maps from 3D shape analysis to 
#' graph alignment. Instead of directly matching nodes (which is hard), it matches 
#' functions defined on the nodes. By finding a robust mapping between basis functions 
#' on each graph, GRASP deduces the underlying node-to-node alignment as a by-product.
#'
#' @param data Input data object containing graph domains
#' @param ... Additional arguments passed to specific methods. See 
#'   \code{\link{grasp.hyperdesign}} for details on method-specific parameters 
#'   such as \code{preproc}, \code{ncomp}, \code{q_descriptors}, \code{sigma}, 
#'   \code{lambda}, \code{use_laplacian}, and \code{solver}
#'
#' @details
#' GRASP operates through five conceptual blocks:
#' \itemize{
#'   \item \strong{Block A}: Spectral basis construction using eigendecomposition 
#'     of graph Laplacians
#'   \item \strong{Block B}: Multi-scale node descriptors using heat kernel or 
#'     Personalized PageRank
#'   \item \strong{Block C}: Base alignment using Stiefel manifold optimization
#'   \item \strong{Block D}: Functional correspondence mapping
#'   \item \strong{Block E}: Final node-to-node assignment via linear assignment
#' }
#'
#' The algorithm provides a global, multi-scale view of graph structure that is 
#' robust to noise and structural perturbations. For the pairwise case, complexity 
#' is dominated by the final assignment step (typically O(n³)).
#'
#' Key parameters:
#' \itemize{
#'   \item \code{ncomp}: Number of eigenvectors used (dimension of spectral basis)
#'   \item \code{q_descriptors}: Number of descriptor functions for multi-scale analysis
#'   \item \code{sigma}: Diffusion parameter controlling descriptor time scales
#'   \item \code{lambda}: Regularization balancing eigen-structure vs descriptor alignment
#' }
#'
#' @return A \code{multiblock_biprojector} object containing:
#' \itemize{
#'   \item \code{s}: Aligned spectral embeddings for all nodes
#'   \item \code{v}: Primal vectors for out-of-sample projection
#'   \item \code{assignment}: Node-to-node correspondence vector
#'   \item \code{rotation}: Orthogonal rotation matrix between spectral bases
#'   \item \code{mapping_matrix}: Diagonal functional correspondence matrix
#'   \item Additional metadata for reconstruction and validation
#' }
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign graph data
#' library(multidesign)
#' 
#' # Create synthetic graph domains
#' set.seed(123)
#' domain1 <- list(
#'   x = matrix(rnorm(100), 50, 2),  # Node features
#'   design = data.frame(node_id = 1:50)
#' )
#' domain2 <- list(
#'   x = matrix(rnorm(100), 50, 2),
#'   design = data.frame(node_id = 1:50)  
#' )
#' 
#' # Create hyperdesign
#' hd <- hyperdesign(list(domain1 = domain1, domain2 = domain2))
#' 
#' # Run GRASP alignment with default parameters
#' result <- grasp(hd, ncomp = 20, q_descriptors = 50)
#' 
#' # Access alignment results
#' node_assignment <- result$assignment
#' aligned_embeddings <- result$s
#' 
#' # Use different parameters for robustness
#' result_robust <- grasp(hd, ncomp = 30, sigma = 1.0, lambda = 0.2)
#' }
#'
#' @references
#' Bogo, F., Romero, J., Loper, M., & Black, M. J. (2016). FAUST: Dataset and 
#' evaluation for 3D mesh registration. In Proceedings of the IEEE Conference 
#' on Computer Vision and Pattern Recognition (pp. 3794-3801).
#'
#' Ovsjanikov, M., Ben-Chen, M., Solomon, J., Butscher, A., & Guibas, L. (2012). 
#' Functional maps: a flexible representation of maps between shapes. ACM 
#' Transactions on Graphics, 31(4), 1-11.
#'
#' @seealso \code{\link{grasp.hyperdesign}}
#' @export
grasp <- function(data, ...) {
  UseMethod("grasp")
}

#' CONE-Align: Consensus Optimization for Node Embedding Alignment
#'
#' Performs CONE-Align on hyperdesign data structures. Aligns graph embeddings
#' by alternating between orthogonal transformations and node assignments.
#'
#' CONE-Align tackles graph alignment by iteratively refining both orthogonal
#' transformations and node-to-node correspondences. The algorithm alternates
#' between two steps: (1) finding optimal orthogonal rotations given current
#' assignments via Procrustes analysis, and (2) updating assignments given
#' current transformations via linear assignment.
#'
#' @param data Input data object containing graph domains
#' @param ... Additional arguments passed to specific methods. See 
#'   \code{\link{cone_align.hyperdesign}} for details on method-specific 
#'   parameters such as \code{preproc}, \code{ncomp}, \code{sigma}, 
#'   \code{lambda}, \code{use_laplacian}, \code{solver}, \code{max_iter}, 
#'   and \code{tol}
#'
#' @details
#' CONE-Align operates through the following algorithmic blocks:
#' \itemize{
#'   \item \strong{Spectral Embedding}: Compute Laplacian eigenmaps for each graph
#'   \item \strong{Orthogonal Alignment}: Find rotation matrices via Procrustes analysis
#'   \item \strong{Assignment Update}: Solve linear assignment problem for correspondences
#'   \item \strong{Convergence Check}: Iterate until assignments stabilize
#' }
#'
#' The algorithm minimizes the objective:
#' \deqn{\sum_{i,j} ||Q_i^T Z_i - P_{ij} Q_j^T Z_j||_F^2}
#'
#' where \eqn{Z_i} are the embeddings, \eqn{Q_i} are orthogonal transforms,
#' and \eqn{P_{ij}} are permutation matrices.
#'
#' Key parameters:
#' \itemize{
#'   \item \code{ncomp}: Dimension of spectral embeddings
#'   \item \code{sigma}: Bandwidth for embedding computation
#'   \item \code{lambda}: Regularization for numerical stability
#'   \item \code{solver}: Assignment algorithm ("linear" or "auction")
#' }
#'
#' @return A \code{multiblock_biprojector} object containing:
#' \itemize{
#'   \item \code{s}: Aligned embeddings for all nodes
#'   \item \code{v}: Primal vectors for out-of-sample projection
#'   \item \code{assignment}: Node-to-node correspondence assignments
#'   \item \code{rotation}: Orthogonal transformation matrices
#'   \item \code{sdev}: Standard deviations of aligned components
#'   \item Additional metadata for reconstruction and validation
#' }
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign graph data
#' library(multidesign)
#' 
#' # Create synthetic graph domains
#' set.seed(123)
#' domain1 <- list(
#'   x = matrix(rnorm(100), 50, 2),  # Node features
#'   design = data.frame(node_id = 1:50)
#' )
#' domain2 <- list(
#'   x = matrix(rnorm(100), 50, 2),
#'   design = data.frame(node_id = 1:50)  
#' )
#' 
#' # Create hyperdesign
#' hd <- hyperdesign(list(domain1 = domain1, domain2 = domain2))
#' 
#' # Run CONE-Align with default parameters
#' result <- cone_align(hd, ncomp = 10)
#' 
#' # Access alignment results
#' node_assignment <- result$assignment
#' aligned_embeddings <- result$s
#' 
#' # Use exact solver for optimal results
#' result_exact <- cone_align(hd, ncomp = 10, solver = "linear")
#' }
#'
#' @references
#' Heimann, M., Shen, H., Safavi, T., & Koutra, D. (2018). REGAL: Representation 
#' learning-based graph alignment. In Proceedings of the 27th ACM International 
#' Conference on Information and Knowledge Management (pp. 117-126).
#'
#' @seealso \code{\link{cone_align.hyperdesign}}
#' @export
cone_align <- function(data, ...) {
  UseMethod("cone_align")
}

#' Position-Aware Random Transport (PARROT) Network Alignment
#'
#' Performs PARROT alignment on hyperdesign data structures. Aligns networks
#' using regularized optimal transport with position-aware features and consistency constraints.
#'
#' PARROT tackles network alignment by formulating it as a regularized optimal transport
#' problem. The method incorporates position-aware features through Random Walk with Restart
#' (RWR) descriptors and enforces structural consistency through neighborhood-preserving
#' regularization terms.
#'
#' @param data Input data object containing network domains
#' @param anchors Name of anchor/correspondence variable for semi-supervised alignment
#' @param ... Additional arguments passed to specific methods. See 
#'   \code{\link{parrot.hyperdesign}} for details on method-specific parameters 
#'   such as \code{preproc}, \code{ncomp}, \code{sigma}, \code{lambda}, 
#'   \code{tau}, \code{solver}, \code{max_iter}, and \code{tol}
#'
#' @details
#' PARROT operates through the following algorithmic components:
#' \itemize{
#'   \item \strong{Position-Aware Features}: Compute RWR descriptors capturing network position
#'   \item \strong{Cross-Network Cost}: Build transport cost matrix between networks
#'   \item \strong{Consistency Regularization}: Add structural similarity constraints
#'   \item \strong{Optimal Transport}: Solve regularized transport problem via Sinkhorn
#' }
#'
#' The algorithm minimizes the objective:
#' \deqn{L(S) = \langle C, S \rangle + \lambda_1 \Omega_1(S) + \lambda_2 \Omega_2(S) + \tau H(S)}
#'
#' where \eqn{C} is the position-aware cost matrix, \eqn{\Omega_1, \Omega_2} are consistency
#' regularizers, and \eqn{H(S)} is the entropy regularization term with parameter \eqn{\tau}.
#'
#' Key parameters:
#' \itemize{
#'   \item \code{sigma}: RWR restart probability and diffusion parameter
#'   \item \code{lambda}: Consistency regularization weights
#'   \item \code{tau}: Entropy regularization parameter for Sinkhorn
#'   \item \code{solver}: Transport solver ("sinkhorn" or "exact")
#' }
#'
#' @return A \code{multiblock_biprojector} object containing:
#' \itemize{
#'   \item \code{s}: Aligned network embeddings 
#'   \item \code{v}: Primal vectors for out-of-sample projection
#'   \item \code{alignment_matrix}: Soft alignment/transport plan between networks
#'   \item \code{transport_plan}: Dense transport matrix S
#'   \item \code{sdev}: Standard deviations of aligned components
#'   \item Additional metadata for reconstruction and validation
#' }
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign network data
#' library(multidesign)
#' 
#' # Create synthetic network domains
#' set.seed(123)
#' domain1 <- list(
#'   x = matrix(rnorm(200), 100, 2),  # Node features
#'   design = data.frame(
#'     node_id = 1:100,
#'     anchors = c(1:10, rep(NA, 90))  # First 10 nodes are anchors
#'   )
#' )
#' domain2 <- list(
#'   x = matrix(rnorm(200), 100, 2),
#'   design = data.frame(
#'     node_id = 1:100,
#'     anchors = c(1:10, rep(NA, 90))  # Corresponding anchors
#'   )
#' )
#' 
#' # Create hyperdesign
#' hd <- hyperdesign(list(domain1 = domain1, domain2 = domain2))
#' 
#' # Run PARROT alignment with default parameters
#' result <- parrot(hd, anchors = anchors)
#' 
#' # Access alignment results
#' transport_plan <- result$transport_plan
#' aligned_embeddings <- result$s
#' 
#' # Use different regularization settings
#' result_strong <- parrot(hd, anchors = anchors, lambda = 0.5, tau = 0.1)
#' }
#'
#' @references
#' Wang, S., Chen, Z., Yu, X., Li, T., Yang, J., & Liu, X. (2022). PARROT: 
#' Position-aware regularized optimal transport for network alignment. 
#' In Proceedings of the 28th ACM SIGKDD Conference on Knowledge Discovery 
#' and Data Mining (pp. 1896-1905).
#'
#' @seealso \code{\link{parrot.hyperdesign}}
#' @export
parrot <- function(data, anchors, ...) {
  UseMethod("parrot")
}

#' Multiple Graph CONE-Align: Consensus Optimization for Multi-Domain Node Embedding Alignment
#'
#' Performs CONE-Align on three or more graph domains simultaneously. Extends the
#' pairwise CONE-Align algorithm to handle multiple graphs through iterative
#' alignment to a common reference frame.
#'
#' @param data Input data containing three or more graph domains. Can be:
#'   \itemize{
#'     \item A \code{hyperdesign} object with 3+ domains
#'     \item A list of 3+ matrices (nodes x features)
#'   }
#' @param ... Additional arguments passed to fitting methods
#'
#' @return A \code{multiblock_biprojector} object containing alignment results
#'   for all graphs
#'
#' @seealso \code{\link{cone_align}} for pairwise alignment, 
#'   \code{\link{cone_align_multiple.hyperdesign}} for detailed documentation
#' @export
cone_align_multiple <- function(data, ...) {
  UseMethod("cone_align_multiple")
}

#' Multi-Graph GRASP: Graph Alignment by Spectral Corresponding Functions
#'
#' Performs GRASP alignment on three or more graph domains simultaneously. Extends
#' the pairwise GRASP algorithm to handle multiple graphs through joint diagonalization
#' and shared latent basis alignment.
#'
#' @param data Input data containing three or more graph domains. Can be:
#'   \itemize{
#'     \item A \code{hyperdesign} object with 3+ domains
#'     \item A list of 3+ matrices (nodes x features)
#'   }
#' @param ... Additional arguments passed to fitting methods
#'
#' @return A \code{grasp_multiset} object containing alignment results
#'   for all graphs
#'
#' @seealso \code{\link{grasp}} for pairwise alignment, 
#'   \code{\link{grasp_multiset.hyperdesign}} for detailed documentation
#' @export
grasp_multiset <- function(data, ...) {
  UseMethod("grasp_multiset")
}
