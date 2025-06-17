#' Node Identity Extraction
#'
#' This function extracts node identity features based on the degree distributions
#' of k-hop neighborhoods and optional node attributes.
#'
#' @param A An adjacency matrix representing the graph.
#' @param k The maximum number of hops to consider.
#' @param attributes A matrix of node attributes (optional).
#'
#' @return A matrix of node identity features, where each row corresponds to a node in the graph.
extract_node_identity <- function(A, k, attributes = NULL) {
  G <- graph_from_adjacency_matrix(A, mode = "undirected")
  
  # Compute maximum degree in the graph
  max_degree <- max(degree(G))
  
  # Compute number of logarithmic bins
  num_bins <- ceiling(log2(max_degree))
  
  # Compute degree distributions of k-hop neighborhoods with logarithmic binning
  degree_distributions <- lapply(1:k, function(k) {
    neighborhood_degrees <- neighborhood(G, order = k, mode = "all", mindist = k)
    sapply(neighborhood_degrees, function(x) {
      binned_degrees <- floor(log2(degree(G, v = x) + 1))
      degree_counts <- tabulate(binned_degrees, nbins = num_bins)
      as.vector(degree_counts)
    })
  })
  
  # Concatenate degree distributions and optional attributes
  identity_features <- do.call(cbind, degree_distributions)
  if (!is.null(attributes)) {
    identity_features <- cbind(identity_features, attributes)
  }
  
  return(identity_features)
}

#' xNetMF: Cross-Network Matrix Factorization
#'
#' This function learns node representations using the xNetMF method, which efficiently
#' factorizes an approximated similarity matrix based on the Nyström method.
#'
#' @param X A matrix of node identity features, where each row corresponds to a node.
#' @param p The number of landmark nodes to use.
#' @param d The dimensionality of the learned node representations.
#'
#' @return A matrix of learned node representations, where each row corresponds to a node.
xNetMF <- function(X, p, d) {
  # Randomly select landmark nodes
  landmark_indices <- sample(1:nrow(X), p)
  
  # Compute node-to-landmark similarity matrix C
  C <- exp(-0.5 * as.matrix(dist2(X, X[landmark_indices,])))
  
  # Compute landmark-to-landmark similarity matrix W
  W <- C[landmark_indices,]
  
  # Compute SVD of the pseudoinverse of W
  svd_W <- svd(ginv(W))
  
  # Obtain node embeddings
  embeddings <- C %*% svd_W$u %*% diag(sqrt(svd_W$d[1:d]))
  
  return(embeddings)
}

#' Node Representation Alignment
#'
#' This function aligns node representations across two graphs using a k-d tree.
#'
#' @param emb1 A matrix of node representations for the first graph, where each row corresponds to a node.
#' @param emb2 A matrix of node representations for the second graph, where each row corresponds to a node.
#'
#' @return A vector indicating the aligned node indices from the second graph for each node in the first graph.
align_embeddings <- function(emb1, emb2) {
  # Build a k-d tree for efficient nearest neighbor search
  kd_tree <- nn2(data = emb2, query = emb1, k = 1)
  
  # Extract aligned node indices
  alignment <- kd_tree$nn.idx[,1]
  
  return(alignment)
}

#' REGAL: REpresentation-based Graph ALignment
#'
#' This function implements the REGAL algorithm for aligning nodes across two graphs
#' based on their learned representations. It consists of three main steps:
#' 1. Node identity extraction
#' 2. Efficient similarity-based representation learning using xNetMF
#' 3. Fast node representation alignment using a k-d tree
#'
#' @param A1 An adjacency matrix representing the first graph.
#' @param A2 An adjacency matrix representing the second graph.
#' @param k The maximum number of hops to consider for node identity extraction.
#' @param p The number of landmark nodes to use in xNetMF.
#' @param d The dimensionality of the learned node representations.
#' @param attributes1 A matrix of node attributes for the first graph (optional).
#' @param attributes2 A matrix of node attributes for the second graph (optional).
#'
#' @return A list with two elements:
#' - embeddings: A matrix of aligned node representations, where each row corresponds to a node in G1 or G2.
#' - alignment: A vector indicating the aligned node indices from G2 for each node in G1.
regal_align <- function(A1, A2, k = 2, p = 100, d = 128, attributes1 = NULL, attributes2 = NULL) {
  # Step 1: Node Identity Extraction
  identity1 <- extract_node_identity(A1, k, attributes1)
  identity2 <- extract_node_identity(A2, k, attributes2)
  
  # Step 2: Efficient Similarity-based Representation Learning
  embeddings <- xNetMF(rbind(identity1, identity2), p, d)
  
  # Step 3: Fast Node Representation Alignment
  alignment <- align_embeddings(embeddings[1:nrow(A1),], embeddings[-(1:nrow(A1)),])
  
  return(list(embeddings = embeddings, alignment = alignment))
}