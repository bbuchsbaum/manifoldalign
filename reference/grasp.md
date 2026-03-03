# Graph Alignment by Spectral Corresponding Functions (GRASP)

Performs GRASP alignment on graph data structures. GRASP tackles the
"unrestricted" graph alignment problem where only graph structures are
given, finding optimal node-to-node correspondences between graphs by
matching functions defined on nodes.

Performs GRASP alignment on hyperdesign data structures. Projects graph
data from multiple domains into aligned spectral spaces.

## Usage

``` r
grasp(data, ...)

# S3 method for class 'hyperdesign'
grasp(
  data,
  preproc = multivarious::center(),
  ncomp = 30,
  q_descriptors = 100,
  sigma = 0.73,
  lambda = 0.1,
  use_laplacian = TRUE,
  solver = "linear",
  ...
)

# Default S3 method
grasp(data, ...)
```

## Arguments

- data:

  A hyperdesign object containing multiple graph domains

- ...:

  Additional arguments (currently unused)

- preproc:

  Preprocessing function to apply to the data (default: center)

- ncomp:

  Number of spectral components to extract (default: 30)

- q_descriptors:

  Number of diffusion-time descriptors (default: 100)

- sigma:

  Diffusion parameter for descriptor computation (default: 0.73)

- lambda:

  Regularization parameter for base alignment (default: 0.1)

- use_laplacian:

  Whether to use Laplacian normalization (default: TRUE)

- solver:

  Assignment method: "linear" for exact assignment (default)

## Value

A `multiblock_biprojector` object containing:

- `s`: Aligned spectral embeddings for all nodes

- `v`: Primal vectors for out-of-sample projection

- `assignment`: Node-to-node correspondence vector

- `rotation`: Orthogonal rotation matrix between spectral bases

- `mapping_matrix`: Diagonal functional correspondence matrix

- Additional metadata for reconstruction and validation

A multiblock_biprojector object containing the GRASP alignment

## Details

The core idea of GRASP is to transfer functional maps from 3D shape
analysis to graph alignment. Instead of directly matching nodes (which
is hard), it matches functions defined on the nodes. By finding a robust
mapping between basis functions on each graph, GRASP deduces the
underlying node-to-node alignment as a by-product.

GRASP operates through five conceptual blocks:

- **Block A**: Spectral basis construction using eigendecomposition of
  graph Laplacians

- **Block B**: Multi-scale node descriptors using heat kernel or
  Personalized PageRank

- **Block C**: Base alignment using Stiefel manifold optimization

- **Block D**: Functional correspondence mapping

- **Block E**: Final node-to-node assignment via linear assignment

The algorithm provides a global, multi-scale view of graph structure
that is robust to noise and structural perturbations. For the pairwise
case, complexity is dominated by the final assignment step (typically
O(n³)).

Key parameters:

- `ncomp`: Number of eigenvectors used (dimension of spectral basis)

- `q_descriptors`: Number of descriptor functions for multi-scale
  analysis

- `sigma`: Diffusion parameter controlling descriptor time scales

- `lambda`: Regularization balancing eigen-structure vs descriptor
  alignment

## References

Bogo, F., Romero, J., Loper, M., & Black, M. J. (2016). FAUST: Dataset
and evaluation for 3D mesh registration. In Proceedings of the IEEE
Conference on Computer Vision and Pattern Recognition (pp. 3794-3801).

Ovsjanikov, M., Ben-Chen, M., Solomon, J., Butscher, A., & Guibas, L.
(2012). Functional maps: a flexible representation of maps between
shapes. ACM Transactions on Graphics, 31(4), 1-11.

## See also

`grasp.hyperdesign`

## Examples

``` r
# \donttest{
# Example with hyperdesign graph data

# Create synthetic graph domains
set.seed(123)
domain1 <- list(
  x = matrix(rnorm(100), 50, 2),  # Node features
  design = data.frame(node_id = 1:50)
)
domain2 <- list(
  x = matrix(rnorm(100), 50, 2),
  design = data.frame(node_id = 1:50)  
)

# Create hyperdesign
hd <- structure(list(domain1 = domain1, domain2 = domain2), class = "hyperdesign")

# Run GRASP alignment with default parameters
result <- grasp(hd, ncomp = 20, q_descriptors = 50)

# Access alignment results
node_assignment <- result$assignment
aligned_embeddings <- result$s

# Use different parameters for robustness
result_robust <- grasp(hd, ncomp = 30, sigma = 1.0, lambda = 0.2)
# }

# \donttest{
# Example with hyperdesign graph data
library(multidesign)
library(tibble)

# Create synthetic graph domains (node coordinates for graph construction)
set.seed(123)
X1 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 1
X2 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 2

# Create design data frames with node indices (GRASP doesn't use features)
design1 <- tibble(node_id = 1:50)
design2 <- tibble(node_id = 1:50)

# Create multidesign objects
md1 <- multidesign(X1, design1)
md2 <- multidesign(X2, design2)

# Create hyperdesign from multidesign objects
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Run GRASP alignment (purely spectral, uses graph structure from X)
result <- grasp(hd, ncomp = 20, q_descriptors = 50)
# }
```
