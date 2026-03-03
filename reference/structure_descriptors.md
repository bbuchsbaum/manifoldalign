# Structural graph descriptors (HKS, diffusion coordinates)

These helpers compute \*\*intrinsic\*\* node descriptors from a
within-domain spectral basis (eigenvalues/eigenvectors). They are
designed to be reusable across multiple alignment methods (e.g., GRASP,
spectral MNN, and as features for UOT costs).

The main entry points are: - \[hks_time_grid()\] to construct a
diffusion-time grid. - \[compute_hks_descriptors()\] to compute Heat
Kernel Signatures (HKS). - \[compute_diffusion_coordinates()\] to
compute diffusion coordinates at one or more times.

## Value

See individual functions.

## Details

\*\*Basis format\*\*

These functions accept either: - a list with fields \`values\` (length
K) and \`vectors\` (N x K), or - explicit \`values\`/\`vectors\`
arguments.

For Laplacian-based bases, the first eigenvalue is typically 0 and
corresponds to a constant mode. By default, eigenpairs with \`value \<=
eigen_tol\` are dropped.
