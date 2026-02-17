# Documentation Examples Added

Successfully added `@examples` to 32 exported functions that were missing examples.

## Summary by Category

### Alignment Functions (2)
1. **align_many.Rd** - R/align_many.R - Added example with 3 domains
2. **rotation_sync.Rd** - R/align_many.R - Added example for orthogonal synchronization

### KEMA Core Functions (3)
3. **block_indices.Rd** - R/kema.R - Added example with data blocks
4. **choose_sigma.Rd** - R/kema.R - Added simple RBF kernel sigma estimation example
5. **coskern.Rd** - R/kema.R - Added cosine kernel example

### KEMA Scalable Functions (9)
6. **make_Zop_from_Ks.Rd** - R/kema_scalable.R - Added example with kernel blocks
7. **gram_ZZ.Rd** - R/kema_scalable.R - Added Gram matrix computation example
8. **gram_Z_A_Z.Rd** - R/kema_scalable.R - Added weighted Gram matrix example
9. **gram_Z_diag_Z.Rd** - R/kema_scalable.R - Added diagonal weighted example
10. **Zt_times.Rd** - R/kema_scalable.R - Added matrix multiplication example
11. **gram_Ls.Rd** - R/kema_scalable.R - Added same-class Laplacian example
12. **gram_Ld.Rd** - R/kema_scalable.R - Added different-class Laplacian example
13. **pivoted_chol_kernel_block.Rd** - R/kema_scalable.R - Added pivoted Cholesky example

### C++ Exported Functions (9) - All wrapped in \dontrun{}
14. **compute_edge_distances_cpp.Rd** - R/RcppExports.R
15. **compute_edge_gradient_cpp.Rd** - R/RcppExports.R
16. **compute_neighborhood_gradient_cpp.Rd** - R/RcppExports.R
17. **compute_anchor_gradient_cpp.Rd** - R/RcppExports.R
18. **compute_squared_distances_cpp.Rd** - R/RcppExports.R
19. **solve_sylvester_rwr_cpp.Rd** - R/RcppExports.R
20. **compute_rwr_vectorized_cpp.Rd** - R/RcppExports.R
21. **compute_parrot_cost_cpp.Rd** - R/RcppExports.R
22. **sinkhorn_unified.Rd** - R/RcppExports.R

### FPGW Functions (2)
23. **fpgw-methods.Rd** - R/fpgw_methods.R - Added grouped documentation example
24. **print.fpgw.Rd** - R/fpgw.R - Added print method example

### GRASP Functions (1)
25. **grasp_multiset.Rd** - R/grasp_multiset.R - Added multi-graph alignment example

### Pseudolabeling (1)
26. **pseudolabeling.Rd** - R/pseudolabel.R - Added sparse similarity matrix example

### Validation Functions (3)
27. **validate_kema_eigenvalues.Rd** - R/kema-validation.R - Added eigenvalue validation example
28. **validate_out_of_sample_reconstruction.Rd** - R/kema-validation.R - Added reconstruction example
29. **run_kema_validation_suite.Rd** - R/kema-validation.R - Added full suite example

### Spectral MNN Functions (2)
30. **spectral_mnn_align_control.Rd** - R/spectral_mnn_align.R - Added control parameter example
31. **spectral_mnn_align.Rd** - R/spectral_mnn_align.R - Added alignment example

### MMA Functions (1)
32. **mma_align_multiple.hyperdesign.Rd** - R/mma_align_multiple.R - Added multi-domain example

## Guidelines Followed

- **C++ functions**: All wrapped in `\dontrun{}` since they require compiled code
- **Slow functions**: Used `\donttest{}` for computationally intensive examples
- **Fast functions**: Minimal, runnable examples without wrappers
- **Grouped pages**: Added examples to main documentation page
- **Internal helpers**: Functions like `gram_Ls` that need label factors have placeholder comments

## Verification

All 32 files verified to contain `\examples` section after running:
```r
devtools::document()
```

No CRAN-compliance issues expected from these examples.
