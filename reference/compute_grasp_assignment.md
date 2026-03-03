# Blocks D & E: Functional Map and Assignment - Corrected

Blocks D & E: Functional Map and Assignment - Corrected

## Usage

``` r
compute_grasp_assignment(
  basis1,
  basis2,
  desc1,
  desc2,
  M,
  distance_method = "cosine",
  solver_method = "linear",
  alpha = 0.85
)
```

## Arguments

- basis1:

  First domain's spectral basis

- basis2:

  Second domain's spectral basis

- desc1:

  First domain's descriptors

- desc2:

  Second domain's descriptors

- M:

  Rotation matrix from base alignment

- distance_method:

  Distance method: "cosine" or "euclidean"

- solver_method:

  Assignment solver: "linear", "hungarian", or "auction"

## Value

List with assignment and related matrices
