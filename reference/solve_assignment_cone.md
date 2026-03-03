# Solve Assignment Problem for CONE-Align

Solve Assignment Problem for CONE-Align

## Usage

``` r
solve_assignment_cone(Z1, Z2, Q, solver, anchors = NULL)
```

## Arguments

- Z1:

  First embedding matrix

- Z2:

  Second embedding matrix

- Q:

  Orthogonal transformation matrix

- solver:

  Assignment algorithm ("linear" or "auction")

## Value

Vector of assignment indices
