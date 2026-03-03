# Solve Assignment for Multi-Graph GRASP

Computes node assignment between two embedded graphs using cosine or
Euclidean distance.

## Usage

``` r
solve_assignment_multiset(
  E1,
  E2,
  distance = "cosine",
  solver = "auction",
  alpha = 0.5
)
```

## Arguments

- E1:

  First embedding matrix

- E2:

  Second embedding matrix

- distance:

  Distance metric: "cosine" or "euclidean"

- solver:

  Assignment solver: "auction" or "linear"

- alpha:

  Weight on cosine vs. Euclidean distance when distance = "cosine"

## Value

Assignment vector
