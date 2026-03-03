# Apply UOT coupling as a barycentric map (dense, matrix signal)

Matrix version for time series or multi-contrast signals.

## Usage

``` r
uot_apply_map_dense_mat_cpp(
  cost,
  alpha,
  beta,
  fbar,
  gbar,
  epsilon,
  signal,
  delta = 1e-08
)
```

## Arguments

- cost:

  Dense cost matrix (n x m).

- alpha:

  Source masses (length n).

- beta:

  Target masses (length m).

- fbar:

  Translation-invariant source potential (length n).

- gbar:

  Translation-invariant target potential (length m).

- epsilon:

  Entropic regularization parameter (\> 0).

- signal:

  Source signal matrix (T x n).

- delta:

  Stabilizer added to the denominator.

## Value

A numeric matrix (T x m) in template space.
