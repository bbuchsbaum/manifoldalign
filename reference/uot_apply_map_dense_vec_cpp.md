# Apply UOT coupling as a barycentric map (dense, vector signal)

Computes \\\hat s_j = (\sum_i \pi\_{ij} s_i) / (\pi\_{2,j} + \delta)\\
without materializing \\\pi\\, using translation-invariant potentials.

## Usage

``` r
uot_apply_map_dense_vec_cpp(
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

  Source signal vector (length n).

- delta:

  Stabilizer added to the denominator.

## Value

A numeric vector (length m) in template space.
