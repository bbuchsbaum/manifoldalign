# Compute log second marginal and transported feature means (dense)

Given translation-invariant dual potentials (fbar,gbar), compute the log
of the second marginal (column sums) of the entropic UOT coupling and
the feature barycentric projection means per target node.

## Usage

``` r
uot_kl_logpi2_meanF_dense_cpp(cost, alpha, beta, fbar, gbar, epsilon, F = NULL)
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

- F:

  Optional source features (n x D). If provided, returns mean features.

## Value

A list with \`log_pi2\` (length m) and optionally \`mean_F\` (m x D).
