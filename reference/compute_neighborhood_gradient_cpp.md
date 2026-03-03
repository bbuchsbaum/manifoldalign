# Compute neighborhood consistency gradient

Gradient of KL divergence term for neighborhood consistency

## Usage

``` r
compute_neighborhood_gradient_cpp(S_prev, W1, W2, eps = 1e-16)
```

## Arguments

- S_prev:

  Previous transport plan

- W1:

  Transition matrix of network 1

- W2:

  Transition matrix of network 2

- eps:

  Small epsilon for numerical stability

## Value

Gradient matrix
