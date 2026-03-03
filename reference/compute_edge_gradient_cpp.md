# Compute edge consistency gradient

Efficient computation of the gradient of edge consistency regularizer

## Usage

``` r
compute_edge_gradient_cpp(S, A1, A2, X1, X2)
```

## Arguments

- S:

  Current transport plan (n1 x n2)

- A1:

  Adjacency matrix of network 1 (sparse)

- A2:

  Adjacency matrix of network 2 (sparse)

- X1:

  Features of network 1 (n1 x d)

- X2:

  Features of network 2 (n2 x d)

## Value

Gradient matrix (n1 x n2)
