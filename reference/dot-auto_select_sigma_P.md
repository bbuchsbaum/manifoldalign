# Auto-select sigma_P using histogram-spread heuristic

Auto-select sigma_P using histogram-spread heuristic

## Usage

``` r
.auto_select_sigma_P(X, W0, T, M, verbose = FALSE)
```

## Arguments

- X:

  Preprocessed data matrix

- W0:

  Initial weight matrix

- T:

  Target similarity matrix

- M:

  Mask matrix

- verbose:

  Logical for progress reporting

## Value

Numeric scalar giving the optimal sigma_P value
