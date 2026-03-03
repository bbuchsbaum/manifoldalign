# Block indices for data list

Block indices for data list

## Usage

``` r
block_indices(data_list)
```

## Arguments

- data_list:

  List of data blocks

## Value

Matrix with start and end indices for each block

## Examples

``` r
data1 <- list(x = matrix(rnorm(20), 10, 2))
data2 <- list(x = matrix(rnorm(30), 15, 2))
data_list <- list(data1, data2)
indices <- block_indices(data_list)
indices
#>      start end
#> [1,]     1  10
#> [2,]    11  25
```
