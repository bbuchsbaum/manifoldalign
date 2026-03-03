# Benchmark Data for Alignment Methods

A small, deterministic multi-domain dataset used across alignment
examples and tests. It contains three domains that share two latent
classes but exhibit different linear transformations and noise levels.
The object stores the raw matrices/design frames and the class labels so
alignment methods can generate comparable hyperdesigns on demand.

## Usage

``` r
alignment_benchmark
```

## Format

A list with the following components:

- domains:

  List of three elements named \`domain1\`, \`domain2\`, and
  \`domain3\`. Each element contains a centered numeric matrix \`x\`
  (samples x features) and a \`design\` data frame with \`id\`,
  \`condition\`, and \`domain\` columns.

- labels:

  Factor of length 80 giving the latent class for each sample.

## Examples

``` r
data <- alignment_benchmark
lapply(data$domains, function(dom) dim(dom$x))
#> $domain1
#> [1] 80  4
#> 
#> $domain2
#> [1] 80  4
#> 
#> $domain3
#> [1] 80  4
#> 
table(data$labels)
#> 
#> class_A class_B 
#>      40      40 
```
