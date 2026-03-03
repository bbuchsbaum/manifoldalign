# Construct a lightweight hyperdesign from matrices and labels

Builds a minimal 'hyperdesign' object (list of domains with fields \`x\`
and \`design\`) from a list of matrices and optional label vectors. This
avoids requiring the multidesign package for simple use cases.

## Usage

``` r
as_hyperdesign(
  X_list,
  labels = NULL,
  label_name = "label",
  domain_names = NULL
)
```

## Arguments

- X_list:

  List of matrices (samples x features) per domain

- labels:

  Optional list of label vectors aligned to rows of each matrix

- label_name:

  Name for the label column in each domain's design

- domain_names:

  Optional names for the domains

## Value

An object of class 'hyperdesign' suitable for
kema()/generalized_procrustes()

## Examples

``` r
# Create hyperdesign from list of matrices
X1 <- matrix(rnorm(100), 50, 2)
X2 <- matrix(rnorm(100), 50, 2)
X3 <- matrix(rnorm(100), 50, 2)
hd <- as_hyperdesign(list(X1, X2, X3))

# With labels
labels <- list(
  rep(c("A", "B"), each = 25),
  rep(c("A", "B"), each = 25),
  rep(c("A", "B"), each = 25)
)
hd_labeled <- as_hyperdesign(list(X1, X2, X3), labels = labels)
```
