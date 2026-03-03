# Sanitize labels before computing similarity matrices

Ensures there are no NA values passed into simfun by replacing each
missing label with a unique placeholder. Using unique placeholders
prevents all unlabeled samples from being treated as the same class
while still keeping the input compatible with functions like
neighborweights::binary_label_matrix, which cannot accept NA values.

## Usage

``` r
sanitize_labels_for_similarity(labels)
```

## Arguments

- labels:

  Vector of labels that may contain NA entries

## Value

Vector with the same length and names as \`labels\` but with all missing
entries replaced by unique placeholders.
