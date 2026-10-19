# Refine a RAPID-MA fit with guarded class-balanced pseudo-labels

Refine a RAPID-MA fit with guarded class-balanced pseudo-labels

## Usage

``` r
rapid_ma_semisup(
  fit,
  labels = NULL,
  train_mask = NULL,
  mode = c("transductive", "inductive"),
  initial_probabilities = NULL,
  validation_probabilities = NULL,
  control = rapid_ma_semisup_control()
)
```

## Arguments

- fit:

  A fitted \`rapid_ma\` object.

- labels:

  Optional training-label vectors. In inductive mode, labels outside
  \`train_mask\` are discarded before any statistic is computed.

- train_mask:

  Optional logical vectors identifying training rows.

- mode:

  \`"transductive"\` or \`"inductive"\`.

- initial_probabilities:

  Optional named-class probability matrices, one per domain. This
  permits a supervised downstream classifier to seed the sparse relation
  teacher without exposing any additional labels.

- validation_probabilities:

  Optional out-of-fold probability matrices used only when
  \`control\$relation_validation = TRUE\`. Values at observed label rows
  must be produced without training on those rows.

- control:

  Settings from \[rapid_ma_semisup_control()\].

## Value

The updated \`rapid_ma\` object with a \`semisup\` provenance record.
