# Alignment-method lifecycle policy

Return and enforce the package-wide policy used to distinguish numerical
correctness from empirical performance. The policy defines five
lifecycle states, independent evidence levels for correctness, accuracy,
and scalability, immediate quarantine triggers, permitted state
transitions, and deprecation cadence.

## Usage

``` r
alignment_lifecycle_policy()

validate_alignment_lifecycle(
  record,
  previous_state = NULL,
  previous_quarantined = FALSE,
  policy = alignment_lifecycle_policy()
)

assert_alignment_release_ready(records, policy = alignment_lifecycle_policy())
```

## Arguments

- record:

  A named list containing lifecycle metadata for one method.

- previous_state:

  Optional prior lifecycle state. When supplied, the transition to
  \`record\$state\` is validated.

- previous_quarantined:

  Whether the method is being restored from an active quarantine.
  Restoration to stable requires independent review and explicit oracle
  and benchmark evidence.

- policy:

  A policy object returned by \`alignment_lifecycle_policy()\`.

- records:

  One lifecycle record or a list of lifecycle records.

## Value

\`alignment_lifecycle_policy()\` returns a named policy list.
\`validate_alignment_lifecycle()\` returns an
\`alignment_lifecycle_problems\` data frame (zero rows means valid).
\`assert_alignment_release_ready()\` invisibly returns \`TRUE\`, or
raises an error containing the combined problems in its \`problems\`
field.

## Details

\`validate_alignment_lifecycle()\` validates one method record. A record
is a named list with \`method\`, \`state\`, \`recommended\`,
\`quarantined\`, \`api_enabled\`, \`findings\`, and \`evidence\` fields.
Stable records require oracle-level correctness plus comparative
accuracy and scalability evidence. Deprecated and removed records
additionally require a \`deprecation\` list; see
\`docs/ALIGNMENT_LIFECYCLE.md\` for the complete schema.

\`assert_alignment_release_ready()\` validates a list of records and
raises an \`alignment_release_error\` when any lifecycle invariant is
violated. In particular, deprecated, removed, experimental,
legacy/reference, or quarantined methods cannot be presented as
recommended.

## Examples

``` r
record <- list(
  method = "example_aligner",
  state = "stable",
  recommended = TRUE,
  quarantined = FALSE,
  api_enabled = TRUE,
  findings = character(),
  evidence = list(
    correctness = "oracle",
    accuracy = "comparative",
    scalability = "comparative",
    oracle_ids = "tests/testthat/test-example-oracle.R",
    benchmark_ids = "benchmarks/example-v1.csv",
    reviewed_by = "maintainer",
    reviewed_on = "2026-07-22"
  )
)
validate_alignment_lifecycle(record)
#> [1] method   code     severity field    message 
#> <0 rows> (or 0-length row.names)
assert_alignment_release_ready(record)
```
