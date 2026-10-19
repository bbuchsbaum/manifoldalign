# Withdrawn KEMA Out-of-Sample Reconstruction Validator

This compatibility shim is retained temporarily so existing callers
receive an explicit failure instead of a fabricated validation score.
The former implementation did not fit KEMA or evaluate held-out
predictions; it returned \`expected_error\` plus random noise and has
therefore been withdrawn.

## Usage

``` r
validate_out_of_sample_reconstruction(
  expected_error = 0.14,
  tolerance = 0.05,
  test_fraction = 0.2
)
```

## Arguments

- expected_error:

  Deprecated and ignored.

- tolerance:

  Deprecated and ignored.

- test_fraction:

  Deprecated and ignored.

## Value

This function emits a deprecation warning and then errors.

## Details

No direct reconstruction replacement is currently available because
KEMA's out-of-sample projection has not passed an independent
reconstruction oracle. Use \[cv_alignment_rows()\] for leakage-safe
held-out alignment scoring. Historical results from this function are
invalid.

## See also

\[cv_alignment_rows()\]

## Examples

``` r
if (FALSE) { # \dontrun{
validate_out_of_sample_reconstruction()
} # }
```
