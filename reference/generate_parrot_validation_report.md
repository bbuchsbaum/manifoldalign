# Generate PARROT Validation Report

Creates a summary report from validation results

## Usage

``` r
generate_parrot_validation_report(results_df, output_dir, timestamp)
```

## Arguments

- results_df:

  Data frame from run_parrot_validation_suite

- output_dir:

  Directory to save report

- timestamp:

  Timestamp for file naming

## Value

A character string containing the formatted validation report, invisibly
