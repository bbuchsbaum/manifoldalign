# Resolve a hyperdesign object into a normalized domain list

Guarantees a consistent list-of-lists with entries \`x\` and \`design\`
while preserving original names for downstream use.

## Usage

``` r
resolve_hyperdesign(data)
```

## Arguments

- data:

  Hyperdesign object produced by multidesign::hyperdesign() or a
  compatible structure used in tests.

## Value

List with \`domains\`, \`domain_names\`, and \`n_domains\`.
