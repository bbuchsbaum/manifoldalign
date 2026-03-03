# Capabilities for an aligner

Capabilities for an aligner

## Usage

``` r
aligner_capabilities(algo)
```

## Arguments

- algo:

  An aligner object created by new_aligner or an aligner constructor

## Value

A list with elements group (character) and supports_multi (logical)

## Examples

``` r
algo <- new_aligner("test", group = "O")
caps <- aligner_capabilities(algo)
```
