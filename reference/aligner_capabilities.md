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

A list with elements group (character) and supports_multi (logical). Use
\[alignment_method_info()\] for the complete capability and lifecycle
contract; registry meta-tests verify that these adapter fields agree
with that record.

## Examples

``` r
algo <- new_aligner("test", group = "O")
caps <- aligner_capabilities(algo)
```
