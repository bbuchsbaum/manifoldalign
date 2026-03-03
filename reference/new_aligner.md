# Generic alignment adapter interface and transform helpers

This file defines a minimal, method-agnostic interface for lifting
pairwise manifold alignment algorithms to the multi-domain setting.
Adapters for concrete methods (e.g., GRASP, PARROT, GPA, KEMA) should
implement these generics for their aligner classes. Create a simple
aligner descriptor

## Usage

``` r
new_aligner(name, group, supports_multi = FALSE, ...)
```

## Arguments

- name:

  Character name of the algorithm

- group:

  One of "O" (orthogonal), "GL" (linear), "perm" (assignment), or
  "kernel" (embedding-only)

- supports_multi:

  Logical; if TRUE, adapter should also implement a native multi-domain
  entry point via fit_many()

- ...:

  Additional fields stored on the object

## Value

An object of class c("\<name\>\_aligner", "aligner")

## Examples

``` r
algo <- new_aligner("my_method", group = "O")
```
