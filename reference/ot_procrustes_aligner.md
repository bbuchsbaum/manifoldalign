# OT + Procrustes aligner adapter for align_many()

Provides a pairwise interface implementing the orthogonal-invariant OT
objective of Alvarez-Melis et al. (AISTATS 2019), specialized to
learning an orthogonal map between feature spaces. Construct an
OT-Procrustes aligner descriptor

## Usage

``` r
ot_procrustes_aligner()
```

## Value

An object of class \`c("ot_procrustes_aligner", "aligner")\`.

## Details

This aligner estimates an orthogonal transform between two domains by
alternating between an entropic OT coupling (Sinkhorn) and a Procrustes
update (SVD) under the current coupling. It is intended as a
lightweight, initialization-robust alternative to adversarial alignment
for feature spaces that are identifiable only up to global orthogonal
transforms.

## Examples

``` r
algo <- ot_procrustes_aligner()
aligner_capabilities(algo)
#> $group
#> [1] "O"
#> 
#> $supports_multi
#> [1] FALSE
#> 
```
