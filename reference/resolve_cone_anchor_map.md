# Resolve anchor correspondences for CONE-Align

Returns a mapping from source indices to target indices using label
vectors where shared non-NA values indicate anchored correspondences.

## Usage

``` r
resolve_cone_anchor_map(vec1, vec2)
```

## Arguments

- vec1:

  Anchor labels for source domain (length n1)

- vec2:

  Anchor labels for target domain (length n2)

## Value

A list with \`map\` (integer vector length n1, NA for unanchored) and
\`values\` of shared anchors; or NULL if no anchors found.
