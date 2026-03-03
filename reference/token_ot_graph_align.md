# Token-level Optimal Transport Graph Alignment

Align two graphs (node sets) by representing each node as a \*\*bag of
context tokens\*\* and aligning nodes via: 1) multi-view node embeddings
(late-fusion similarity), and 2) token-level entropic optimal transport
(Sinkhorn) between node bags.

## Usage

``` r
token_ot_graph_align(data, ...)

# S3 method for class 'hyperdesign'
token_ot_graph_align(
  data,
  anchors = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = token_ot_graph_align_control(),
  ...
)

# S3 method for class 'list'
token_ot_graph_align(data, ...)
```

## Arguments

- data:

  A \`hyperdesign\` object (2 domains) or a list of two matrices.

- ...:

  Additional arguments (reserved for future extensions).

- anchors:

  Optional unquoted column name in each domain's design table providing
  anchor identifiers (shared values indicate known correspondences). Use
  \`NULL\` for unsupervised alignment.

- preproc:

  Optional preprocessing (default: \`multivarious::center()\`).

- ncomp:

  Integer number of components stored in \`s\` (scores) and \`v\`
  (loadings). Used for the baseline embedding view when available.

- control:

  Control settings from \[token_ot_graph_align_control()\].

## Value

A \`multiblock_biprojector\` with subclass \`token_ot_graph_align\`.
Extras include \`transport_plan\` (sparse coupling), \`assignment\`
(greedy or Hungarian), and \`diagnostics\`.

## Details

This implementation is designed to stay sparse and efficient by
generating a top-k candidate set per source node (ANN / kNN) and running
Sinkhorn only on the candidate bipartite graph.
