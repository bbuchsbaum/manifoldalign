# RAPID-MA bounded node-matching validation

Status: multi-view matcher improvement verified; exact correspondence remains
an open accuracy problem.

## Bottom line

On 20 untouched synthetic cross-modality seeds, the new RAPID-MA matcher
reached 14.13% Hits@1 at 100% coverage. The legacy prototype-bucket matcher
reached 3.50% Hits@1 at 65.60% coverage, and nearest-neighbor retrieval in the
RAPID embedding reached 1.57% Hits@1. The multi-view matcher therefore improves
overall exact retrieval by 4.04x over the legacy matcher and 9.02x over the
embedding-only path.

The result is materially better, but 14.13% is not solved node correspondence.
A position-only diagnostic reaches 18.53% Hits@1, showing that most of the
available structural signal is now used while a meaningful gap remains.

## Frozen results

The fixture has 300 source nodes with known mates among 360 target nodes,
different source and target feature widths, nonlinear target features, noisy
positions, mixed numeric/categorical attributes, missing attributes, class
imbalance, and only four labels per class. Development used seeds 7101--7103.
Implementation commit `8a1992b6fe7f8287b1444e79e5a23327298aa33a`
was then evaluated on untouched seeds 7201--7220.

Hits@1 below is normalized over all 300 source nodes; unmatched rows count as
incorrect. "Matched Hits@1" conditions on rows for which a method returned a
valid target.

| Method | Runs | Hits@1 | Matched Hits@1 | Coverage | Fit (s) | Match (s) |
|---|---:|---:|---:|---:|---:|---:|
| Position kNN diagnostic | 20/20 | **18.53%** | 18.53% | 100.00% | 0 | 0.0019 |
| RAPID multi-view, independent | 20/20 | **14.13%** | 14.13% | **100.00%** | 0.386 | 0.242 |
| RAPID multi-view, one-to-one | 20/20 | 14.02% | **15.93%** | 87.72% | 0.386 | 0.252 |
| RAPID prototype, one-to-one | 20/20 | 3.50% | 5.32% | 65.60% | 0.386 | 0.242 |
| RAPID embedding kNN | 20/20 | 1.57% | 1.57% | 100.00% | 0.386 | NA |
| Spectral-MNN embedding kNN | 20/20 | 0.58% | 0.58% | 100.00% | 0.044 | NA |
| Token-OT embedding kNN | 20/20 | 0.52% | 0.52% | 100.00% | 0.055 | NA |
| SSMA embedding kNN | 20/20 | 0.27% | 0.27% | 100.00% | 0.032 | NA |
| LeMA embedding kNN | 20/20 | 0.05% | 1.25% | 4.00% | 0.378 | NA |

LeMA receives the paired labeled-landmark correspondence expected by its
reference formulation, but its source embedding is valid for only those 12
landmarks. Its conditional Hits@1 is 1.25%; multiplying by 4% coverage gives
the 0.05% all-source result in the table.

The RAPID matcher returns only its selected target, so its reported MRR equals
Hits@1 (MRR@1). Embedding and position diagnostics retain bounded rankings up
to 100 candidates. Their MRR@100 values are available in the raw artifact and
must not be compared directly with RAPID's MRR@1 as though the rank caps were
the same.

## What changed

The legacy matcher generated candidates only from the inverted sparse
node-to-prototype index. Nodes in a broad class prototype frequently did not
share a candidate bucket with their true mate, even when their spatial
positions were close.

The new matcher unions bounded candidates from:

- the inverted sparse prototype index;
- the learned latent embedding;
- retained positions;
- retained encoded attributes;
- an optional retained structural signature.

Only this bounded union is reranked with normalized latent, structural,
position, attribute, and prototype costs. The default requests 12 neighbors
per comparable view and retains at most 32 candidates per source. The held-out
run retained exactly 9,600 candidate edges for 300 source rows, never a
300-by-360 dense pair matrix. The returned matcher object averaged 0.020 MB.

Candidate construction was also changed from repeated vector and bucket
growth to grouped linear-size assembly. Tests verify that doubling the fixture
keeps candidate edges bounded by `source_nodes * candidate_cap`, matcher state
grows by less than 3x, and no dense cross-domain allocation is reported.

Two assignment modes are explicit:

- `one_to_one` preserves the pair-transform behavior and prevents reuse of a
  target; competition for target nodes lowers coverage on rectangular data;
- `independent` provides retrieval semantics, permits target reuse, and gives
  full source coverage. It is the appropriate mode for Hits@1 retrieval.

Anchors remain fixed and protected in both modes. The one-to-one mode remains
the API default for backward compatibility.

## Interpretation and limits

The position-only diagnostic is not a learned alignment method. It establishes
how much correspondence information the fixture's noisy coordinates make
available to a bounded nearest-neighbor index. The 14.13% multi-view result is
76% of that diagnostic Hits@1 while also using alignment, prototypes, and
attributes.

This benchmark shows that explicit structure fixes a real candidate-recall
failure. It does not establish high exact-match accuracy on a real sensor pair,
graphs without comparable coordinates, partial-overlap identity matching, or
million-node matcher throughput. The classification and blocked-inductive
reports should not be read as evidence of exact correspondence.

## Reproduction and artifacts

```sh
RAPID_MA_PROFILE=full \
  Rscript inst/benchmarks/rapid_ma_matching.R docs/rapid-ma-matching
```

- `rapid-ma-matching-raw.csv`: every seed and method;
- `rapid-ma-matching-summary.csv`: means and standard deviations;
- `rapid-ma-matching-metadata.txt`: frozen seeds, commit, controls, platform,
  and source hashes.

SHA-256:

```text
886bcee27443559fff90d370368e3dc6e482405a67cd9b9d16d01d23af984bc0  rapid-ma-matching-raw.csv
d58322554702792449a3b442643f7c25b7564922d875aa7794194aeeb25c505e  rapid-ma-matching-summary.csv
0d2fe446841fb2cda5b21cf1ecbea6d1eb2168e8ed69934c762dbae1aa13b743  rapid-ma-matching-metadata.txt
```
