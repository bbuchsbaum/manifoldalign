# Alignment tournament v2: interpretation and decisions

This document interprets the frozen final tournament. It is subordinate to the
machine-readable artifacts in this directory: [`raw_results.csv`](raw_results.csv),
[`controls.csv`](controls.csv), [`summary_by_regime.csv`](summary_by_regime.csv),
and [`validation_problems.csv`](validation_problems.csv). The final manifest has
MD5 `2f79e4b34867edba288239fbc427a437`; all 588 method-case rows are present and
the validation-problem table is empty.

## Bottom line

RAPID-MA is the current package choice when the task truly requires all of the
following: caller-supplied sparse topology, nonnegative edge weights, spatial
positions, mixed node attributes, unequal graph sizes, and partial overlap. It
completed all 42 final cases and is the only comparator in this tournament that
consumed the supplied relations. That is an input-contract and robustness win,
not a universal accuracy win.

RAPID-MA is not currently the best closed-world matcher, and its matching
accuracy is not yet satisfactory. Spectral MNN won the tiny exact oracle
regime. More importantly, a position-only control beat RAPID-MA in every
non-oracle regime, often by a wide margin. The current fusion stage therefore
destroys useful anchor-calibrated spatial signal on these cases.

For downstream classification, the conclusion is more favorable. RAPID-MA's
frozen source-label probe beat the no-alignment source-only classifier in every
non-oracle regime. It still remained far below the target-only ceiling, and no
method recovered unseen/private classes.

## Matching accuracy and completion

All values below are intent-to-run Hits@1: explicit errors count as zero and
are never dropped. `ok/intended` is the completion denominator.

| Regime | Leading applicable methods |
|---|---|
| Oracle, closed | Spectral MNN 0.750 (6/6); CONE 0.583 (6/6); RAPID-MA 0.472 (6/6); GPCA and GW 0.403 (6/6 each) |
| Small, closed | RAPID-MA 0.214 (6/6); CONE 0.168 (4/6); LRA 0.074 (6/6); PARROT 0.069 (6/6); Spectral MNN 0.067 (6/6) |
| Small, partial | RAPID-MA 0.096 (9/9); GlobalGeo 0.009 (2/9); FPGW 0.0046 (9/9) |
| Medium, closed | RAPID-MA 0.0265 (6/6); Spectral MNN 0.0117 (6/6); Token-OT 0.0052 (6/6); SSMA 0.0025 (6/6) |
| Medium, partial | RAPID-MA 0.0258 (9/9); no other method was semantically and resource applicable |
| Large, partial | RAPID-MA 0.00433 (6/6); no other method was semantically and resource applicable |

The exact oracle is deliberately tiny and should not be mistaken for a
realistic accuracy ranking. The non-oracle rows show a different problem:
RAPID-MA is usually the strongest package method among those allowed to run,
but absolute correspondence accuracy falls rapidly with graph size and
perturbation severity.

## Negative controls expose negative transfer

The position-only control estimates a rigid transform from training anchors
and matches positions without observing topology or attributes. Its Hits@1 is:

| Regime | Position only | RAPID-MA |
|---|---:|---:|
| Small, closed | 0.7014 | 0.2139 |
| Small, partial | 0.7330 | 0.0963 |
| Medium, closed | 0.0642 | 0.0265 |
| Medium, partial | 0.0551 | 0.0258 |
| Large, partial | 0.0309 | 0.0043 |

This is the tournament's most actionable result. A production RAPID-MA update
should never underperform the calibrated position view when positions and
anchors are available. The next algorithm version should use development-only
anchor holdouts to gate view weights, retain a position-first candidate path,
and fall back to the position solution when multi-view fusion does not improve
the held-out anchor loss. The current final cases must not be used to tune that
revision; a new versioned final receipt is required.

Raw-attribute nearest-neighbour and topology-descriptor controls were weak at
medium and large scale. The spatial control is therefore not a disguised
all-purpose identity leak; it isolates a genuinely strong view in this
generator family.

## Classification transfer

Classification is a common frozen source-label downstream probe applied to
each alignment output, not a method-specific private classifier. Intent-to-run
overall accuracy is compared with source-only transfer and the target-only
ceiling:

| Regime | RAPID-MA | Source only | Target-only ceiling |
|---|---:|---:|---:|
| Small, closed | 0.6019 | 0.1829 | 0.9954 |
| Small, partial | 0.4458 | 0.1620 | 0.9784 |
| Medium, closed | 0.5073 | 0.1844 | 0.9878 |
| Medium, partial | 0.4563 | 0.2211 | 0.9931 |
| Large, partial | 0.5059 | 0.3102 | 0.9961 |

RAPID-MA adds useful cross-domain transfer in every non-oracle regime. It is
not always the strongest closed-world classifier: PARROT reached 0.8588 and
LRA 0.7060 on small closed cases; Token-OT reached 0.6701 and SSMA 0.6635 on
medium closed cases. These task-specific wins justify retaining those methods
even when their node-correspondence scores are poor. Target-only performance
near 0.98--1.00 also shows substantial remaining headroom. Unseen-class
accuracy is zero throughout, so none of these results supports open-set label
discovery claims.

## Reliability, efficiency, and scope

RAPID-MA completed 42/42 intended cases. On the six large cases (11,300 or
22,600 total nodes), runtime was about 14.4--36.8 seconds, sampled peak RSS was
about 0.94--2.07 GB, and retained candidate support was 352,000 or 704,000
edges. This validates the current tournament implementation through 22,600
total nodes with full scoring, not through one million nodes.

The older scaling receipt records a completed 220,000-node run and a timed-out
1,000,002-node attempt. Those historical facts remain useful, but the
million-node gate is not passed. Likewise, the spatially blocked inductive
receipt failed its predeclared improvement gate. This tournament is
transductive alignment followed by a frozen source-label probe; it does not
replace a blocked out-of-sample evaluation.

GlobalGeo completed 8/21 intended cases. Thirteen cases failed explicitly
because the landmark graph was disconnected, including all 6 small
closed-world cases and 7 of 9 small partial-overlap cases. CONE completed 10/12
intended cases; two equal-size closed-world cases failed because its Procrustes
inputs had unequal latent dimensions. These failures are visible in the raw
rows and are not evidence of silent wrong answers, but the corresponding
capability claims should be narrowed until repaired.

FPGW completed 21/21 applicable cases and retained its exact fixed-mass oracle
status, but small-case median runtime was about 12--13 seconds and matching was
weak on this attributed-graph task. It remains a coherent metric-measure OT
method; it should not be described as consuming observed graph topology.

## Lifecycle consequences

- Keep RAPID-MA experimental. It has the broadest observed-graph contract and
  complete tournament execution, but matching, blocked-inductive, and
  million-node promotion gates remain open.
- Keep Spectral MNN, SSMA, Token-OT, LRA, GPCA, GRASP, and PARROT as
  experimental feature-derived or feature-manifold comparators in their stated
  regimes. PARROT and LRA retain distinctive classification evidence.
- Keep GW and FPGW as metric-measure OT methods. Their low matching accuracy on
  this task is not a correctness defect.
- Reclassify GlobalGeo as a legacy/reference method until its disconnected
  landmark-graph contract is repaired; withdraw partial-overlap support claims.
- Reclassify CONE as legacy/reference until latent-dimension reconciliation
  passes the frozen cases. Do not deprecate it solely for being slower or less
  accurate: it was the second-best oracle matcher and second-best small closed
  matcher on an intent-to-run basis.
- Keep KEMA experimental after its numerical-fidelity repair; deprecate only
  the duplicate `kema_orig()` entry point as already scheduled.
- Keep Linear Similarity Embedding as a single-domain target-similarity
  estimator. It is correctly not applicable to this cross-domain tournament.
- Keep `multiscale_manifold_align()` deprecated. Its prior solver and
  completion defects remain the package's method-level incoherence case.

The package-wide sunset review is the authoritative lifecycle decision record;
this section records the tournament evidence feeding that review.

## Reproduction and audit

Run from the package root:

```sh
Rscript inst/benchmarks/run_alignment_tournament.R \
  --profile=final \
  --output=docs/evidence/alignment-tournament-v2 \
  --resume=true
```

The runner isolates each method-case pair in an external process on Unix,
enforces predeclared resource envelopes, checkpoints every typed row, and
records peak RSS, parameters, provenance, session information, and recursive
checksums for every published receipt file. Resumable binary checkpoints and
the debug log remain local operational state and are deliberately outside the
published checksum manifest. An official receipt requires 588 rows, the frozen
manifest hash, and zero rows in `validation_problems.csv`.
