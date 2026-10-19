# Testing strategy for `manifoldalign`

`manifoldalign` treats numerical correctness, comparative accuracy, and
scalability as separate claims. A green package test suite establishes only the
claims represented by its oracles and invariants; it does not establish that a
method wins an accuracy benchmark.

The machine-readable source of truth is
`tests/testthat/ALIGNMENT_TEST_MATRIX.yml`. Its `ci_contract` section is checked
against the canonical method registry on every PR. The contract names each
method's test families, PR and nightly files, mandatory tests, release evidence,
and known gaps.

## Evidence tiers

| Tier | Trigger | Required work | Evidence produced |
|---|---|---|---|
| PR | Every push and pull request | Deterministic contracts, invariants, metamorphic checks, adversarial cases, regressions, and tiny exact or differential oracles | Per-test and per-method CSVs; warnings and non-mandatory skips are visible |
| Nightly | Daily schedule or manual dispatch | PR tier plus gated benchmarks, 20-seed synthetic accuracy and matching comparisons, dependency-backed paths, and moderate scaling with sampled peak RSS | Test status plus regenerated accuracy, matching, scaling, metadata, and memory receipts |
| Release | Version tag or manual dispatch | Nightly tier plus long stress settings, full scaling and ablation evidence, and optional leakage-safe real-scene evaluation | Full artifact bundle; unavailable real data are recorded explicitly |

Run a tier locally from the package root:

```sh
Rscript tools/run-alignment-test-tier.R pr
Rscript tools/run-alignment-test-tier.R nightly
Rscript tools/run-alignment-test-tier.R release
```

Use `--output-dir=/path/to/evidence` to keep local artifacts outside the
checkout. Nightly and release runs are intentionally expensive.

## Mandatory tests

A mandatory test is addressed by exact `path::test title`. CI fails if it:

- fails or errors;
- is skipped, including a skip caused by a missing optional dependency; or
- disappears because its file or title was renamed without updating the
  contract.

Warnings do not currently fail the tier, but they are emitted with method and
severity so that recurring numerical warnings cannot hide in aggregate test
output. Experimental methods may have gaps, but each gap must be named. A
method cannot be marked stable unless it has all required computational
families, a mandatory oracle, and release evidence.

## Computational families

- `contract`: types, shapes, indexing, sparsity, and error behavior.
- `oracle`: a tiny exact result or an independently implemented differential
  reference.
- `invariant`: marginals, residuals, orthogonality, monotone objectives,
  bijections, or other necessary mathematical properties.
- `metamorphic`: permutation, relabeling, domain-order, translation, rotation,
  or scale behavior where the method should be equivariant or invariant.
- `adversarial`: disconnected graphs, isolates, duplicates, empty overlap,
  ill-conditioning, extreme weights, and unequal sizes.
- `regression`: a fixture for a previously observed defect.
- `approximation`: approximate versus exact agreement, monotone quality as the
  budget grows, and graceful degradation.
- `benchmark`: comparative predictive or matching evidence.
- `stress`: long-running scale, stability, and memory checks.

Stable methods require the first seven families. `benchmark` and `stress`
remain required by the lifecycle policy before a comparative accuracy or scale
claim can be promoted.

## Anti-tautology rule

Assertions such as `n_correct >= 0`, `length(x) >= 0`, or an accuracy lower
bound without a meaningful upper bound are forbidden. They cannot distinguish
a correct implementation from a broken one. The matrix names the forbidden
patterns, and the meta-test scans the suite for them.

## CI artifacts

Each tier writes under `artifacts/alignment-ci/<tier>/`:

- `alignment-contract-problems.csv`;
- `alignment-test-results.csv` with exact test, method, mandatory flag, status,
  severity, counts, and duration;
- `alignment-test-summary.csv` aggregated by method, severity, and status;
- `alignment-contract-status.csv` with known gaps and release availability;
- `alignment-test-summary.md` for the GitHub job summary.

Scheduled workflows also write benchmark raw data, summaries, and provenance
metadata below `benchmarks/`. Artifacts are uploaded even when a job fails.

## Frozen package tournament

The final attributed-graph tournament is a release artifact, not an ordinary
unit test. It runs 14 methods against all 42 frozen cases in external processes
and writes a complete 588-row method-case grid. Final cases have zero tuning
budget. Unsupported and not-applicable tasks remain visible; explicit errors
remain in intent-to-run denominators.

Regenerate it from the package root with:

```sh
Rscript inst/benchmarks/run_alignment_tournament.R \
  --profile=final \
  --output=docs/evidence/alignment-tournament-v2 \
  --resume=true
```

The official receipt is valid only when `validation_problems.csv` has zero rows
and every published entry in `checksums.csv` matches. Resumable binary
checkpoints and the debug log are local operational state rather than release
artifacts. The interpretation must report negative controls and failed
promotion gates, not only the strongest method row.

## Adding or promoting a method

1. Add or update its canonical registry record.
2. Add exactly one `ci_contract.method_records` entry.
3. Declare only families demonstrated by existing tests.
4. Make at least one dependency-independent oracle or invariant mandatory.
5. Add adversarial and regression cases for every discovered failure mode.
6. Name every remaining gap while the method is experimental.
7. Run the PR tier and inspect the method-level artifact, not just the exit
   status.
8. Before stable promotion, run nightly and release evidence and satisfy the
   lifecycle gates.

## Current limitations

- No method is stable or recommended yet.
- FPGW now has a mandatory base-R exhaustive-flow oracle covering compiled
  balanced and partial transport, an independent four-index objective,
  weighted unequal marginals, fixed mass, and target-permutation equivariance.
  Its larger `lpSolve` suite remains nightly-only. The former `lambda > 0`
  route is withdrawn; only classical FGW and explicit fixed-`rho` partial FPGW
  are supported. It completed 21/21 applicable final tournament cases, but is
  slow and has weak matching on this attributed-graph task; scale evidence is
  still missing.
- The attributed-graph v2 tournament is complete: 588/588 typed rows and zero
  validation problems. RAPID-MA completed 42/42 and has the broadest supplied-
  graph contract, but its node matching loses to the position-only control in
  every non-oracle regime. Position-gated fusion, blocked-inductive accuracy,
  and the million-node gate remain named gaps. The checked
  [interpretation](../docs/evidence/alignment-tournament-v2/INTERPRETATION.md)
  and [raw receipt](../docs/evidence/alignment-tournament-v2/raw_results.csv)
  are the source for these claims.
- GlobalGeo is legacy/reference after 8/21 completion and no longer claims
  partial-overlap support. Pairwise and multi-domain CONE are legacy/reference
  after the pairwise path completed 10/12 and exposed two latent-dimension
  failures. `multiscale_manifold_align()` remains deprecated.
- Chikusei real-scene data are external, about 1 GB, and CC BY 3.0. CI downloads
  them only after explicit acceptance in a manual release workflow. Otherwise
  the release artifact records real-scene evidence as not requested or not
  available.

See `docs/ORACLE_CI_POLICY.md` for the enforcement and promotion rules,
`docs/BENCHMARK_PROTOCOL.md` for the accuracy protocol, and
`docs/ALIGNMENT_METHOD_SUNSET_REVIEW.md` for the evidence-backed decision on
every exported aligner.

## Source audits and installed-package tests

`R CMD check` runs the runtime tests under `tests/testthat` against the installed
package. Checkout-specific audits of the CI contract, tracked comparison table,
and benchmark scripts live under `tools/tests`. The required
`Rscript tools/run-alignment-test-tier.R pr` command runs both suites, records
both in its evidence tables, and fails on either suite's errors. This separation
does not skip source audits when the package is checked; they remain a required
CI job alongside the installed-package check.
