# Alignment method lifecycle policy

This policy controls when `manifoldalign` may describe an alignment method as
recommended, when a correctness finding must quarantine it, and when an older
method may be deprecated or removed. It deliberately separates numerical
correctness from benchmark performance. A method can be correct but inaccurate
for a particular regime; poor accuracy alone is not evidence of a wrong result.

The executable policy is returned by `alignment_lifecycle_policy()`. Package and
release checks use `validate_alignment_lifecycle()` and
`assert_alignment_release_ready()` rather than interpreting prose labels.

## Lifecycle states

| Metadata token | Meaning | May be recommended? | Public API? |
|---|---|---:|---:|
| `stable` | Correctness oracles and comparative accuracy/scale evidence pass. | yes | yes |
| `experimental` | The contract is explicit and smoke-tested, but stable evidence is incomplete. | no | yes |
| `legacy_reference` | Retained for comparison or reproducibility, not current guidance. | no | yes, if safe |
| `deprecated` | Replacement and removal schedule are published. | no | warning or disabled stub |
| `removed` | Public implementation is absent; only archival metadata may remain. | no | no |

Lifecycle transitions are conservative. Stable methods may move to
experimental, legacy/reference, or deprecated. Experimental methods may become
stable, legacy/reference, or deprecated. Legacy/reference methods must return
through experimental before stable. Deprecated methods may be repaired back to
experimental or removed. Removed is terminal. Repeating the current state is a
metadata update, not a transition.

## Evidence dimensions

Correctness, accuracy, and scalability are independent fields:

| Dimension | Ordered levels |
|---|---|
| correctness | `none`, `smoke`, `invariant`, `oracle` |
| accuracy | `none`, `single_fixture`, `comparative`, `release` |
| scalability | `none`, `smoke`, `comparative`, `release` |

Stable requires `oracle`, `comparative`, and `comparative`, respectively. It
also requires named oracle and benchmark artifacts plus reviewer and review-date
provenance. An oracle must check the mathematical quantity claimed by the
method: for example, an exact tiny solver, a differential implementation, known
permutation recovery, coupling marginals, objective agreement, or a numerical
invariant with a justified tolerance. A green smoke test is not an oracle.

Benchmark evidence must identify the applicable regime, fixture generation,
seed, decode rule, accuracy metric, runtime/memory metric, dependency/backend,
and artifact version. Accuracy evidence never substitutes for correctness
evidence, and correctness evidence does not establish comparative accuracy.

## Immediate quarantine

The following findings require `quarantined = TRUE` immediately:

- `fabricated_validation`: the validator does not compute the quantity it claims;
- `failed_numerical_invariant`: a material residual, feasibility, orthogonality,
  marginal, or equivalent invariant fails;
- `silent_unsupported_input`: an advertised input is ignored or accepted without
  implementing its semantics;
- `wrong_result_risk`: there is credible evidence that supported calls can
  return materially wrong results.

`low_accuracy` is explicitly not a quarantine trigger. It starts a benchmark or
dominance review.

An active quarantine record must identify its triggers, reason, repair deadline,
and repair plan. Fabricated validation, failed numerical invariants, and general
wrong-result risks disable the affected public API until repaired. Silent
unsupported inputs may retain the supported path only after the unsupported path
fails explicitly and is no longer advertised.

Restoration from quarantine to stable requires all triggers to be resolved,
fresh oracle and comparative benchmark artifact IDs, and an independent reviewer
with an ISO review date. Merely deleting a failing test or loosening its tolerance
does not count as repair evidence.

## Deprecation and removal

A deprecated or removed record must contain:

- `reason`;
- a distinct `replacement`;
- `warning_since`, `warning_through`, and `removal_at` package versions;
- a concrete `migration` path;
- an `archive` location for reproduction.

The warning and migration path must span at least one minor release.
`removal_at` must be the next major version after `warning_since`. Help pages and
runtime warnings must name both the replacement and removal version. An API that
fabricates evidence is disabled immediately, but its erroring compatibility stub
and migration metadata may remain through the normal warning window.

Deprecation solely for comparative performance requires reproducible evidence of
material domination in at least three applicable regimes and confirmation that
the method has no unique supported capability. A unique capability blocks a
performance-only deprecation; the method may instead become legacy/reference
while its niche is documented. Correctness quarantine and deprecation do not
need a dominance finding. An implementation that fails explicitly and lacks a
coherent supported contract may use the `incoherent_implementation` reason; this
does not turn an overt solver error into a silent wrong-result finding.

Archived scripts or an internal legacy path should preserve published analyses,
but they must not keep a misleading recommended public API alive.

## Record schema

Each registry record contains this lifecycle fragment:

```r
list(
  method = "method_name",
  state = "experimental",
  recommended = FALSE,
  quarantined = FALSE,
  api_enabled = TRUE,
  findings = character(),
  evidence = list(
    correctness = "invariant",
    accuracy = "comparative",
    scalability = "smoke",
    oracle_ids = "tests/testthat/test-method-invariants.R",
    benchmark_ids = "benchmarks/method-v1.csv",
    reviewed_by = "reviewer",
    reviewed_on = "2026-07-22"
  )
)
```

Deprecated methods add `deprecation`; quarantined methods add `quarantine`; a
stable restoration adds `restoration`. The executable validator reports typed
problem codes and the release assertion aggregates them across the registry.
Release checks fail if metadata is malformed, a forbidden transition is
requested, stable evidence is insufficient, an immediate trigger is not
quarantined, an unsafe quarantined API remains enabled, or any non-stable or
quarantined method is still marked recommended.

## Current package review

The evidence-backed decision for every exported aligner is recorded in
[`ALIGNMENT_METHOD_SUNSET_REVIEW.md`](ALIGNMENT_METHOD_SUNSET_REVIEW.md).
GlobalGeo, CONE-Align, and CONE-Align Multiple are currently legacy/reference;
`multiscale_manifold_align()` is deprecated; all other retained methods remain
experimental rather than recommended. The review cites completion denominators,
oracles, unique capabilities, and repair gates separately from accuracy.
