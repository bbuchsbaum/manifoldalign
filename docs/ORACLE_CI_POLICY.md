# Computational oracle and CI policy

This policy prevents a method from acquiring a strong lifecycle label merely
because its code executes. It applies to every canonical alignment method in
`alignment_method_registry()`.

## Contract authority

Three artifacts have distinct responsibilities:

1. `R/alignment_registry.R` defines the method, input, lifecycle, and claimed
   evidence.
2. `tests/testthat/ALIGNMENT_TEST_MATRIX.yml` defines the computational test
   families and exact CI paths supporting those claims.
3. `docs/BENCHMARK_PROTOCOL.md` defines leakage-safe comparative accuracy and
   scalability evaluation.

`tests/testthat/test-alignment-oracle-ci.R` fails on drift between the first two.
`tools/run-alignment-test-tier.R` executes the relevant suite and emits the
observed results. A registry claim, a test declaration, and an observed pass
are all required; none substitutes for the others.

## Stable-method gate

A stable method must demonstrate:

- contracts and shapes;
- a tiny exact or independently implemented differential oracle;
- mathematical invariants appropriate to its output;
- applicable metamorphic behavior;
- adversarial boundary behavior;
- regression coverage for every known correctness defect; and
- approximation quality or graceful degradation where it uses an approximate
  solver.

It must also name at least one mandatory PR oracle and at least one release
evidence artifact. The lifecycle policy separately requires comparative
accuracy and scalability evidence. The meta-test mutation-checks this rule, so
it remains enforced even while the current registry contains no stable method.

## Mandatory-oracle semantics

Mandatory paths use `relative/file.R::exact test_that title` identifiers. The
PR tier must observe every identifier. Missing, failed, errored, or skipped
mandatory tests are `critical` and produce a nonzero runner exit status.

Direct optional-dependency skips are prohibited inside a mandatory test. If an
algorithm needs an optional package for its public adapter, its core numerical
oracle should normally exercise a dependency-independent kernel. The
coupled-diagonalization gradient oracle compares the analytic gradient with
central finite differences without constructing a `multidesign` object. The
FPGW oracle independently enumerates every feasible tiny integer flow and
checks the compiled transport solver, four-index objective, mass/capacity
constraints, and permutation equivariance without `lpSolve`.

## Status and severity

| Observation | Severity | Tier result |
|---|---|---|
| Failure or error | `critical` | Fail |
| Missing or skipped mandatory test | `critical` | Fail |
| Warning | `warning` | Pass with visible debt |
| Non-mandatory skip in PR | `info` | Pass with visible coverage gap |
| Non-mandatory skip in nightly/release | `warning` | Pass with visible coverage gap |
| Passed test | `ok` | Pass |

Known gaps are emitted separately. They do not fail experimental methods, but
they block stable promotion. A release run also records `not_available` when a
method has no declared release artifact.

## Accuracy and scale are separate gates

An oracle establishes implementation fidelity on a controlled problem. It does
not establish predictive dominance. Comparative claims require the frozen
benchmark manifests, controls, held-out scoring, multiple seeds, failure-rate
reporting, and resource receipts described in `BENCHMARK_PROTOCOL.md`.

Likewise, a passing toy test is not a scale claim. Nightly and release evidence
sample peak RSS in isolated processes and retain timeout and memory-limit
outcomes as data rather than dropping them.

## Real-data boundary

The release workflow supports the leakage-safe Chikusei blocked-inductive
benchmark, but the data are not bundled. A manual caller must both request the
benchmark and explicitly accept its CC BY 3.0 terms. A tag-triggered release
does not imply consent. When the data are not run, the uploaded
`real-data-status.csv` states that fact.

## Review checklist

Before accepting a new evidence claim, verify that:

- the method's contract record and registry record agree;
- the mandatory test actually ran and did not skip;
- the oracle is independent enough to catch a shared implementation defect;
- no assertion is true by construction, such as a nonnegative correctness
  count;
- warnings and failure rates are included in the artifact;
- benchmark tuning never reads frozen final labels; and
- accuracy, runtime, peak RSS, and unsupported/not-available outcomes are all
  reported rather than selectively omitted.
