# Repository instructions

These instructions apply to the entire repository. `AGENTS.md` is the shared
source of truth for coding agents; `CLAUDE.md` imports this file. Put durable,
repository-wide guidance here, not transient task status or personal settings.

## Package contract

- `manifoldalign` is an R package for manifold alignment and domain adaptation.
- Unless an API explicitly says otherwise, matrices use rows for samples and
  columns for features. Preserve row identity and domain membership across
  fitting, transformation, validation, and benchmarking.
- Sparse `Matrix` objects are normal package values. Do not coerce them to base
  matrices merely to satisfy a test; accept both representations when the API
  contract permits both, and preserve sparsity in performance-sensitive code.
- `adjoin` is the CRAN graph-construction dependency formerly named
  `neighborweights`. New code, tests, documentation, and imports must use
  `adjoin`. Treat remaining `neighborweights` references as migration debt. If
  a touched path still uses the old name, migrate its coherent code, test,
  roxygen, and namespace slice together.
- Alignment methods have explicit lifecycle states. Do not call a method stable
  or recommended unless `alignment_method_registry()` and the gates in
  `docs/ALIGNMENT_LIFECYCLE.md` support that claim.

## Start in the right place

- `R/align_framework.R` defines the shared aligner interface.
- `R/alignment_registry.R` and `R/alignment_lifecycle.R` own method capability
  and lifecycle policy.
- `R/kema.R` and `R/kema_orig.R` own the fidelity-gated KEMA interfaces and the
  original generalized eigenproblem.
- `R/fpgw*.R`, `R/rapid_ma*.R`, and the other method-named files own their
  numerical implementations and adapters.
- `tests/testthat/ALIGNMENT_TEST_MATRIX.yml` is the machine-readable CI
  contract. `tests/TESTING_STRATEGY.md` explains the evidence tiers.
- `docs/ORACLE_CI_POLICY.md` and `docs/BENCHMARK_PROTOCOL.md` define the
  correctness and comparative-evidence gates.

Search for the owning implementation and its tests before changing an API.
Prefer extending an existing abstraction over adding a parallel registry,
result format, or policy layer.

## Work safely

- Start with `git status --short --branch` and inspect the relevant diff. The
  worktree may contain unrelated user changes; do not overwrite, revert,
  reformat, or stage them.
- Keep edits scoped to the requested behavior. Do not commit, push, change
  dependencies, or regenerate expensive evidence unless the user asks or the
  task clearly requires it.
- Read the exact failure first and reproduce it through a public or package
  entry point. Check class, shape, names, finiteness, and sparse/dense type
  before proposing a numerical explanation.
- Fix causes at the owning boundary. Do not weaken assertions, tolerances,
  lifecycle gates, warnings, or skips to make a check green.
- Preserve backwards compatibility only when it is part of the documented
  contract. Reject withdrawn or ignored arguments explicitly rather than
  silently accepting them.

## R package workflow

Use roxygen comments in `R/` as the source for public documentation and
namespace declarations. Do not hand-edit `NAMESPACE` or generated `man/*.Rd`
files. After changing exports, signatures, parameters, or roxygen imports, run:

```sh
Rscript -e 'devtools::document()'
```

Run the narrowest relevant test first, then broaden in proportion to the
change:

```sh
Rscript -e 'testthat::test_local(filter = "kema")'
Rscript -e 'devtools::test()'
Rscript tools/run-alignment-test-tier.R pr
```

The PR evidence tier is required for changes to alignment behavior, method
contracts, registries, or computational tests. Nightly and release tiers are
intentionally expensive; run them only when their evidence is in scope.

For package-level validation, check the built source package rather than only
the loaded checkout:

```sh
Rscript -e 'rcmdcheck::rcmdcheck(args = "--as-cran", error_on = "warning")'
```

Report focused tests, the full suite, the PR evidence tier, and `R CMD check`
separately. If a check was not run, say so.

## Computational and scientific changes

- State the estimand, input contract, comparator, and failure semantics before
  changing a numerical method or benchmark.
- Test numerical code against a tiny exact or independently implemented oracle
  where feasible. Add relevant invariants, metamorphic cases, adversarial
  inputs, and a regression test for the observed defect.
- Use fixed seeds for stochastic tests and record seeds, parameters, package
  versions, and checksums in published evidence.
- Keep unsupported, unavailable, skipped, timed-out, and errored outcomes
  distinct. Never drop failed cases from a denominator or reinterpret missing
  predictions as successes.
- Prevent label, correspondence, split, and row-identity leakage. Fit
  preprocessing only on the permitted training partition.
- Treat accuracy, numerical correctness, runtime, and memory as separate
  claims. A favorable metric does not override failed validation, warnings, or
  an invalid comparator.
- Preserve protocol-frozen manifests and receipts. Regenerate evidence only
  with its documented runner, then validate row completeness, problem tables,
  provenance, and checksums before changing an interpretation.

## Documentation and generated output

- Update reader-facing sources when public behavior changes: roxygen comments,
  `README.md`, vignettes, and the hand-authored policy files under `docs/`.
- Build the site with `Rscript -e 'pkgdown::build_site()'` when website output is
  in scope. Do not hand-edit generated HTML, search indexes, or site assets.
- Do not present a local test, benchmark, site build, or package check as hosted
  CI, CRAN, or published-site evidence.

## Completion checklist

Before handing off a change:

1. Review `git diff --check`, `git status --short`, and the final scoped diff.
2. Run the focused checks required by the changed behavior.
3. Regenerate and review documentation only when its source changed.
4. Confirm that no unrelated files or generated artifacts were altered.
5. Summarize behavior changed, validation run, and any remaining uncertainty.
