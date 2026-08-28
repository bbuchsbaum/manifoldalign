# Alignment tournament v2 evidence

Profile: `final`; final manifest MD5: `2f79e4b34867edba288239fbc427a437`.

This receipt contains 42 typed method-case rows across 42 cases and 1 methods.
Failures and timeouts remain in intent-to-run denominators; unsupported tasks remain visible but are not scored as algorithm failures.
Matching uses native top-1 followed by a declared transport or common-embedding ranking. Classification is a frozen source-label downstream probe, not the method's private classifier.
The checked interpretation and lifecycle implications are in `INTERPRETATION.md`.
Checksums cover published receipt files; resumable binary checkpoints and the debug log remain local operational state.

## Status counts

| Status | Count |
|---|---:|
| ok | 12 |
| unsupported | 30 |

## Regeneration

```sh
Rscript inst/benchmarks/run_alignment_tournament.R --profile=final --output=artifacts/kema-eigencore-benchmark-2026-08-27/tournament-final --resume=true
```

See `validation_problems.csv`; an official receipt is valid only when it has zero rows.
