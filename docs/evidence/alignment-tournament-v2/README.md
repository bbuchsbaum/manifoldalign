# Alignment tournament v2 evidence

Profile: `final`; final manifest MD5: `2f79e4b34867edba288239fbc427a437`.

This receipt contains 588 typed method-case rows across 42 cases and 14 methods.
Failures and timeouts remain in intent-to-run denominators; unsupported tasks remain visible but are not scored as algorithm failures.
Matching uses native top-1 followed by a declared transport or common-embedding ranking. Classification is a frozen source-label downstream probe, not the method's private classifier.
The checked interpretation and lifecycle implications are in `INTERPRETATION.md`.
Checksums cover published receipt files; resumable binary checkpoints and the debug log remain local operational state.

## Status counts

| Status | Count |
|---|---:|
| error | 15 |
| not_applicable | 42 |
| ok | 207 |
| unsupported | 324 |

## Regeneration

```sh
Rscript inst/benchmarks/run_alignment_tournament.R --profile=final --output=/Users/bbuchsbaum/code/manifoldalign/docs/evidence/alignment-tournament-v2 --resume=true
```

See `validation_problems.csv`; an official receipt is valid only when it has zero rows.
