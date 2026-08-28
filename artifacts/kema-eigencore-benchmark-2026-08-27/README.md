# KEMA/eigencore benchmark check (2026-08-27)

## Conclusion

The repaired KEMA is substantially faster than the historical implementation
for the two RBF variants, but it is not a universal performance win.

- RBF KEMA: large runtime improvements on the historical feature challenges,
  with small average accuracy gains.
- Linear KEMA: slower, with a small average centroid-accuracy regression.
- Frozen tournament: identical predictive scores and lower peak memory, but a
  modest runtime regression on the small supported cases.
- Numerical contract: every supported current fit passed the independent
  residual and B-orthogonality gates.
- Unsupported scope is unchanged: partial/open-set and large-scaling cases are
  not evidence for KEMA.

The historical and current challenge runs exercise different KEMA
implementations, not merely interchangeable eigensolver calls. The challenge
speedups therefore support a claim about the repaired KEMA path as a whole;
they do not isolate eigencore as the sole cause.

## Paired historical feature challenges

The comparison pairs the existing fast-profile results in
`dev-scripts/benchmarks/toy_feature_aligners_results.csv` with the current
KEMA-only run by method, scenario, and seed. There are 45 pairs: three KEMA
variants, five scenarios, and three seeds. Status agrees in all 45 pairs; 36
pairs are successful and the nine errors are all the known partial-open-set
shape failure.

| Variant | Successful pairs | Geometric runtime speedup (old/current) | Mean latent-correlation change | Mean centroid-accuracy change |
|---|---:|---:|---:|---:|
| Linear | 12 | 0.41x | +0.003 | -0.0276 |
| Full RBF | 12 | 100.7x | +0.079 | +0.0143 |
| Reduced RBF | 12 | 65.2x | +0.026 | +0.0486 |

The reduced-RBF arithmetic mean includes a historical 984-second outlier, so
the geometric summaries above are preferred. Its median paired speedup is
96.2x; full RBF's median is 287.8x. Linear KEMA is about 2.42x slower by the
inverse geometric ratio.

The complete per-challenge comparison is in
[`comparison/toy_summary.csv`](comparison/toy_summary.csv), and every paired
seed is in [`comparison/toy_paired.csv`](comparison/toy_paired.csv).

## Frozen alignment tournament

KEMA was rerun alone over the immutable 42-case final manifest (manifest MD5
`2f79e4b34867edba288239fbc427a437`) and paired with the checked historical
tournament artifact.

- Status is identical in all 42 cases: 12 successful and 30 unsupported.
- Top-1 correspondence accuracy is exactly identical in all 12 successful
  pairs.
- Classification overall accuracy is exactly identical in all 12 successful
  pairs.
- Mean internal runtime is 0.209 seconds versus 0.153 seconds historically
  (current is 1.37x slower by the ratio of means).
- Mean external process time is 0.292 seconds versus 0.279 seconds
  historically (current is 4.6% slower by the ratio of means).
- Mean peak RSS is 292.4 MB versus 323.0 MB historically, a 30.6 MB or 9.5%
  reduction.

This is preservation of tournament accuracy and an improvement in memory, not
an improvement in tournament predictive performance. See
[`comparison/tournament_paired.csv`](comparison/tournament_paired.csv) and the
new [`tournament-final/run_receipt.txt`](tournament-final/run_receipt.txt).

## Expanded current feature profile

The full feature profile adds `many_domain_consensus` and `highdim_stress` to
the five fast scenarios. Across seven scenarios, three variants, and three
seeds, 54 of 63 fits succeeded. All 54 successful fits passed the eigensolver
certificate and the independent KEMA fidelity gate. The nine failures were
again confined to `partial_open_set`.

The new larger challenges completed quickly:

| Challenge | Variant | Mean runtime (s) | Mean latent correlation | Mean centroid accuracy |
|---|---|---:|---:|---:|
| Four-domain consensus | Linear | 0.057 | 0.845 | 0.574 |
| Four-domain consensus | Full RBF | 0.115 | 0.100 | 0.840 |
| Four-domain consensus | Reduced RBF | 0.072 | 0.263 | 0.932 |
| High-dimensional stress | Linear | 0.065 | 0.650 | 1.000 |
| High-dimensional stress | Full RBF | 0.068 | 0.703 | 0.993 |
| High-dimensional stress | Reduced RBF | 0.059 | 0.479 | 0.799 |

There is no frozen historical full-profile artifact for those two additional
challenges, so these rows establish current robustness and cost, not an
old-versus-new speedup. Full results and provenance are in
[`toy-feature-full/summary.csv`](toy-feature-full/summary.csv) and
[`toy-feature-full/receipt.txt`](toy-feature-full/receipt.txt).

## Attribution and remaining limits

A partial replay of the clean `5a02dfc` KEMA/PRIMME path reproduced the old
behavior: tiny RBF fits took tens of seconds and returned warning-only fits with
relative residuals around 0.9. It was stopped before repeating the recorded
984-second outlier. Since those old eigenpairs do not meet the current fidelity
contract, a raw eigencore-versus-PRIMME timing would not be a common-answer
comparison unless both solvers operate on the same deflated quotient and pass
the same certificate.

The current `operator_exact` path still materializes the quotient matrices.
The six large-scaling tournament cases therefore remain outside KEMA's current
resource envelope; no large-scale performance claim is made here. The governed
real-data release datasets were also not part of this local comparison.

## Reproduction and verification

Current feature runs use eigencore 1.0.3 and fixed seeds 1, 2, and 3:

```sh
LC_ALL=C LANG=C R_LIBS=/tmp/eigencore-cran-kema-lib \
  Rscript artifacts/kema-eigencore-benchmark-2026-08-27/run_toy_feature_kema.R \
  full artifacts/kema-eigencore-benchmark-2026-08-27/toy-feature-full 1,2,3
```

The frozen tournament was run with:

```sh
LC_ALL=C LANG=C R_LIBS=/tmp/eigencore-cran-kema-lib \
  Rscript inst/benchmarks/run_alignment_tournament.R \
  --profile=final --methods=kema \
  --output=artifacts/kema-eigencore-benchmark-2026-08-27/tournament-final \
  --resume=false
```

Validate the evidence-level invariants with:

```sh
LC_ALL=C LANG=C \
  Rscript artifacts/kema-eigencore-benchmark-2026-08-27/verify_results.R
```
