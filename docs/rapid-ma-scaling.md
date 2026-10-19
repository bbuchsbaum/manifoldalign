# RAPID-MA complete-pipeline scaling

> **Regenerated 2026-10-19** with the documented runner and frozen budgets
> (`RAPID_MA_PROFILE=full`, ten-second common-task cap, four-minute cap for the
> RAPID-only sizes) at commit `a6e1cd63b031b5c1f074813a48c5ca07d6479cb8`, on a
> 4-core x86_64 Linux machine (R 4.5.3, OpenBLAS, 15.7 GB) with no other jobs
> running. All 84 tasks are recorded with status `ok`, `timeout`, or `failed`.
> This commit includes exact RAPID-MA speedups (a blocked exact kNN search
> for views with five or more columns and leaner sparse algebra; default fits
> are bit-identical on the benchmark fixtures), the relation-safety gate fix,
> KEMA's preconditioned and matrix-free solvers, and deterministic `irlba`
> starts. The superseded run used a 14-logical-core arm64 macOS machine with
> 36 GB; its summary is kept as `rapid-ma-scaling-summary-previous-macos.csv`.
> Because the caps are wall-clock limits, which tasks finish depends on the
> hardware: the tables below are a within-machine comparison.
>
> An earlier regeneration on this machine at commit `91f23f0`, before the
> speedups, took 180.5 seconds at 50,000 common nodes and exceeded the
> four-minute cap at 100,000; a same-machine control with the pre-fix commit
> `78fbc05` was within 5% of that run
> (`rapid-ma-scaling-control-prefix-summary.csv`).

Status: 220,000-node synthetic run complete on this machine in 149 seconds;
the 1,000,002-node gate failed on the previous machine and was not repeated.

## Bottom line

RAPID-MA's complete validated-teacher pipeline scales approximately linearly
over the measured range. On this machine it fit 100,000 source plus 120,000
target nodes (220,000 total) in 149.4 seconds, retained a 462.8 MB result, and
reached a 1.75 GB sampled process peak. The log-log time slope from 2,500
through 100,000 common nodes is 0.98; the retained-object slope is 0.99. At
50,000 common nodes it took 69.0 seconds, 2.6 times faster than before the
speedups (180.5 seconds) with the same accuracy.

The 10,000-common-node task still exceeded the ten-second common cap
(10.01 seconds).

The speed claims need qualification:

- RAPID-MA is slower than SSMA at every size where both finish (7.45 versus
  5.61 seconds at 5,000 common nodes). It is slower than Spectral MNN up to
  2,500 common nodes and faster at 5,000 (7.45 versus 8.37). None of the
  three finishes 10,000 common nodes in ten seconds.
- It is slightly faster than Token-OT Graph from 1,000 common nodes (1.71
  versus 1.82 seconds; 3.87 versus 4.50 at 2,500), and Token-OT exceeds ten
  seconds from 5,000.
- RAPID-MA is faster than dense LeMA at every size (1.15 versus 1.33
  seconds at 250 common nodes, 1.48 versus 4.82 at 500), and LeMA exceeds ten
  seconds from 1,000 onward.
- KEMA completes 250 and 500 common nodes (1.24 and 5.27 seconds, OA 0.992 and
  0.988), exceeds ten seconds from 1,000 to 5,000, and at 10,000 is refused
  before allocation because even its matrix-free `operator_exact` backend
  needs dense per-domain kernel blocks (about 22.8 GB).
- Multiscale completes through 2,500 common nodes and exceeds ten seconds
  from 5,000 (it finished 5,000 in 7.98 seconds in the earlier regeneration;
  its eigensolves now check convergence and use deterministic starts).
- No direct 10x runtime comparison was completed at 100,000 common nodes.

## Protocol

- Same deterministic synthetic fixture and fixed 12-label-per-domain budget at
  every size.
- Common-node sizes: 250, 500, 1,000, 2,500, 5,000, and 10,000. The target has
  20% additional unmatched nodes, so total nodes are 2.2 times the common size.
- Three repetitions through 1,000; one repetition at larger sizes.
- RAPID-only extensions at 50,000 and 100,000 common nodes.
- Every task creates its fixture and runs in an isolated process. RSS is sampled
  every 20 ms. Common tasks have a ten-second cap; RAPID large tasks have a
  four-minute cap. Timeout and failure rows remain in the data.
- R 4.5.3 on x86_64 Linux, 4 logical cores, 15.7 GB system memory,
  `OMP_NUM_THREADS=1`.
- Hits@1 is exact. Matching MRR is the declared bounded MRR@100 lower bound;
  the evaluator does not materialize an all-pairs distance matrix.

## RAPID-MA scaling

| Common nodes | Total nodes | Runtime (s) | Retained result (MB) | Sampled peak RSS (MB) | OA | Macro-F1 |
|---:|---:|---:|---:|---:|---:|---:|
| 250 | 550 | 1.145 | 1.57 | 306.8 | 0.9826 | 0.9757 |
| 500 | 1,100 | 1.482 | 2.66 | 273.3 | 0.9654 | 0.9615 |
| 1,000 | 2,200 | 1.713 | 5.05 | 312.9 | 0.9888 | 0.9850 |
| 2,500 | 5,500 | 3.866 | 12.01 | 373.0 | 0.9719 | 0.9682 |
| 5,000 | 11,000 | 7.445 | 23.54 | 396.5 | 0.9531 | 0.9577 |
| 10,000 | 22,000 | >10 (timeout) | NA | 490.4 | NA | NA |
| 50,000 | 110,000 | 69.010 | 231.65 | 1,095.7 | 0.9892 | 0.9867 |
| 100,000 | 220,000 | 149.378 | 462.79 | 1,751.9 | 0.9670 | 0.9669 |

Absolute peak RSS includes inherited R and loaded-package state. Retained state
is the cleaner scaling signal and grows essentially linearly.

## Runtime comparison

| Method | 250 | 500 | 1,000 | 2,500 | 5,000 | 10,000 |
|---|---:|---:|---:|---:|---:|---:|
| RAPID-MA teacher | 1.145 | 1.482 | 1.713 | 3.866 | 7.445 | timeout |
| Dense LeMA | 1.334 | 4.815 | timeout | timeout | timeout | timeout |
| KEMA | 1.235 | 5.270 | timeout | timeout | timeout | refused |
| SSMA | **0.312** | **0.422** | **0.532** | **1.548** | **5.610** | timeout |
| Spectral MNN | 0.363 | 0.402 | 0.724 | 2.084 | 8.371 | timeout |
| Token-OT Graph | 1.209 | 1.258 | 1.819 | 4.504 | timeout | timeout |
| Multiscale | 0.497 | 0.571 | 0.786 | 2.770 | timeout | timeout |

Times through 1,000 are three-run means; larger entries are single runs.
"Timeout" means the isolated end-to-end task exceeded ten seconds on this
machine, not that the method can never complete with a larger budget.
"Refused" is KEMA's pre-allocation memory guard (recorded as `failed` with
its message in the raw file).

## Dense-memory contrast

The LeMA reference's learned graph is dense over target nodes plus paired
source landmarks. With 100,000 common nodes and 20% target-only expansion, that
graph has about 120,012 rows. One double-precision matrix of that shape requires
about 115.2 GB; retaining both `W` and `L` alone requires about 230.4 GB before
distance workspaces and embeddings, far beyond either machine.

At 10,000 common nodes, LeMA had reached 3.65 GB sampled RSS when the
ten-second timeout stopped it on this machine. RAPID's sparse retained result
is 462.8 MB at 100,000 common nodes, with a 1.75 GB sampled peak.

## Scaling contracts and limitations

Automated tests double fixture size and check that retained state grows by less
than 3x, transport nonzeros remain bounded by `q * n + K`, relation nonzeros
remain bounded by `degree_cap * n`, and label budget does not grow with size.

What this establishes:

- complete-pipeline operation at 220,000 total nodes on a 4-core machine;
- near-linear runtime (slope 0.98) and retained state (slope 0.99);
- RAPID-MA is faster than dense LeMA at every measured size;
- a crossover over Spectral MNN between 2,500 and 5,000 common nodes, but
  none over SSMA within the ten-second budget on this machine;
- no hidden quadratic matching evaluator in the benchmark.

What it does not establish:

- a completed million-node full fit (the frozen attempt timed out);
- distributed, disk-backed, or GPU execution;
- a direct 10x runtime win over LeMA at 100,000 common nodes;
- real-scene accuracy at large scale;
- stable accuracy from a fixed 12-label budget as sample size grows.

The following diagnosis comes from the previous machine and was not
repeated here. The measured million-node attempt changes the diagnosis. Memory did not
approach the machine ceiling: peak RSS was 7.20 GB. A same-control sampling
profile at 100,000 common nodes attributes 57.74% of sampled time to bounded
relation construction, including 54.41% in `RANN::nn2`, and 26.39% to sparse
UOT evaluation. These two paths account for 84.13% of sampled time.

The million task has 4.545x as many common nodes as the 100,000-common-node
profile but consumed more than 8.45x its profiled fit time. Because the larger
task is right-censored, the corresponding exponent of 1.410 is a lower bound.
The next implementation work is reusable approximate relation indexes (plus a
linear spatial-grid path) and warm-started UOT with cached candidate costs.
Both need frozen small-scale accuracy/oracle checks before the exact same
million-node gate is repeated.

## Artifacts

- `rapid-ma-scaling-raw.csv`: every task, timeout/failure, metric, and sampled
  peak.
- `rapid-ma-scaling-summary.csv`: per-method/per-size aggregation.
- `rapid-ma-scaling-metadata.txt`: hardware and budgets (the runner was
  executed outside the package root, so commit and source hashes are in
  `rapid-ma-regeneration-provenance.txt`).
- `rapid-ma-scaling-control-prefix-summary.csv`: same-machine control with
  the pre-fix commit.
- `rapid-ma-scaling-summary-previous-macos.csv`: the superseded summary from
  the previous machine.
- `rapid-ma-million-{raw,summary}.csv` and
  `rapid-ma-million-metadata.txt`: the frozen 1,000,002-node timeout
  (previous machine, not repeated).
- `rapid-ma-million-profile.txt`: the compact 100,000-common-node function
  profile and filing-ready bottleneck targets.

SHA-256:

```text
720cfec700d90c036f3ef4330996af54714d8c3fe3459a29ae9fa1a1ccd3e989  rapid-ma-scaling-raw.csv
66fc7cbadfa7795421ed4366d44eec442ab4f3ec5e63e76476914163d666f32f  rapid-ma-scaling-summary.csv
9f24c8bdd3efbc644f01d1d4fc6b5f81d59683a66cd582d3778ffcf4b6da875e  rapid-ma-scaling-metadata.txt
985e87c03c1b7412a24cde3e17fb6b2e1db5f1812cbd5f39358d4b66d3393428  rapid-ma-scaling-control-prefix-summary.csv
e71be1565d88563813976e28b39a98eb768647ebf793c691a3a82840140839d9  rapid-ma-scaling-summary-previous-macos.csv
5daae21f9402b7c594278cb1fb1bd14a513950816718075436f3426c917c5a96  rapid-ma-regeneration-provenance.txt
db14e271af34a2ce1d110d2e729ca7466cac29a2a7b98c5d96e4de53e7c362de  rapid-ma-million-raw.csv
4753eca1622f6033e15204558f756c3a68d9e0511a8b6b7b2770b7406df86328  rapid-ma-million-summary.csv
3c6c6114d0f2732bd67a2700db69253f51c63a2ee19e260e5ae714c6d5370cc7  rapid-ma-million-metadata.txt
9c2074ba81eba487bfb51c16813ae7baa1e218c8d4861dcfe5012ddf698951d7  rapid-ma-million-profile.txt
```
