# Governed real-data alignment benchmarks

The real-data suite complements the frozen synthetic graph-v2 protocol. Its
job is to answer two different questions without mixing them:

1. Do alignment methods preserve known node correspondence when observed
   attributes, weighted edges, direction, partial overlap, and unequal graph
   sizes come from real networks?
2. Does cross-sensor transfer survive a geographically blocked, inductive
   remote-sensing split?

Protocol `1.0.0` freezes task definitions in
`inst/benchmarks/real_alignment_manifest.csv` and external artifacts in
`inst/benchmarks/real_alignment_artifacts.csv`. Results without those task,
revision, checksum, split, and evidence-scope fields are not comparable to the
governed suite.

## Frozen tasks

| Task | Source -> target | Signal under stress | Evaluation boundary |
|---|---|---|---|
| `so2sat_lcz42_sar_optical_culture10` | Sentinel-1 SAR validation-west -> Sentinel-2 optical testing-east | Sensor statistics can conflict; geography changes | Fit only on validation-west paired sensors and labels; open testing-east labels only to the scorer |
| `aucs_coauthor_lunch` | AUCS coauthor -> lunch layers | Actor attributes are informative | Hash-frozen training anchors; score known identities including natural target-private actors |
| `monastery_like_dislike` | Sampson `like1` -> `dislike` | Directed, rank-weighted topology is the main signal | Hash-frozen training anchors; score known identities with natural unequal active-node sets |

The AUCS and Sampson tasks use observed layers from `multinet` 4.3.4. They are
not synthetic rewires. Actor names define ground truth, but are retained only
in audit and scoring fields; model matrices do not contain actor IDs. Nodes
inactive in a selected layer are not padded back into that graph. This gives
the tasks their natural unequal sizes and partial overlap.

So2Sat-LCZ42 v4.2 contains paired Sentinel-1 and Sentinel-2 patches across 17
local-climate-zone classes. The Culture-10 geography split is retained:
validation cities are the development region and testing cities are the frozen
east-region final evaluation. The feature extractor summarizes the actual
sensor patches; it does not synthesize observations or infer labels.
Means, standard deviations, and extrema use all finite values in each patch;
quantiles use a fixed evenly spaced 512-value grid to keep full extraction
linear-time and deterministic.

## Provenance and redistribution

| Dataset | Frozen source | License | Package policy |
|---|---|---|---|
| So2Sat-LCZ42 v4.2 | Hugging Face revision `b5c817486899935e864832b93086cc87f3eee473` | CC BY 4.0, plus applicable Copernicus terms | Raw HDF5 artifacts stay in a user cache; only a small attributed feature fixture is redistributed |
| AUCS and Sampson monastery | CRAN `multinet` 4.3.4 source archive | Apache License 2.0 | Selected `.mpx` files are redistributed with exact source checksums |

Every downloaded artifact has an expected byte count and SHA-256 digest. The
acquisition script uses resumable downloads and rejects a file before it enters
the verified cache if either check differs.

## Acquisition

Review the dataset terms, then acquire only the So2Sat validation and testing
artifacts needed by the governed task:

```sh
export MANIFOLDALIGN_ACCEPT_SO2SAT_LICENSE=yes
Rscript inst/benchmarks/prepare_real_alignment_data.R \
  --output-dir=/path/to/alignment-real-cache \
  --download-so2sat \
  --so2sat-artifacts=validation,validation_geo,testing,testing_geo
```

Network fixtures can be regenerated from the same versioned script:

```sh
Rscript inst/benchmarks/prepare_real_alignment_data.R \
  --output-dir=/path/to/alignment-real-cache \
  --extract-network-fixtures-to=inst/extdata/alignment-real
```

## CI fixture versus full accuracy evidence

The packaged `so2sat_lcz42_ci_features.csv` contains real patches selected
deterministically from every observed class in both frozen splits. It is
class-balanced and intentionally small. It proves schema, identity,
determinism, and leakage contracts in pull-request CI. Its accuracy values are
diagnostics only: they must never be reported as dataset-level accuracy or used
to rank methods.

Release-grade remote-sensing accuracy uses every validation and testing row:

```sh
Rscript inst/benchmarks/run_real_alignment_benchmarks.R \
  --output-dir=artifacts/alignment-real-release \
  --data-dir=/path/to/alignment-real-cache \
  --full-so2sat \
  --strict-data
```

Use `--rebuild-so2sat-features` to regenerate the full feature artifact even if
one is already present. The runner verifies the four raw input artifacts,
derives features with `--per-class=all`, records the derived SHA-256 digest,
and emits separate `ci_schema_only` and `full_frozen_split` control rows.
Extractor v2 yields 24,119 validation rows and 24,188 testing rows. Its frozen
feature receipt is 19,531,044 bytes with SHA-256
`9eeb19fce6bb45b1fd4d3c685d03b6f4b8704dd6ec1cd6de059c7aa699fdce78`.
The full feature CSV is a generated release artifact and is not packaged.

## Leakage boundary

The So2Sat case object deliberately separates `fit` and `evaluation` views.
Development data include paired validation-west SAR/optical features and
labels. The fit view receives testing-east optical features, positions, cities,
and row IDs for prediction, but never receives final target labels. The scorer
alone owns testing-east labels and verifies identical final row IDs before
computing accuracy, balanced accuracy, and macro F1.

For networks, training anchors are a deterministic actor-hash subset. Remaining
known identities are dev/final truth. Actor IDs are prohibited in model feature
columns and final labels are unavailable. Supplied relation tables preserve
edge endpoints, direction, and weights from the governed `.mpx` files.

## Required controls

The runner emits controls before package-method comparisons:

- raw-attribute nearest neighbour or source-only classification;
- topology-descriptor-only correspondence when an observed graph exists;
- position-only classification when spatial coordinates and a same-city
  development reference exist;
- row-index correspondence, which exposes accidental ordering leakage;
- target-only classification as an upper-context control for So2Sat;
- deterministic chance; and
- typed `unsupported` rows when a control is incoherent for a task.

Controls are not promoted as alignment methods. They establish whether a task
contains usable attribute, topology, or spatial signal and whether an alleged
gain is better than ordering or chance.

## Evidence and failure semantics

Each run writes registries, case metadata, raw controls, dataset status,
session information, and SHA-256 receipts. A missing or corrupt source artifact
is `dataset_unavailable` or `dataset_invalid`; malformed method output or a
control exception is `algorithm_error`. The runner never turns a dataset
failure into zero accuracy, and it never treats an unrequested full So2Sat run
as completed evidence. `--strict-data` makes dataset failures nonzero for
release automation.

Key outputs are:

- `real-dataset-status.csv`
- `real-case-metadata.csv`
- `real-network-controls.csv`
- `real-so2sat-controls.csv`
- `real-so2sat-input-validation.csv` for full runs
- `run-receipt.csv`

## Limitations

The packaged networks are small correctness and task-realism benchmarks, not
large-scale timing workloads. So2Sat patch summaries make a reproducible
cross-sensor benchmark practical but do not replace end-to-end image models.
The CI subset is deliberately unsuitable for comparative accuracy claims.
Large-scale runtime remains governed by the synthetic graph-v2 scaling tier,
while method ranking must report both frozen synthetic and full real-data
results.

## Sources

- So2Sat-LCZ42 repository and split documentation:
  <https://github.com/zhu-xlab/So2Sat-LCZ42>
- Frozen So2Sat dataset artifacts:
  <https://huggingface.co/datasets/zhu-xlab/So2Sat-LCZ42>
- `multinet` package and data documentation:
  <https://cran.r-project.org/package=multinet>
