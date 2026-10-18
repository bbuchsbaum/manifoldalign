# Alignment benchmark protocol v2

Status: frozen protocol inputs, version 2.0.0, reviewed 2026-07-22.

Protocol v2 replaces the former graph benchmark's paired noisy feature clouds.
Every qualifying graph case now contains observed weighted relations and enough
truth to audit correspondence, overlap, supervision, perturbation severity, and
data leakage. A result produced from `X` alone may remain a useful comparator,
but it is not evidence that a method supports supplied graph topology.

## Case contract

Each `graph_alignment_case_v2` has two domains. Each domain contains:

- a finite numeric feature matrix `X`;
- a sparse edge-list relation table with `from`, `to`, nonnegative `weight`,
  `relation`, and `directed` fields;
- spatial positions and mixed numeric, categorical, logical, and missing node
  attributes;
- observed training labels, immutable row identities, train/dev/final node
  roles, and complementary shared/private-node masks.

The case-level truth contains a source-to-target map, an identity-bearing
correspondence table, complete generator labels, overlap/private masks, and the
anchor subset. Source-private nodes map to `NA`; shared correspondence is a
one-to-one map onto target-shared nodes. Anchors must be true correspondences
drawn from train nodes.

Every case records the generator version, seed, full configuration hash, case
content hash, preprocessing specification, requested perturbations, and
observed perturbation severity. `validate_graph_alignment_case_v2()`
recomputes the hashes and checks all structural invariants.

## Frozen scenarios

The installed manifests are immutable protocol inputs:

| Manifest | Cases | Tuning budget | MD5 |
|---|---:|---:|---|
| `graph_alignment_v2_development_manifest.csv` | 12 | Positive, recorded per case | `6030e2e2d49b0473bba78ed283ea2fe6` |
| `graph_alignment_v2_final_manifest.csv` | 42 | 0 | `2f79e4b34867edba288239fbc427a437` |

The expected hashes are also compiled into the package. Loading a changed
manifest fails until an intentional protocol revision updates both the file and
the expected hash.

The 42 final cases are paired inputs: every method is evaluated on the same
case seeds and truth. They cover each combination needed by the protocol:

| Tier | Purpose | Final cases |
|---|---|---:|
| `oracle` | Tiny cases for exact mapping and solver checks | 6 |
| `small_accuracy` | Repeated accuracy comparisons; 9 partial-overlap and 6 closed-world cases | 15 |
| `medium_robustness` | Severe combined perturbations; 9 partial-overlap and 6 closed-world cases | 15 |
| `large_scaling` | Resource and failure-ceiling evidence | 6 |

Each tier includes community/SBM, geometric/lattice, and hub/scale-free graph
families. The 30 accuracy and robustness cases satisfy the predeclared final-
seed floor. Closed-world cases keep classic methods with no partial-overlap
contract in a nontrivial comparison; partial-overlap cases test honest open-set
behavior without handing those methods generator masks. Final seeds never
appear in the development manifest.

## Perturbation coverage

Manifest configurations independently control:

- edge deletion, addition, and rewiring;
- edge-weight noise and rescaling;
- directed relations, heterophily, and a separately named misleading
  relation;
- coordinate jitter and rigid rotation/translation;
- numeric/categorical attribute noise and missingness;
- unequal domain sizes, partial overlap, class imbalance, domain-private
  classes, and anchor scarcity.

Paired severity cases use common random draws for features, row identities,
anchors, and node splits. Changing edge severity therefore changes topology
without silently changing the feature clouds. Tests enforce the requested and
observed deletion ordering and structural-edge Jaccard response.

## Splits and tuning boundary

Development cases have an explicit positive tuning budget. Final cases have a
budget of zero. Official scoring distinguishes `purpose = "tune"` from
`purpose = "report"`; attempting to tune against a final case raises
`graph_benchmark_final_label_leakage`.

Within a case, preprocessing is fixed as
`train_numeric_center_scale_v1`. Numeric feature centering and scaling are fit
on train nodes only. Mixed attributes retain explicit missing values for
methods that natively support them. Target truth, final labels, and generator
correspondence are forbidden as features or preprocessing inputs.

## Required controls

`benchmark_graph_alignment_controls_v2()` evaluates the following controls on
the same case and evaluation role:

| Control | Task | Uses observed topology |
|---|---|---|
| Raw-attribute nearest neighbour | Node correspondence | No |
| Position-only matching | Node correspondence | No |
| Topology-descriptor-only matching | Node correspondence | Yes |
| Anchor-only prediction | Node correspondence | No |
| Row-identity prediction | Node correspondence | No |
| Source-only classifier | Target classification | No |
| Target-only classifier | Target classification | No |
| Seeded chance prediction | Node correspondence | No |

Control predictions carry deterministic hashes. The source-only classifier is
the no-alignment transfer baseline and is fit only from observed source labels.
The target-only classifier is fit only from observed target training labels;
it is a classification ceiling check, not an alignment method.

## Comparator and result semantics

`benchmark_graph_alignment_methods()` now generates v2 development cases and
records protocol version, case ID/hash, split, tier, family, preprocessing,
adapter, input class, and observed-edge count in every row. Adapter preflight
is driven by the canonical method registry and refuses unsupported signed,
directed, weighted, multirelation, unequal-size, partial-overlap, supervision,
or domain-count contracts before a solver runs.

The input boundary is intentionally strict:

| Comparator | Protocol-v2 input disposition | Supplied relations consumed |
|---|---|---:|
| RAPID-MA | Supplied graph, positions, and mixed attributes | Yes |
| Token-OT, Spectral MNN, SSMA, KEMA, GlobalGeo, LRA, GPCA, CONE, GRASP, PARROT | Feature-derived graph or feature manifold | No |
| GW/FPGW | Metric-measure geometry | No |
| LSE | Single-domain target-similarity estimator; cross-domain result is `unsupported` | No |
| LeMA (Hong et al. reference) | External reference; no reproducible package adapter yet, so result is `not_available` | No claim |

Constructing a k-nearest-neighbour graph from `X` does not turn a method into
an observed-graph method. In particular, CONE, GRASP, and PARROT are never
credited with consuming caller-supplied edges. A deterministic
`adapter_input_hash` covers exactly the data delivered to each method. The
edge-sensitivity conformance test mutates relation weights while holding node
data fixed and requires only RAPID-MA's adapter hash to change.

Every successful prediction is normalized into an identity-bearing table with
source row ID, source role, anchor/private masks, true target row ID, predicted
target row ID, and correctness. The adapter must preserve row identity and the
anchor/non-anchor/private partitions. Top-1 accuracy uses all eligible shared
source nodes as its denominator, so missing predictions count as errors rather
than disappearing. Results separately report coverage, anchor and non-anchor
accuracy, and private-node rejection accuracy.

Method-specific search spaces are versioned and hashed before evaluation.
Development cases receive their recorded search budget; final cases receive
zero tuning trials. All methods share the same per-tier wall-clock ceilings:
30 seconds for oracle, 120 seconds for small accuracy, 600 seconds for medium
robustness, and 1,800 seconds for large scaling. The interactive dispatcher
enforces a best-effort in-process elapsed-time limit and records trials used,
convergence, and typed timeout failures. The official tournament must execute
one method-case pair per external process to obtain comparable peak RSS and
hard timeout enforcement; interactive rows therefore set `peak_rss_mb = NA`
and `memory_measurement = "external_process_required"`.

Typed outcomes are `ok`, `unsupported`, `not_available`, and `error`.
Accuracy is valid only for `ok` rows. Missing anchors, missing optional
dependencies, malformed output, and execution failure retain distinct status
codes and cannot be converted into zero accuracy.

This protocol establishes trustworthy inputs and controls. It does not itself
declare a winning method; that requires the frozen package tournament and
real-data benchmark receipts tracked by the evidence epic.

## Official final tournament

The package-wide final tournament is regenerated from the package root with:

```sh
Rscript inst/benchmarks/run_alignment_tournament.R \
  --profile=final \
  --output=docs/evidence/alignment-tournament-v2 \
  --resume=true
```

The checked receipt contains 588 rows (42 cases by 14 methods), typed outcomes
for every method-case pair, external process runtime and peak RSS, raw controls,
the frozen tuning ledger, conditional and intent-to-run summaries, paired
effects, Pareto tables, scaling fits, session information, and recursive
checksums. Its `validation_problems.csv` has zero rows. The interpretation is
checked beside the raw receipt and must state RAPID-MA's matching,
blocked-inductive, and million-node limitations as well as its supported
regimes; a broad input contract is not reported as a universal accuracy win.
