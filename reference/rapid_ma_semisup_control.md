# Control RAPID-MA semi-supervised teacher refinement

Control RAPID-MA semi-supervised teacher refinement

## Usage

``` r
rapid_ma_semisup_control(
  rounds = 3L,
  ema_decay = 0.8,
  confidence_threshold = 0.8,
  minimum_threshold = 0.6,
  balance_strength = 0.08,
  max_pseudo_fraction = 0.5,
  min_margin = 0.15,
  uncertainty_power = 1,
  min_weight = 0.05,
  relation_consistency = TRUE,
  relation_validation = FALSE,
  relation_validation_margin = 1/6,
  relation_validation_max_relations = 6L,
  propagation_steps = 2L,
  propagation_strength = 0.35,
  min_relation_agreement = 0.5,
  relation_gate_min = 0.05,
  label_weight_multiplier = 1.5,
  consistency_weight = 0.5,
  relation_weight = 0.1,
  rollback_tolerance = 1e-06
)
```

## Arguments

- rounds:

  Maximum pseudo-label refresh rounds.

- ema_decay:

  EMA teacher retention in \`\[0, 1)\`.

- confidence_threshold:

  Base class confidence threshold.

- minimum_threshold:

  Lower bound for minority-class thresholds.

- balance_strength:

  Amount by which rare classes receive lower thresholds.

- max_pseudo_fraction:

  Maximum fraction of eligible nodes selected per round.

- min_margin:

  Minimum top-one versus top-two probability margin.

- uncertainty_power:

  Exponent applied to confidence-above-threshold weights.

- min_weight:

  Minimum accepted pseudo-label weight.

- relation_consistency:

  Require agreement with trusted relation neighbours.

- relation_validation:

  Use explicitly supplied out-of-fold probabilities to guard against
  misleading relation channels.

- relation_validation_margin:

  Minimum balanced-accuracy improvement over the full relation set
  before selecting a safer subset.

- relation_validation_max_relations:

  Maximum highest-weight relations considered when constructing bounded
  single and pair candidates.

- propagation_steps:

  Fixed number of sparse relation-propagation steps.

- propagation_strength:

  Convex weight assigned to relation neighbours.

- min_relation_agreement:

  Minimum neighbour probability for the candidate class.

- relation_gate_min:

  Ignore relations below this learned gate weight.

- label_weight_multiplier:

  Multiplier for pseudo-label transport costs.

- consistency_weight, relation_weight:

  Nonnegative refinement score weights.

- rollback_tolerance:

  Allowed score increase before rollback.

## Value

A \`rapid_ma_semisup_control\` object.
