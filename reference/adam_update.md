# Adam optimizer update step

Adam optimizer update step

## Usage

``` r
adam_update(
  par,
  grad,
  m_state,
  v_state,
  t,
  lr = 0.001,
  b1 = 0.9,
  b2 = 0.999,
  eps = 1e-08
)
```

## Arguments

- par:

  Current parameter matrix

- grad:

  Gradient matrix

- m_state:

  First moment state (NULL on first call)

- v_state:

  Second moment state (NULL on first call)

- t:

  Iteration number (1-indexed)

- lr:

  Learning rate

- b1:

  Exponential decay rate for first moment

- b2:

  Exponential decay rate for second moment

- eps:

  Small constant for numerical stability

## Value

A list with updated parameters and optimizer state
