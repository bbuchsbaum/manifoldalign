# Coupled Diagonalization TODO (from paper cross-check)

Context: review against Eynard et al. (Section 3.3, subspace parametrization + optimization).
Scope: `R/coupled_diagonalization.R`.

## Confirmed (no action required)

- Current formulation matches the subspace-parametrized coupled diagonalization objective structure.
- Current diagonalization and coupling gradients are internally consistent with the implemented objective.
- QR re-orthonormalization is a valid Stiefel retraction approach (even without `fmincon`).

## P0: Correctness and robustness

- [ ] Fix Armijo condition mismatch for normalized direction.
  - Location: `R/coupled_diagonalization.R:452`, `R/coupled_diagonalization.R:479`.
  - Issue: update uses `direction <- grad / ||grad||`, but Armijo RHS uses `||grad||^2`.
  - Decide one consistent form:
    - Keep normalized direction and use RHS proportional to `||grad||`, or
    - Use unnormalized direction `-grad` and keep RHS with `||grad||^2`.
  - Acceptance: line search accepts/rejects based on the mathematically matched Armijo model.

- [ ] Scale user `lambda_target` consistently with per-domain eigenvalue scaling.
  - Location: scaling at `R/coupled_diagonalization.R:343-348`; target handling at `R/coupled_diagonalization.R:126-132` and `R/coupled_diagonalization.R:359-367`.
  - Issue: `Lambda[[i]]` is scaled by `scale_factor`, but user-supplied `lambda_target[[i]]` is not.
  - Acceptance: when `alpha_match > 0`, objective units are consistent for both auto and user targets.

- [ ] Resolve semantics of `alpha_match` to avoid accidental off-diagonal double-counting.
  - Location: `R/coupled_diagonalization.R:638-640`, `R/coupled_diagonalization.R:683-685`, `R/coupled_diagonalization.R:709`.
  - Issue: current `match_residual <- Mi` includes off-diagonal terms, so total off-diagonal weight becomes `1 + alpha_match`.
  - Decision:
    - If intentional, document explicitly in roxygen/docs.
    - If not intentional, change match term to only diagonal mismatch.
  - Acceptance: docs and code match intended loss decomposition.

## P1: Convergence quality improvements

- [ ] Switch block update direction to Stiefel tangent (Riemannian) gradient.
  - Location: `R/coupled_diagonalization.R:440-453` and helper `R/coupled_diagonalization.R:609-651`.
  - Change: project Euclidean gradient `G` via `G - A sym(A^T G)` before line search.
  - Acceptance: all updates remain in tangent direction before QR retraction.

- [ ] Add numerical gradient checks for small synthetic cases.
  - New tests under `tests/testthat/`.
  - Validate diagonal, coupling, and `alpha_match` terms (including with/without scaling).
  - Acceptance: finite-difference checks pass within tolerance.

## P2: Major runtime improvements (same objective)

- [ ] Refactor line search to incremental objective updates (single-block delta).
  - Hotspot: full recomputation in backtracking at `R/coupled_diagonalization.R:472-476`.
  - Keep state: per-domain diagonal cost, `B_i`, `sumB`, norms, and current total.
  - Acceptance: line-search trial evaluates only changed block and updated coupling summary.

- [ ] Add exact Gram reduction for coupling to remove `q` from inner loop.
  - Current coupling uses explicit `B_i = F_i A_i` at `R/coupled_diagonalization.R:687-706`.
  - Precompute `G_ij = t(F_i) %*% F_j` once, then use trace/Gram formulas for cost and gradient.
  - Acceptance: coupling cost/gradient no longer require repeated `q x k'` multiplies per trial step.

- [ ] Support cached eigendecompositions / Laplacian basis reuse across sweeps.
  - Add optional inputs for precomputed `Ubar`, `Lambda`, `FiUbar` (or a cache object).
  - Acceptance: tuning `mu_coupling`, `alpha_match`, `ncomp` can skip graph/eigendecomposition rebuild.

## P3: Implementation cleanup

- [ ] Replace `sweep(Ai, 1L, lambda_i, '*')` with rowwise multiply where safe.
  - Locations: `R/coupled_diagonalization.R:621`, `R/coupled_diagonalization.R:670`.

- [ ] Reduce per-iteration allocations by reusing temporary matrices/vectors in optimizer loop.

- [ ] Profile after P2 changes; consider `RcppEigen` for cost/gradient helpers only if R overhead remains dominant.

## Suggested execution order

1. P0 Armijo + `lambda_target` scaling + `alpha_match` semantics.
2. Add/expand tests to lock objective + gradients.
3. P2 incremental line search objective.
4. P2 Gram reduction.
5. P1 Riemannian gradient direction.
6. P2 cache hooks for repeated runs.
7. P3 micro-optimizations / optional C++ move.

## Regression checklist after each stage

- [ ] Orthogonality preserved: each `A_i^T A_i \approx I`.
- [ ] Objective non-increasing across accepted steps.
- [ ] Convergence on existing examples remains stable.
- [ ] Numerical equivalence vs baseline objective (where algorithmic form should be exact).
- [ ] Runtime benchmark recorded for moderate `m`, `n`, `q`, `k'`, `k`.
