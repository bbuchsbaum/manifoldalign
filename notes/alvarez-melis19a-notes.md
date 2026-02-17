# Notes: Alvarez-Melis, Jegelka, Jaakkola (AISTATS 2019)

**Paper:** *Towards Optimal Transport with Global Invariances* (AISTATS 2019), David Alvarez-Melis, Stefanie Jegelka, Tommi S. Jaakkola.

These are section-by-section reading notes on `alvarez-melis19a.pdf` (10 pages).

---

## 0. One-paragraph summary

Classic discrete optimal transport (OT) can recover correspondences between two point sets *if* cross-set distances are meaningful. With learned representations (e.g., word embeddings), the two spaces are often only identifiable up to a *latent global transformation* (rotation/reflection/etc.), so naïve cross-space costs like `||x - y||` can be meaningless and OT can fail. The paper proposes **OT with global invariances**: jointly optimize a transport coupling `π` and a global map `P` (from a user-chosen invariance class) so that correspondences are computed after accounting for the latent transformation. The resulting objective admits efficient alternating updates: `π` via (entropic) OT/Sinkhorn; `P` via a small SVD-based closed form over Schatten-norm constraint sets. A key special case (`p = ∞`) recovers orthogonal (Procrustes) invariance; another (`p = 2`) connects directly to the Gromov–Wasserstein objective (for cosine similarity + squared loss). They show synthetic point-cloud and unsupervised bilingual lexicon induction results competitive with prior work, often at lower compute and with less fragile initialization (via entropy annealing).

---

## 1. Motivation + problem statement (Intro)

- OT provides:
  1) a way to compute **soft correspondences** between sets (the coupling), and
  2) a **distributional distance** between empirical measures (Wasserstein cost) useful as a training signal (e.g., WGANs).
- Limitation highlighted: OT is “locally greedy” (optimizes individual mass movements) and implicitly assumes **cross-space distances are well-defined**.
- In learned representation settings (e.g., word embeddings), coordinate systems can differ by a **global transform**:
  - embeddings are typically only identifiable up to rotation/reflection (and sometimes more),
  - thus pairwise distances *across* embedding runs/languages can be ill-defined.
- Goal: incorporate a *latent* global transform directly into OT, so couplings become **invariant to such transforms**.

Main contributions listed:
- A formulation of discrete OT that incorporates global geometric invariances.
- Efficient algorithms for the resulting class of problems.
- Application to unsupervised word translation with comparable accuracy at lower runtime and without ad-hoc initialization.

---

## 2. Related work (high level)

- Prior “OT + alignment” shape/feature matching: alternating between a soft assignment/coupling and a global transform estimate (e.g., Softassign Procrustes).
- Word translation alignment:
  - OT/Wasserstein + Procrustes-style alignment approaches exist, but often rely on heuristic initialization (e.g., adversarial training) and can be brittle.
- Gromov–Wasserstein (GW):
  - avoids cross-space distances altogether by matching **intra-space relational structure**,
  - useful when spaces are unregistered, but does not produce a direct cross-space ground metric.
- This paper frames its approach as a compromise:
  - still uses a cross-space cost, but makes it invariant to a structured latent transform,
  - yields a global map `P` as part of optimization (useful for out-of-sample mapping).

---

## 3. Background: alignment + OT (Section 3)

### 3.1 Notation

- Point sets: `X = {x^(i)}_{i=1..n}` and `Y = {y^(j)}_{j=1..m}`.
- Empirical measures:
  - `μ = Σ_i p_i δ_{x^(i)}`, `ν = Σ_j q_j δ_{y^(j)}`.
- Transport polytope (couplings with given marginals):
  - `Π(p,q) = { π ∈ R_+^{n×m} : π 1 = p, π^T 1 = q }`.
- Frobenius inner product: `⟨A,B⟩ = Σ_{i,j} A_{ij} B_{ij}`.

### 3.2 Alignment from paired samples: Procrustes

- If points are paired (supervised) and dimensions match, one can learn a global map `T` by minimizing a matrix discrepancy, often Frobenius norm.
- Orthogonal Procrustes:
  - `min_{P ∈ O_d} ||X - P Y||_F^2`
  - closed form via SVD of `X Y^T` (`P* = U V^T`).
- Emphasis: Procrustes requires paired columns → inherently supervised.

### 3.3 Correspondences between aligned spaces: OT

- Discrete OT (Kantorovich relaxation) solves:
  - `min_{π ∈ Π(p,q)} ⟨π, C⟩`, where `C_{ij} = c(x^(i), y^(j))`.
- Entropic regularization adds `H(π)` to enable fast Sinkhorn scaling, producing dense couplings and efficient computation.
- Core issue: when `X` and `Y` are not registered, a naïve cross-space cost `c(x,y)=||x-y||` can be meaningless.

---

## 4. OT with global invariances (Section 4)

### 4.1 Joint objective: coupling + transform

They propose to **jointly** solve for:
- a coupling `π ∈ Π(p,q)` (local correspondences), and
- a global transform `f ∈ F` (space registration / invariance),

by minimizing transport cost after transforming one set:

`min_{π ∈ Π(p,q)} min_{f ∈ F} ⟨π, C(X, f(Y))⟩`  (Eq. 7)

with costs `C(X,f(Y))_{ij} = c(x^(i), f(y^(j)))`.

They focus on **linear** transforms `f(y) = P y` and define invariance classes by Schatten norms:

`F_p := { P ∈ R^{d×d} : ||P||_{S_p} ≤ k_p }`  (Eq. 8)

where `||P||_{S_p} = ||σ(P)||_p` is the `ℓ_p` norm of singular values.

Interpretation of `p` (via singular-value geometry; cf. Fig. 1):
- `p = 1` (nuclear/trace norm ball): encourages *sparse singular values* → low-rank-ish transforms (e.g., projections / subspace relations).
- `p = ∞` (spectral/operator norm ball): extreme points correspond to (partial) isometries; with radius 1 it naturally encodes rigid transforms (orthogonal matrices).
- `p = 2` (Frobenius): “energy” constraint; connects to GW under certain similarity choices.

### 4.2 Squared Euclidean ground cost and algebraic simplification

They specialize to:
- `c(x, Py) = ||x - Py||_2^2`.

Expanding the square yields an objective that can be written (after collecting constants) as a maximization involving:
- a correlation term `⟨π, X^T P Y⟩`,
- normalization terms involving `||x||^2` and `||Py||^2`.

Under either of two conditions (Lemma 4.1):
1) all `P ∈ F` are angle-preserving, *or*
2) `Y` is (weighted) whitened and all `P ∈ F` have equal Frobenius norm,

the dependence of the `||Py||^2` term can be removed, and the joint problem becomes:

`max_{π ∈ Π(p,q)} max_{P ∈ F} ⟨ X π Y^T, P ⟩`  (Eq. 10)

This is the key “clean” form enabling simple alternating optimization.

### 4.3 Closed-form `P` update via SVD (Lemma 4.2)

Let `M := X π Y^T` with SVD `M = U diag(σ) V^T`.

Maximizing a linear objective `⟨P, M⟩` over a Schatten-`p` norm ball has a solution that:
- uses the same singular vectors (`U, V`), and
- chooses singular values on the boundary according to dual norm geometry.

Concretely (schematically):
- `P* = U diag(s) V^T`,
- `s` is chosen to satisfy `||s||_p ≤ k_p` and maximize `s^T σ`.

This generalizes the classic Procrustes result (`p = ∞` case).

### 4.4 Salient special cases

#### Case A: `p = ∞` (spectral norm ball → orthogonal invariance)

- Schatten-∞ is the **spectral norm** `||P||_2 = max σ(P)`.
- Choosing `k_∞ = 1` ensures `I ∈ F_∞`.
- Combined with the geometric conditions in Lemma 4.1, this recovers the orthogonal group (`F_∞ = O_d`).

Resulting objective after eliminating `P`:
- `max_{π ∈ Π(p,q)} || X π Y^T ||_*` (nuclear norm) (Eq. 13)

Geometric reading they highlight:
- If `p,q` are uniform and `X̂ = X π` denotes the barycentric-transported source, then `X̂ Y^T` is essentially a cross-covariance matrix.
- Optimizing OT in this form encourages couplings that maximize cross-feature correlations after transportation.

Practical update for `P` given `π`:
- If `M = U Σ V^T`, then `P* = U V^T` (Procrustes).

#### Case B: `p = 2` (Frobenius → connection to Gromov–Wasserstein)

- Schatten-2 is the **Frobenius norm**.
- They choose radius `k_2 = ||I||_F = sqrt(d)` so `I ∈ F_2`.

Eliminating `P` yields:
- `max_{π ∈ Π(p,q)} || X π Y^T ||_F` (Eq. 14)

Lemma 4.3 (connection to GW):
- Consider GW with cosine similarities and squared loss.
- Then minimizing GW is equivalent to maximizing `||X π Y^T||_F` over couplings (up to constants), i.e., Eq. 14.

Practical update for `P` given `π` (Frobenius ball):
- maximizing `⟨P, M⟩` s.t. `||P||_F ≤ sqrt(d)` yields `P* = sqrt(d) * M / ||M||_F`.

#### Case C: `p = 1` (nuclear norm)

- Mentioned as encouraging sparse singular spectra (low-rank structure).
- The paper notes details for this case are discussed in the appendix (not in the 10-page main text).

---

## 5. Optimization procedure (Section 4.3)

### 5.1 Alternating maximization

Given the simplified form (Eq. 10), solve by alternating:

1) **Update `P` given `π`:** closed form from Lemma 4.2 (cost is an SVD of a `d×d` matrix → `O(d^3)`).
2) **Update `π` given `P`:** standard OT on a cost matrix induced by `P` (exact LP or entropic Sinkhorn).

### 5.2 Entropic regularization annealing (stabilizing non-convexity)

- Joint problem is not jointly concave → alternating optimization can be initialization-sensitive.
- The paper’s practical fix is to **anneal** the entropic regularization strength:
  - start at a large regularization (smoother objective, less non-concave),
  - decay each iteration: `ε_t = α ε_{t-1}` with `α ∈ (0,1)`,
  - stop at a minimum `ε` and/or when objective converges.
- Reported defaults: `ε_0 = 1`, `α = 0.95`.

### 5.3 Pseudocode sketch (implementation-oriented)

Below is a compact sketch of the alternating method described in the paper (main text; “Algorithm 1” is referenced as being in the supplement).

```text
Inputs:
  X ∈ R^{d×n}, Y ∈ R^{d×m}  (columns are points)
  weights p ∈ Δ^n, q ∈ Δ^m
  Schatten choice p_schatten ∈ {∞, 2, 1, ...} with radius k_p
  entropic regularization ε0, decay α, min ε_min

Initialize:
  P ← I (or any reasonable starting transform)
  ε ← ε0

Repeat until convergence:
  # (A) coupling update (given current P)
  Define cost C_{ij} = ||x_i - P y_j||^2  (or equivalent dot-product form)
  π ← Sinkhorn(p, q, C, ε)

  # (B) transform update (given current π)
  M ← X π Y^T
  Compute SVD: M = U diag(σ) V^T
  Update P:
    if p_schatten = ∞ (orthogonal):
      P ← U V^T
    else if p_schatten = 2 (Frobenius ball radius sqrt(d)):
      P ← sqrt(d) * M / ||M||_F
    else:
      choose s maximizing s^T σ subject to ||s||_p ≤ k_p
      P ← U diag(s) V^T

  # (C) anneal entropy
  ε ← max(ε_min, α * ε)
```

Notes:
- The Sinkhorn step can be replaced by exact OT if desired (slower).
- The transform update always uses only a `d×d` SVD (cheap when `d ≈ 300` for word embeddings).
- The paper mentions a separate incremental/scalable variant for very large problems (Appendix F).

---

## 6. Experiments (Section 5)

### 6.1 Synthetic point clouds

Setup:
- Sample source point cloud `X ⊂ R^3`.
- Sample `P` from an invariance family `F_p` and generate `Y = P X` (plus optional Gaussian noise).
- Evaluate matching accuracy by converting the soft coupling `π` to a hard matching: `j*(i) = argmax_j π_{ij}`.

Findings:
- When the latent transform is orthogonal, adding orthogonal invariance (`p = ∞`) recovers correct correspondences; classic OT can fail (greedy proximity).
- With added noise, the invariance-matched method performs best; they also compare against an “oracle” that knows the transform, and note cases where optimizing over `P` can help mitigate noise.

### 6.2 Unsupervised word translation (bilingual lexicon induction)

Motivation:
- Monolingual embedding spaces have similar relational geometry but are only defined up to global orthogonal transforms → orthogonal invariance is natural.

Dataset:
- Conneau et al. (2018) fastText embeddings + dictionaries; they evaluate En↔{Es, Fr, De, It, Ru}.

Retrieval / evaluation:
- Nearest-neighbor retrieval uses CSLS (cross-domain similarity local scaling).

Results (Table 1; selected highlights):
- Their `p = ∞` invariant OT + CSLS achieves accuracy comparable to strong unsupervised baselines.
- Reported runtime is much lower than adversarial approaches (minutes per language pair on CPU).
- They explicitly report their method *without* the iterative refinement step used in some prior work, so the fairest comparison is against adversarial + CSLS without refinement.

Optimization dynamics (Fig. 4):
- Typical pattern: early iterations mostly adjust `P`, then a phase where both `P` and `π` change strongly and objective drops, then convergence.
- The unsupervised objective correlates with translation accuracy, supporting early stopping/model selection without labeled validation.

---

## 7. Takeaways / how to use this idea

- If you can define a plausible global invariance family `F` (orthogonal, low-rank, etc.), you can “register” spaces *while* computing correspondences via OT.
- Alternating updates are simple and fast:
  - `π` via Sinkhorn,
  - `P` via SVD of `X π Y^T`.
- Entropy annealing is a pragmatic approach to reduce dependence on brittle initialization and hyperparameter sensitivity.
- The `p = 2` specialization provides a clean conceptual bridge to GW: it clarifies when “GW-like” matching emerges from an OT-with-transform view.

