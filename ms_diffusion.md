# Multiscale manifold alignment (Wang & Mahadevan) — diffusion wavelets scalability notes

## Part 1 — replacing diffusion-wavelet QR + avoiding expensive operator construction

### Problem focus: the real bottleneck

Diffusion Wavelets (as used in Wang & Mahadevan multiscale manifold alignment) builds a hierarchy by repeatedly compressing dyadic powers of an operator using a *modified QR* routine (paper Fig. 2). A naïve/full QR on the dense operator blocks is prohibitive.

Key observation: diffusion wavelets does **not** need a full QR factorization. At each scale it only needs:

- an orthonormal basis for the **column space / range** of the operator power (to tolerance `ε`), and
- a small **compressed** representation of the operator in that basis.

So: replace “big QR” with any scalable procedure that yields (i) an approximate range basis and (ii) a small projected operator.

---

## 1) What the QR is actually doing (and what we can replace)

In the diffusion-wavelets loop (paper Fig. 2), the QR subroutine is used as a **rank-revealing, ε-truncated orthogonalization**.

- Input: a matrix (block)  
  \(A = [T^{2^j}]_{\phi_j\phi_j}\)
- Output: a basis \(Q\) whose columns span \(\mathrm{range}(A)\) up to precision `ε`  
  (keep directions that are not numerically dependent / not below `ε`)
- Plus a compressed operator representation in that basis, e.g.  
  \(B_j = Q_{j+1}^T\, T^{2^j}\, Q_{j+1}\)

Conceptually the algorithm builds nested subspaces:

- \(V_{j+1}\subseteq V_j\) (nested)
- \(p_{j+1}\ll p_j\) in real data: dimensions collapse quickly (typically <10 levels)

Therefore, any method that reliably produces:

- an orthonormal basis \(Q_{j+1}\) for \(\mathrm{range}(T^{2^j})\) (or the appropriate transformed operator), and
- a small projected operator \(B_j\),

is a viable substitute for the “modified QR”.

---

## 2) Major algebraic speedup: rewrite \(T\) to avoid pseudoinverses + huge feature-space products

The paper defines (Algorithm step 1):

\[
T = F^{+}\, Z L Z^T\, (F^T)^{+}
\]

with \(FF^T = Z D Z^T\), and \(F\) is constructed via SVD.

Let:

- \(A := Z D^{1/2}\)
- thin SVD: \(A = U\Sigma V^T\)

Then:

- \(Z D Z^T = U\Sigma^2 U^T\)
- choose \(F = U\Sigma\)
- pseudoinverses:
  - \(F^{+} = \Sigma^{-1}U^T\)
  - \((F^T)^{+} = U\Sigma^{-1}\)

Also:

- \(Z = A D^{-1/2} = U\Sigma V^T D^{-1/2}\)

Compute the “pseudoinverse sandwich” terms:

- \(F^{+}Z = \Sigma^{-1}U^T(U\Sigma V^T D^{-1/2}) = V^T D^{-1/2}\)
- \(Z^T(F^T)^{+} = (D^{-1/2}V\Sigma U^T)(U\Sigma^{-1}) = D^{-1/2}V\)

So:

\[
\boxed{\,T = V^T\, D^{-1/2}\, L\, D^{-1/2}\, V\,}
\]

Why this matters:

- avoids explicitly forming \(ZLZ^T\) in feature space,
- avoids explicitly forming pseudoinverses,
- makes it natural to apply \(T\) **matrix-free**, leveraging sparse \(L\) (kNN graph Laplacian) and dense/small \(V\).

### Matrix-free matvec for \(T\)

To compute \(y = T x\):

1. \(u \leftarrow Vx\)
2. \(u \leftarrow D^{-1/2}u\)
3. \(u \leftarrow L u\)  (sparse matvec)
4. \(u \leftarrow D^{-1/2}u\)
5. \(y \leftarrow V^T u\)

This is exactly the structure you want for randomized range finding and iterative eigensolvers.

---

## 3) Replace diffusion-wavelet QR with scalable equivalents

### Option A — implement ε-truncated / rank-revealing QR (not full QR)

Interpret the “modified QR” as:

- incremental Gram–Schmidt (or Householder),
- drop candidate columns whose residual norm after orthogonalization is `< ε`,
- optionally with pivoting (choose the “best” next column).

This reduces cost from “full QR” to “build a basis of size \(p_{j+1}\)”, roughly scaling with the retained rank rather than \(n^3\).

Caveat: still assumes efficient access to **columns** (or blocks) of the operator. For a huge dense \(T\), column access is still expensive unless you have fast operator structure.

### Option B — randomized range finding (usually best for big data)

At each scale you need a basis for \(\mathrm{range}(A)\). Randomized linear algebra does exactly that:

1. Draw \(\Omega \in \mathbb{R}^{n\times s}\) (Gaussian or Rademacher), with sketch size \(s \approx 2\!-\!4\times\) expected rank
2. Form \(Y = A\Omega\) using only matvecs / block-matvecs
3. Thin QR: \(Y = QR\)
4. \(Q\) approximates \(\mathrm{range}(A)\)

Why it helps:

- QR is only on a tall-skinny \(n\times s\) matrix with \(s\ll n\).
- Heavy lifting is matvecs, which can use the matrix-free \(T = V^T D^{-1/2} L D^{-1/2} V\) chain.

#### Preserving multiscale / nested structure

Enforce nesting via *small* compressed operators:

- At scale \(j\), have \(Q_j\) spanning the current subspace.
- Form the compressed operator:
  \[
  B_j = Q_j^T\, A_j\, Q_j
  \]
  where \(A_j\) represents the dyadic power (or diffusion operator at that scale).
- Compute eig/SVD of \(B_j\) (small \(p_j\times p_j\)).
- Keep components above threshold → define \(p_{j+1}\).
- Update:
  \[
  Q_{j+1} = Q_j U_j
  \]

This matches the diffusion-wavelets “compress, square, compress…” spirit.

#### Adaptive rank selection

Avoid guessing \(p_j\):

- Use a few random test vectors \(g\) to estimate:
  \[
  \frac{\|(I-QQ^T)Ag\|}{\|Ag\|}
  \]
- If too large: increase sketch size \(s\) and/or add 1–2 subspace iterations.

Rough scaling:

- \(O(\text{cost}(A\cdot\text{block})\times s) + O(ns^2)\)  
  instead of \(O(n^3)\).

### Option C — PSD case: pivoted Cholesky / Nyström

If the operator is symmetric PSD (paper: \(L\) PSD; \(T\) PSD under the construction), use rank-revealing PSD approximations:

- pivoted / incomplete Cholesky with pivoting, or
- Nyström / landmark column sampling.

These yield \(A \approx GG^T\), from which an orthonormal basis can be derived cheaply.

Most applicable when:

- you work with a symmetric version of the operator,
- you can access diagonal / selected columns (or approximate them).

### Option D — skip diffusion wavelets: use eigen-subspace equivalence (symmetric case)

Paper Theorem 3 (symmetric case): diffusion-wavelet scaling functions at level \(k\) span the same space as the eigenvector subspace used by manifold projections (up to rotation).

Engineering implication:

- compute the needed part of the spectrum of \(T\) with an iterative solver (Lanczos / LOBPCG),
- define “scales” by eigenvalue thresholds / spectral gaps (coarse → fine).

This eliminates QR-in-the-loop entirely and uses only matvecs.

### Option E — graph-native multiscale: algebraic multigrid (AMG)

Diffusion wavelets is building a multiresolution hierarchy for a diffusion-like operator. AMG is a mature scalable alternative on sparse graph Laplacians:

- build coarse levels (graph coarsening),
- use AMG prolongation/restriction as scaling bases,
- compute alignment at each coarse level.

Best when the joint Laplacian is huge and sparse (QR is dead on arrival).

---

## 4) Practical implementation plan (scalable path)

### Phase 0 — choose the hierarchy strategy

Rule of thumb:

- \(N \lesssim 20\text{k}\), moderate density → randomized range finder (Option B)
- huge + sparse graph Laplacian → AMG (Option E)
- symmetric PSD and simplest pipeline → eigen/LOBPCG + spectral-gap scales (Option D)
- PSD kernel-like + diagonal/column access → pivoted Cholesky/Nyström (Option C)

Default starting point if unsure: **Option B** (drop-in replacement for QR with matrix-free matvecs).

### Phase 1 — build sparse graphs (often the real cost driver)

Keep the joint Laplacian \(L\) sparse:

1. Build kNN graphs for \(X\) and \(Y\) (typical \(k \sim 10\!-\!50\))
2. Symmetrize weights (mutual kNN or \(W\leftarrow (W+W^T)/2\))
3. Construct \(L_x, L_y\)
4. Add cross-domain correspondences via the \(\Omega\) blocks (paper’s \(L\) construction)

Goal: \(\mathrm{nnz}(L) = O(k(m+n))\), not \(O((m+n)^2)\).

### Phase 2 — build/apply \(T\) without forming huge matrices (use the SVD rewrite)

Compute a thin / randomized SVD of:

- \(A = Z D^{1/2}\)

Truncate aggressively:

- keep singular values above a tolerance, or
- keep top \(\tilde r\).

Then form/apply:

- \(T = V^T D^{-1/2} L D^{-1/2} V\)

Depending on \(\tilde r\):

- form \(T\) explicitly if \(\tilde r\) is a few thousand, or
- keep \(T\) matrix-free if \(\tilde r\) is larger.

### Phase 3 — replace diffusion-wavelets QR with randomized nested compression (recommended)

Instead of “QR(T)” at the finest level:

- choose sketch size \(s\) (start 128–512),
- compute \(Y = (T^\star)\Omega\), where \(T^\star\) is the operator you waveletize.

#### Avoid instability from \(T^+\)

The paper runs diffusion wavelets on \(T^+\). In practice, waveletize a spectrally transformed operator that emphasizes small eigenvalues of \(T\) without explicitly forming \(T^+\), e.g.:

- \(M = (I + \tau T)^{-1}\)  (shift-invert diffusion-like), or
- \(M = I - \alpha T\) with \(\alpha \le 1/\lambda_{\max}\) (keeps spectrum in \([0,1]\) when PSD).

Then perform DWT-like squaring/compression on \(M\): large eigenvalues of \(M\) correspond to small eigenvalues of \(T\).

#### Nested level updates

At each level \(j\):

1. Form small \(B_j = Q_j^T M^{2^j} Q_j\) (or update by squaring the small matrix)
2. Eig/SVD of \(B_j\)
3. Keep components above `ε` → defines \(p_{j+1}\)
4. \(Q_{j+1} = Q_j U_j\)

All “heavy” decompositions are then on small \(p_j\times p_j\) matrices.

### Phase 4 — compute mapping functions without materializing huge factors

Paper’s mapping functions:

\[
\begin{bmatrix}\alpha_k\\ \beta_k\end{bmatrix} = (F^T)^+ [\phi_k]_{\phi_0}
\quad\text{with}\quad
(F^T)^+ = U\Sigma^{-1}
\]

If \(U\) is too large to store:

- from \(A = Z D^{1/2} = U\Sigma V^T\) we have \(U = A V \Sigma^{-1}\)
- so:
  \[
  U\Sigma^{-1} = A V \Sigma^{-2} = (Z D^{1/2})\, V\, \Sigma^{-2}
  \]

Thus a practical mapping computation:

\[
\boxed{\;
\Gamma_k :=
\begin{bmatrix}\alpha_k\\ \beta_k\end{bmatrix}
= (Z D^{1/2})\, V\, \Sigma^{-2}\, Q_k
\;}
\]

Then split:

- \(\alpha_k\): top \(p\) rows
- \(\beta_k\): bottom \(q\) rows

This is attractive when \(Z\) is sparse (e.g., one-hot / BoW); everything after \(Z\) is dense but smaller.

### Phase 5 — validation and selecting a useful scale

Practical level selection criteria:

- held-out correspondences: nearest-neighbor accuracy in latent space,
- stability across random seeds / sketches,
- spectral-gap heuristics (especially if using eigen-based scales).

In practice, only a handful of levels tend to be meaningful.

---

## 5) “Speedup stack” summary (first things to do)

Most actionable sequence:

1. Sparse kNN graphs → sparse \(L\)
2. Thin/randomized SVD of \(A = Z D^{1/2}\)
3. Use \(T = V^T D^{-1/2} L D^{-1/2} V\) (avoid pseudoinverse products)
4. Replace diffusion-wavelets QR with randomized range finding + nested compression
5. Compute mappings via \(\Gamma_k = (Z D^{1/2}) V \Sigma^{-2} Q_k\) (avoid storing \(U\))

---

## Part 2 — complete scalable algorithm (single/multiset) + QR-free implementations

### What changes from the paper implementation

Paper pipeline (Wang & Mahadevan) is:

1. Build joint \(L, D, Z\)
2. Build \(T = F^+ Z L Z^T (F^T)^+\)
3. Run DWT on \(T^+\) (modified QR at each scale)
4. Compute mappings via \((F^T)^+ \Phi_k\)

Scalable replacement keeps steps (1)(2)(4), but replaces step (3):

- **Method A (preferred for symmetric PSD):** iterative eigensolver hierarchy (Theorem 3-equivalent subspaces)
- **Method B:** randomized/sketched DWT (thin QR on sketches only)

### End-to-end scalable pipeline (2-set or S-set)

1) **Within-domain sparse graphs**

- Build symmetric kNN per domain.
- Use sparse storage (CSR/CSC).
- Keep \(\mathrm{nnz}(L)\) near \(O(kN)\).

2) **Joint graph/Laplacian**

- Assemble block \(W_{\text{joint}}\):
  - diagonal blocks = within-domain \(W_s\),
  - off-diagonal blocks = sparse correspondences \(C_{st}\) with weights \(\mu_{st}\).
- Build \(L = D_{\text{joint}} - W_{\text{joint}}\).

3) **Operator compression core**

- For each domain \(s\), compute truncated SVD of \(X^{(s)} D_s^{1/2}\):
  \[
  X^{(s)} D_s^{1/2} \approx U_s \Sigma_s V_s^T
  \]
- Build block-diagonal \(V = \mathrm{blkdiag}(V_1,\dots,V_S)\), \(D^{-1/2}\).
- Use compressed operator:
  \[
  T = V^T D^{-1/2} L D^{-1/2} V
  \]
- Apply \(T\) matrix-free with sparse matvec in the middle.

4) **Multiscale basis without full QR**

- **Method A (spectral):**
  - solve \(T u_i = \lambda_i u_i\) for small non-zero \(\lambda_i\),
  - define dyadic scales \(t_j\), keep vectors by \(\exp(-t_j \lambda_i)\ge \varepsilon\),
  - obtain nested \(\Phi_j\) by truncation.
- **Method B (randomized DWT):**
  - use \(A_0 \approx (T+\delta I)^{-1}\) (regularized inverse),
  - range sketch \(Y=A(\Omega)\), thin QR on \(Y\),
  - compress to small \(B\), eig-truncate, square in reduced space, repeat.

5) **Mapping extraction**

- Paper form stays:
  \[
  \Gamma_k = (F^T)^+ \Phi_k
  \]
- Practical low-rank equivalent:
  \[
  \Gamma_k = (Z D^{1/2}) V \Sigma^{-2} \Phi_k
  \]
- Split \(\Gamma_k\) by domain block for domain-specific maps.

### Multiset-specific construction notes

- Replace 2-block structures by \(S \times S\) blocks:
  - \(Z = \mathrm{blkdiag}(X^{(1)},\dots,X^{(S)})\),
  - \(D = \mathrm{blkdiag}(D_1,\dots,D_S)\),
  - \(L\) from within + cross-domain correspondence blocks.
- Keep cross-domain edges sparse (anchor/star topology preferred over dense all-pairs).
- Complexity is controlled by chosen ranks \(r_s\), not raw feature dims.

### Practical defaults for implementation

- Start with **Method A** (`backend="spectral"`) for symmetric graphs.
- Use moderate per-domain truncation (`r_s` 32–256 initially).
- Use kNN \(k\) in 10–30 range, symmetrized heat weights.
- For very large \(N\), prioritize:
  - approximate kNN construction,
  - block matvecs for dense feature multiplies,
  - avoiding any explicit dense \(N\times N\) matrices.

---

## Part 1A - codebase conformity notes (implementation planning)

### Existing interface patterns to follow

- Public method entry points are S3 generics with `hyperdesign`/`list`/`default`
  methods (e.g., `kema`, `spectral_mnn_align`, `global_geo_align`).
- The package uses an adapter contract for pairwise/multiset composition:
  `new_aligner()`, `fit_pair()`, `fit_many()`, `relative_transform()`,
  `pair_loss()`, and `latent_dim()`.
- Multi-domain native methods are wrapped by adapter constructors in
  `R/adapters_multi.R` (e.g., `kema_aligner`, `spectral_mnn_aligner`).

### Return object conventions

- Preferred result shape is a `multiblock_biprojector` with an additional
  subclass (e.g., `"spectral_mnn_align"`, `"global_geo_align"`).
- Shared helper `new_alignment_result(...)` standardizes this structure and
  attaches method-specific extras.
- Core fields expected across methods:
  - `s` (scores / embedding),
  - `v` (loadings / primal map),
  - `sdev`,
  - `preproc`,
  - `block_indices`,
  plus method extras (`rotations`, `anchors`, `eigenvalues`, etc.).

### Control object conventions

- New configurable methods typically define:
  - `foo_control(...)` returning a typed list (`class = "foo_control"`),
  - `resolve_foo_control(control)` for defaults + validation + named-list merge.
- Validation style is strict and explicit (`stop(..., call. = FALSE)`), with
  unknown control fields rejected.

### Linear algebra / scalability patterns already present

- Matrix-free operator pattern is already used in `R/kema_scalable.R`:
  operators are represented as lists with:
  - dimensions (`n`, `r`),
  - forward matvec `mv(x)`,
  - transpose matvec `tmv(y)`.
- Block/batched gram builders (`gram_ZZ`, `gram_Z_A_Z`, `gram_Z_diag_Z`) avoid
  materializing huge matrices.
- Pivoted Cholesky helper (`pivoted_chol_kernel_block`) already exists and can
  support PSD low-rank alternatives.
- `lowrank_align` already supports explicit vs operator solver backends, which
  is a strong precedent for adding a solver backend switch for multiscale KEMA.

### Reusable spectral/diffusion utilities

- `R/structure_descriptors.R` provides reusable, tested primitives:
  - Laplacian construction (`compute_graph_laplacian`),
  - eigensolver wrapper (`compute_laplacian_basis`),
  - diffusion coordinates (`compute_diffusion_coordinates`),
  - HKS descriptors.
- This suggests implementing multiscale hierarchy code as an internal operator
  module, while reusing existing spectral helper style and validation.

### Recommended integration shape for multiscale implementation

- Keep `kema(...)` API stable; add optional control/backend arguments rather than
  replacing current defaults.
- Introduce a control object (tentative):
  `kema_multiscale_control(...)` with fields like:
  - `backend = c("exact", "randomized", "eigen", "amg")`,
  - `sketch_size`, `oversample`, `n_iter`,
  - `eps_rank`, `max_levels`,
  - `operator_transform` (`"shift_invert"` / `"linear"`),
  - `sparse_knn_k`, `verbose`.
- Add backend dispatch inside KEMA solver stage:
  - baseline exact path unchanged,
  - multiscale path uses matrix-free operator and nested compression.
- Keep output class as `multiblock_biprojector` + subclass extras
  (e.g., hierarchy metadata, per-level dimensions, diagnostics).

### Test alignment with existing style

- Follow current testthat expectations:
  - class checks for subclass + `multiblock_biprojector`,
  - finite outputs and dimension checks,
  - backend consistency checks on toy data (as done in low-rank/operator tests).
- Add small smoke tests first, then numerical-regression tests behind optional
  benchmark gates where needed.
