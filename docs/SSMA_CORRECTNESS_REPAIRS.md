# SSMA numerical correctness repairs (2026-11-03)

These changes repair the existing objective and API behavior. They do not establish useful cross-task alignment or qualify historical scientific outputs under the repaired implementation.

- Domain-score optimization removes feasible-set constants before scaling and testing the tangent gradient. Positive rescaling of the cost and block-scalar identity shifts no longer turn a nonstationary solution into immediate convergence. Original objective values, removed constants, and scaling are reported.
- Full-rank domain blocks receive exact conditional orthogonal updates, including reflections, before stationarity is accepted. This addresses disconnected orientation components. General multi-domain global optimality remains unproven.
- Global SSMA retains valid zero-cost eigenvectors. The old blanket near-zero filter could remove perfectly aligned nonconstant scores. There is no implicit constant-mode deflation. eig_tol now bounds materially negative numerical eigenvalues instead of discarding zero-cost directions.
- Singleton label sampling preserves the actual row index. Explicit indices must be finite positive in-range integers connecting different domains. Invalid links raise errors instead of being truncated or removed. Nonnegative finite weights retain the existing zero-weight omission policy.
- The shared safe_compute wrapper preserves successfully computed values after warnings. A custom warning callback can handle and muffle the warning without unwinding the computation.
- Held-out retrieval uses the original target-domain IDs and counts partial missing/invalid predictions as misses. An entirely unavailable result remains NA.
- The PR evidence runner reports missing mandatory checks even when another test fails without a test title. Those failures remain critical.

Regression coverage includes public API witnesses, cost-scale and shift invariance, the exact two-domain full-rank Procrustes oracle over 20 fixed seeds, exact-correspondence zero modes in both eigen backends, label/index integrity, warning recovery, and held-out denominators. The new SSMA regression cases are mandatory PR evidence.

Historical fits must retain their original implementation provenance. Recomputing them under this implementation is a new analysis run. Weak-anchor selection remains a concern about the objective, and coordinate-dependent WM noise/pooling remains a downstream statistical issue; neither is changed here.
