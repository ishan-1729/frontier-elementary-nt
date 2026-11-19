# approx.binomial: binomial (Padé) approximants
Parent: approx     Children: none

## What it does (Mordell, k>0)
1. Pigeonhole and Bézout place every solution of x²+d=y³ in one of finitely many cubic equations
   with a fixed coordinate.
2. In each equation, every large solution is compared with one fixed solution.
3. Explicit binomial approximants then produce two integer pairs, at least one of them nonzero. If the
   solution were too large, both pairs would have squared length below 1, which is impossible.

## Theorems using it
[mordell_positive_k](../theorems/mordell_positive_k/README.md).

## Provenance
The approximants are classical (Padé or hypergeometric). The proof establishes their properties directly
and makes no novelty claim; see
[PLUS_K_PROOF.md](../../shared/prior/ghost_parametrization/team_work/plus_k_completion/PLUS_K_PROOF.md).
