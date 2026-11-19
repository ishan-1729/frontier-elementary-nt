# Mordell, k>0: finiteness of x²+d=y³
Status: LEAN     Methods: [approx.binomial](../../methods/approx.binomial.md)
Generalizes: none     Generalized by: with [mordell_negative_k](../mordell_negative_k/README.md), finiteness for every k≠0

## Statement
For each fixed integer d>0, the equation x²+d=y³ has finitely many integer solutions.

## Proof 1 (`approx.binomial`)
- **Written proof:** [plus_k_completion/PLUS_K_PROOF.md](../../../shared/prior/ghost_parametrization/team_work/plus_k_completion/PLUS_K_PROOF.md).
  The adversarial audit is [AUDIT.md](../../../shared/prior/ghost_parametrization/team_work/plus_k_completion/AUDIT.md).
- **Lean:** `PlusK.mordell_plus_finite` in
  [PlusKComplete.lean](../../../shared/prior/ghost_parametrization/team_work/plus_k_completion/PlusKComplete.lean).
  - The audit file is [PlusKFinalAudit.lean](../../../shared/prior/ghost_parametrization/team_work/plus_k_completion/PlusKFinalAudit.lean).
  - The record is [FINAL_VERIFICATION.json](../../../shared/prior/ghost_parametrization/team_work/plus_k_completion/FINAL_VERIFICATION.json),
    dated 2026-09-25: 31 new modules, standard axioms only.
- **Provenance:** the approximants are classical binomial approximants, and no novelty is claimed.

## Open questions
- A ghost-method proof for k>0. [ghost.unit](../../methods/ghost.unit.md) needs a single unit, but here the
  cubic is totally real. Draft: [PlusK_Ghost_Descent_Arithmetic_Edition_Tangent_Revision.tex](../../../shared/prior/ghost_parametrization/team_work/pro4_artifacts/PlusK_Ghost_Descent_Arithmetic_Edition_Tangent_Revision.tex);
  the state of the problem is in [CURRENT_STATE.md](../../../shared/prior/ghost_parametrization/team_work/CURRENT_STATE.md).
