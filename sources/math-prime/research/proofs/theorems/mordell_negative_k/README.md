# Mordell, k<0: finiteness of x²+k=y³
Status: LEAN
Methods: [ghost.unit](../../methods/ghost.unit.md) + [lifting.binomial](../../methods/lifting.binomial.md)
Generalizes: none     Generalized by: with [mordell_positive_k](../mordell_positive_k/README.md), finiteness for every k≠0

## Statement
For each fixed integer k<0, the equation x²+k=y³ has finitely many integer solutions.
Equivalently, C²=A³+d with d=−k>0 and (A,C)=(y,x).

## Proof 1
- **Paper:** [papers/mordell_negative_k.tex](../../../shared/prior/elementary_diophantine_paper/papers/mordell_negative_k.tex),
  and the section "Negative Mordell coefficient" of [main.tex](../../../shared/prior/elementary_diophantine_paper/main.tex).
- **Earlier drafts:**
  - [MinusK_Ghost_Descent_Arithmetic_Edition_Tangent_Revision.tex](../../../shared/prior/ghost_parametrization/team_work/pro4_artifacts/MinusK_Ghost_Descent_Arithmetic_Edition_Tangent_Revision.tex);
  - the written reconstruction [team2/ROOT_POSITIVE_MORDELL_FINITE_WRITTEN.md](../../../shared/prior/ghost_parametrization/team_work/team2/ROOT_POSITIVE_MORDELL_FINITE_WRITTEN.md),
    where "positive c" means Y²=X³+c.
- **Lean:** `Paper.Minus.negative_k_finite`, with the return lemmas `Paper.ZeroTwo.sequence_zero_bound_small`
  and `Paper.ZeroTwo.matrix_zero_bound_small`, in
  [elementary_diophantine_paper](../../../shared/prior/elementary_diophantine_paper/README.md).
- **Audits:**
  - [audit_20260927/negative_k_review.md](../../../shared/prior/audit_20260927/negative_k_review.md), which reviews the draft;
  - [AUDIT_AND_OPTIMIZATIONS.md](../../../shared/prior/elementary_diophantine_paper/AUDIT_AND_OPTIMIZATIONS.md).
- **Scope:** finiteness, plus a count of returns on each orbit. It gives no algorithm and no bound on positions.
- **Provenance:** the team's own audit calls it an elementary reconstruction of Skolem's mechanism (see
  [ghost.unit](../../methods/ghost.unit.md)).
