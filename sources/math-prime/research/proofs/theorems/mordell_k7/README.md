# Mordell, k=7: x²+7=y³
Status: LEAN
Methods: [ghost.rational](../../methods/ghost.rational.md) + [norm.pairs](../../methods/norm.pairs.md) + [lifting.somos](../../methods/lifting.somos.md) (prime 5)
Generalizes: none     Generalized by: none. The ghost route for every k>0 stayed open; finiteness for k>0 is
[mordell_positive_k](../mordell_positive_k/README.md), proved by another method.

## Statement
For integers x, y: x²+7=y³ if and only if (x,y) ∈ {(±1,2), (±181,32)}.

## Proof 1
- **Final paper:** [papers/mordell_k7.tex](../../../shared/prior/elementary_diophantine_paper/papers/mordell_k7.tex)
  (with a PDF next to it). The same proof is the section "Seven below a cube" of
  [main.tex](../../../shared/prior/elementary_diophantine_paper/main.tex).
- **Earlier paper version:** [team2/ROOT_CUBIC_PROOF.tex](../../../shared/prior/ghost_parametrization/team_work/team2/ROOT_CUBIC_PROOF.tex),
  "An Elementary Solution of the Mordell Equation x²+7=y³", 16 pages, 24 September 2026.
- **Steps:**
  1. *Ghosts.* The ghosts are the reduced rational points (A,B,C) of Y²=X³−7, and the marker is B=1.
  2. *Parents.* The maps Φ to Y²=X³+189 and Ψ back, with T=ΨΦ, have parents built with `norm.pairs` in two pair rings.
  3. *One index.* The chord S through (2,1) satisfies TS=S³T, so every rational point is exactly one Q_n=Sⁿ(O) with n≠0.
  4. *Marker.* The reduced denominators |h_n| follow a Somos-type recurrence, and 3 divides B_{3m}.
  5. *Exclusion.* A factored stride identity gives Somos-4 subsequences s_k, and 5-adic exact lifting rules out |s_k|=1 for every k≠0.
- **Lean:** `Paper.Seven.integer_solutions`, built in
  [elementary_diophantine_paper](../../../shared/prior/elementary_diophantine_paper/README.md) (209 modules).
  An earlier `CubicFinal.lean` check is recorded in
  [audit_20260927/formal_receipt.json](../../../shared/prior/audit_20260927/formal_receipt.json).
- **Audit:** [audit_20260927/k7_review.md](../../../shared/prior/audit_20260927/k7_review.md).
- **Curve:** y²=x³−7 is LMFDB 5292.i1, with Mordell–Weil group Z. This matches the single index.
