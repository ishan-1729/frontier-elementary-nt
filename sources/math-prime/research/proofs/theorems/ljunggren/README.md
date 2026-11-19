# Ljunggren: x²+1=2y⁴
Status: LEAN, two independent proofs
Methods: proof 1: [ghost.rational](../../methods/ghost.rational.md) + [lifting.somos](../../methods/lifting.somos.md) (prime 3); proof 2: [approx.normsign](../../methods/approx.normsign.md)
Generalizes: none     Generalized by: [quartic_t_family](../quartic_t_family/README.md), as member t=1 (2x⁴−y²=1)

## Statement
For positive integers x, y: x²+1=2y⁴ if and only if (x,y)=(1,1) or (239,13).

## Proof 1: ghost parametrization and descent (`ghost.rational` + `lifting.somos`)
- **Manuscript:** [ghost_descent/ljunggren_elementary_proof_descent.pdf](ghost_descent/ljunggren_elementary_proof_descent.pdf),
  "Descent and finite congruences for x²+1=2y⁴", 21 September 2026, 15 pages.
- **Source:** [ghost_descent/new_elementary_proof.tex](ghost_descent/new_elementary_proof.tex), a copy of
  [Archive/ljunggren_working_2026_09/research_new/new_elementary_proof.tex](../../../../Archive/ljunggren_working_2026_09/research_new/new_elementary_proof.tex)
  (sha256 2a9e70c4…).
- **Steps:**
  1. *Nodes.* The enlargement is r⁴+q²=2s⁴ with gcd(r,s)=1 and q≡1 (mod 4); the original solutions are the nodes with r=1.
  2. *Descent.* The signed descent uses Pythagoras twice, a²+2b²=c² and the Vierzahlensatz, with height s.
  3. *One index.* The coordinates (−r²/s², rq/s³) lie on Y²=X³−2X. The chord and tangent maps A and D satisfy DA=A²D, so every node is one point P_l.
  4. *Marker.* The Somos-4 sequence b_{l+2}b_{l−2}=1764b_{l+1}b_{l−1}−3107b_l², with seeds (13,1,1,13), carries the reduced
     numerators; the marker is f_l=−1.
  5. *Exclusion.* The reflection b_{−l−1}=b_l preserves f_l, so a marked index can be taken odd, l=2j−1.
     - There the marker reads G(j)=b_{2j}b_{2j−2}−9b_{2j−1}²=−4f_{2j−1}=4.
     - A 3-adic finite congruence proves G(j)=4 only for j∈{0,1}, that is l=±1 (Theorem 5.1 of the manuscript).
     - Over all integer indices, f_l=−1 therefore holds exactly at l∈{−2,−1,0,1}, which are ±1 and their
       reflections. This gives y=|b_l|∈{1,13}.
     - Corrected on 2026-10-26. The step used to read "excludes every index except l=±1", which is true only after
       the reflection to odd indices. The audit in the separate math-treasure repository
       (`research/ghost-classification-audit/README.md`) noted it. Recomputed here: f₋₂=f₋₁=f₀=f₁=−1 and
       f₂=−1803649.
- **Lean:**
  - The entry is [lean/parent_NewElementary.lean](../../../../lean/parent_NewElementary.lean), with theorems
    `Ljunggren.NewElementary.classification`, `classification_iff` and `classification_nat`.
  - The supporting modules are `lean/parent_*.lean`, `lean/grids_Conrad*.lean`, `lean/discrete_PointDomain.lean` and
    `lean/Ljunggren/correspondence_Somos*.lean`.
  - The dependency audit is [lean/parent_NewProofAudit.lean](../../../../lean/parent_NewProofAudit.lean).
- **Verification record:**
  - The manuscript's Appendix B reports a rebuild of 34 local modules and a recursive audit of
    10,856 declarations.
  - These modules lie outside the import graph that `lean/verify.sh` checks (see the
    [README](../../../../README.md) and [ACCEPTED.md](../../../../ACCEPTED.md)).
  - The build was not re-run during the 2026-10-14 reorganization.
- **Working history:** [Archive/ljunggren_working_2026_09/research_new/](../../../../Archive/ljunggren_working_2026_09/research_new/).
- **Provenance:** Appendix B (the descent attributed to Conrad, Elkies 2014, Hone–Swart 2008, Poonen 2014).
  Priority for the complete argument has not been established.

## Proof 2: norm-sign approximation (`approx.normsign`)
- **Exposition** (moved on 2026-10-14 from `results/ljunggren/`): [norm_sign/elementary_proof.md](norm_sign/elementary_proof.md),
  [elementary_norm_sign_proof.md](norm_sign/elementary_norm_sign_proof.md),
  [elementary_norm_sign_lemma.md](norm_sign/elementary_norm_sign_lemma.md),
  [elementary_arithmetic_branch.md](norm_sign/elementary_arithmetic_branch.md),
  [norm_sign_adversarial_audit.md](norm_sign/norm_sign_adversarial_audit.md),
  [LaTeX](norm_sign/ljunggren_elementary_proof.tex), [PDF](norm_sign/ljunggren_elementary_proof.pdf),
  [Overleaf package](norm_sign/ljunggren_overleaf.zip), [check script](norm_sign/verify_elementary_proof.py).
- **Paper version** (Wronskian certificate, finished by one congruence modulo 13):
  [papers/ljunggren.tex](../../../shared/prior/elementary_diophantine_paper/papers/ljunggren.tex).
- **Lean:**
  - [lean/Ljunggren/CompleteNormSign.lean](../../../../lean/Ljunggren/CompleteNormSign.lean) proves
    `Ljunggren.NormSignProof.ljunggren_iff`. It is verified by `lean/verify.sh`; the records are
    [lean/FINAL_AUDIT.md](../../../../lean/FINAL_AUDIT.md) and
    [docs/supervisor/evidence/T006/07_verify.txt](../../../../docs/supervisor/evidence/T006/07_verify.txt).
  - The paper package proves `Paper.Ljunggren.positive_integer_iff`.
- **Acceptance:** [ACCEPTED.md](../../../../ACCEPTED.md) lists it as THEOREM (historical, pre-supervisor),
  with snapshot tag `accepted/ljunggren-v1`.
- **Audit:** [audit_20260927/REPORT.md](../../../shared/prior/audit_20260927/REPORT.md) and
  [ljunggren_review.md](../../../shared/prior/audit_20260927/ljunggren_review.md).

## Remarks
- Squares in the Pell sequence at odd index: H_n²−2P_n²=−1, so P_n=y² gives H_n²+1=2y⁴. Ljunggren
  then leaves P_n∈{1,169}. This corollary is not written up separately.
- The working history of both proofs is in [Archive/ljunggren_working_2026_09/](../../../../Archive/ljunggren_working_2026_09/).
