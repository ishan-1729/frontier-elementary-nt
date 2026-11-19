# Quartic family: (t+1)x⁴ − t·y² = 1
Status: members range from LEAN to OPEN (table below); for t=6 and t=12 the odd-Pell-index half is WRITTEN-AUDITED
(Proof 2)
Methods: Proof 1 (the member table): [ghost.rational](../../methods/ghost.rational.md) +
[lifting.somos](../../methods/lifting.somos.md); Proof 2 (t=6 and t=12, odd Pell index):
[approx.normsign](../../methods/approx.normsign.md)
Generalizes: [ljunggren](../ljunggren/README.md), which is the member t=1     Generalized by: none

## Statement for a member t
Find all positive integer solutions of (t+1)x⁴ − t·y² = 1. The pair (1,1) is always a solution.

The Bennett–Togbé–Walsh subfamily is t = m(m+1). Its members also have the solution (2m+1, 4m²+4m+3).
Their 2006 paper proves that these are all the solutions, using hypergeometric approximation; we treat that
result as a benchmark, not as a premise.

The ghost enlargement is (t+1)B⁴ − t·Q² = C⁴ with (B,C) primitive, and the marker is C=1.

## Members
The members' proofs in the table are ghost descents (Proof 1, in the linked files). Proof 2 below covers one half
of t=6 and t=12.

| t | Equation | Positive solutions | Status | Files | Curve y²=x³−t²(t+1)x (Mordell–Weil group) |
|---|---|---|---|---|---|
| 1 | 2x⁴−y²=1, i.e. y²+1=2x⁴ | (1,1), (13,239) | LEAN | [ljunggren](../ljunggren/README.md) | y²=x³−2x (Z⊕Z/2) |
| 2 | 3x⁴−2y²=1 (Bumby; BTW m=1) | (1,1), (3,11) | LEAN: `TeamThree.BumbyFinal.bumby_positive_iff` | [ROOT_BUMBY_FINAL_RECEIPT.json](../../../shared/prior/ghost_parametrization/team_work/team3/ROOT_BUMBY_FINAL_RECEIPT.json), [expert_bumby_final_classification.lean](../../../shared/prior/ghost_parametrization/team_work/team3/expert_bumby_final_classification.lean) | y²=x³−12x (Z⊕Z/2) |
| 3 | 4x⁴−3y²=1 | (1,1) | WRITTEN-AUDITED | [ROOT_T3_WRITTEN_RECEIPT.json](../../../shared/prior/ghost_parametrization/team_work/team3/ROOT_T3_WRITTEN_RECEIPT.json) | y²=x³−36x (Z⊕(Z/2)²) |
| 6 | 7x⁴−6y²=1 (BTW m=2) | expected (1,1), (5,27) | PARTIAL: odd Pell index WRITTEN-AUDITED (Proof 2); even index OPEN | Proof 2; [handover](../../../shared/prior/ghost_parametrization/ghost_parametrization_handover/START_HERE.md) | y²=x³−252x (Z²⊕Z/2) |
| 7 | 8x⁴−7y²=1 | (1,1) | WRITTEN-AUDITED | [ROOT_T7_WRITTEN_RECEIPT.json](../../../shared/prior/ghost_parametrization/team_work/team3/ROOT_T7_WRITTEN_RECEIPT.json), [ROOT_T7_UNIFORM_ASSEMBLY.md](../../../shared/prior/ghost_parametrization/team_work/team3/ROOT_T7_UNIFORM_ASSEMBLY.md) | y²=x³−392x (Z⊕Z/2) |
| 11 | 12x⁴−11y²=1 | (1,1) | WRITTEN-AUDITED | [ROOT_T11_WRITTEN_RECEIPT.json](../../../shared/prior/ghost_parametrization/team_work/team3/ROOT_T11_WRITTEN_RECEIPT.json) | y²=x³−1452x (Z⊕Z/2) |
| 12 | 13x⁴−12y²=1 (BTW m=3) | expected (1,1), (7,51) | PARTIAL: odd Pell index WRITTEN-AUDITED (Proof 2); even index OPEN | Proof 2; [handover](../../../shared/prior/ghost_parametrization/ghost_parametrization_handover/START_HERE.md) | y²=x³−117x (Z²⊕Z/2) |
| 19 | 20x⁴−19y²=1 | open | OPEN (structural reductions) | [ROOT_T19_PILOT_RECEIPT.json](../../../shared/prior/ghost_parametrization/team_work/team3/ROOT_T19_PILOT_RECEIPT.json) and the other `team3/ROOT_T19_*` receipts | not in LMFDB |
| 26 | 27x⁴−26y²=1 | (1,1) | WRITTEN-AUDITED | [ROOT_T26_WRITTEN_RECEIPT.json](../../../shared/prior/ghost_parametrization/team_work/team3/ROOT_T26_WRITTEN_RECEIPT.json) | y²=x³−18252x (Z⊕Z/2) |

The curve data are from LMFDB, read on 2026-10-14; the labels are on the
[ghost.rational](../../methods/ghost.rational.md) card. Every solved member has rank 1. The open members
t=6 and t=12 have rank 2: the method's single index is not enough for them. This is the generalization target
`ghost.rational.multi_index`.

## Further results on the enlarged equation
- For t=2,3,7,11,26 (with the exception t=2, C=1), each fixed C>0 has at most one primitive positive pair.
  These are written results with constructed divisor sections
  ([ROOT_EFFECTIVE_INDEX_RECEIPT.json](../../../shared/prior/ghost_parametrization/team_work/team3/ROOT_EFFECTIVE_INDEX_RECEIPT.json)).
- For t=2 with C=25 and C=167, the unrestricted classification is proved in Lean
  ([ROOT_BUMBY_UNRESTRICTED_RECEIPT.json](../../../shared/prior/ghost_parametrization/team_work/team3/ROOT_BUMBY_UNRESTRICTED_RECEIPT.json)).
- For t=6, the reductions in the handover show that an additional solution would need W_p to be a square at a
  prime p≥17 with p≡1 or 7 (mod 8). Index multiples of 5 and 7 are excluded
  ([Proof.md](../../../shared/prior/ghost_parametrization/ghost_parametrization_handover/current_sources/ghost_targets_first_obstructions/Proof.md)).

## Proof 2: norm-sign approximation, odd Pell index of t=6 and t=12 (`approx.normsign`)
- **Statement.** Every positive solution of 7x⁴ − 6y² = 1 has a Pell index k ≥ 0, defined by
  x²√7 + y√6 = (√7+√6)^(2k+1). If k is odd, then (x,y) = (5,27). Likewise, with
  x²√13 + y√12 = (√13+√12)^(2k+1), a positive solution of 13x⁴ − 12y² = 1 with k odd is (7,51). The even-index
  halves stay open.
- **Origin:** branch `origin/claude/math-prime-cross-domain-g1vwsn`, tip `26c718ac`, folder
  `cloud_agent_claude/08_analysis_pade_thue/` (track 08 of the cloud Claude programme, 18 October 2026). The files
  are copied byte for byte into [source/](source/).
- **Files:**
  - [01_reduction.md](source/08_analysis_pade_thue/theorems/01_reduction.md): the reduction (Theorem 1);
  - [02_pade_norm_sign.md](source/08_analysis_pade_thue/theorems/02_pade_norm_sign.md): the approximants
    (Theorem 2, Corollary 2.8);
  - [03_caseB_t6.md](source/08_analysis_pade_thue/theorems/03_caseB_t6.md): t=6 (Theorem 3);
  - [06_t12_caseB.md](source/08_analysis_pade_thue/theorems/06_t12_caseB.md): t=12 (Theorem 6);
  - [04_caseA_obstruction.md](source/08_analysis_pade_thue/theorems/04_caseA_obstruction.md): the even index;
  - the scripts s01–s04 and s06–s09, x02 and x03, with outputs in
    [evidence/](source/08_analysis_pade_thue/evidence/); s09 is an independent PARI/GP re-implementation;
  - [01_literature.md](source/08_analysis_pade_thue/notes/01_literature.md), and the 16 failed routes to the even
    index in [03_failed_attempts.md](source/08_analysis_pade_thue/notes/03_failed_attempts.md).
- **Steps** (t=6):
  1. *Reduction* (Theorem 1, the same as Togbé–Voutier–Walsh 2005, Prop. 2.1): Pell parametrisation by descent,
     index doubling and coprime splitting. An odd-index solution gives Q_{d₁,d₂}(g,h) = 1 for one of the four
     pairs (d₁,d₂) with d₁d₂ = 6, and x = d₁g² + d₂h².
  2. *Centring at (5,27)* turns the four quartics into Φ_a(U,W) = c with a = 27√6 and c ∈ {625, 2500}, where
     Φ_a(U,W) = U⁴ − 4aU³W − 6U²W² + 4aUW³ + W⁴ and (U,W) ∈ Z×√6Z or √6Z×Z.
  3. *Approximants* (Theorem 2). The ladder n(n+1)Y_{n+1} = 8an(2n+1)Y_n + 4(16n²−1)Y_{n−1} of the repository's
     norm-sign lemma, with S_n, T_n ∈ Z[a], the sign (−1)ⁿΦ_a(S_n,T_n) > 0 for every real a ≥ 1, the bracket of
     the root ξ_a, and e_n < (1/(4a))(4/a)ⁿ, S_n ≤ (16a + 8/a)ⁿ. At a = 27√6 the forms Λ_n = UT_n − WS_n lie in
     Z ∪ √6Z, so Λ_n ≠ 0 implies |Λ_n| ≥ 1, and Λ_n, Λ_{n+1} never vanish together.
  4. *Finite analysis.* Large U by the uniform bounds, the middle range by exact first approximants and convexity,
     so U < 11; a search of the 8 admissible pairs (g,h) with U ≤ 12 leaves only (d₁,d₂,g,h) = (2,3,1,1). For
     t=12: a = 102√3, c = 2401, U < 14, and a search of the 12 admissible pairs with U ≤ 15.
- **What is new relative to this repository.** The norm-sign lemma for every real a ≥ 1 is already the
  repository's ([elementary_norm_sign_lemma.md](../ljunggren/norm_sign/elementary_norm_sign_lemma.md), eq. (2)).
  New are: the integrality S_n, T_n ∈ Z[a] by Lagrange inversion and a second-order ODE (the repository proved
  integrality at a = 239 by a separate construction); the explicit bounds of Corollary 2.8; the centring at an
  irrational parameter; and the finite analysis of Theorems 3 and 6.
- **Lean:** none.
- **Audit:**
  - The branch's reviews, all PASS with fidelity MATCH:
    [08__T1_reduction.md](source/reviews/08__T1_reduction.md),
    [08__T2_norm_sign.md](source/reviews/08__T2_norm_sign.md),
    [08__T3_caseB.md](source/reviews/08__T3_caseB.md) (every constant decided exactly in Q(√6); the tightest is
    f₃(968) = 0.99711… < 1), [08__T6_t12.md](source/reviews/08__T6_t12.md), and
    [08__P4_obstruction.md](source/reviews/08__P4_obstruction.md) (Prop. 4.4 after a required repair, which was made).
  - Our check, 2026-10-26 (rev1): CONSISTENT. It re-verified the identities of Theorems 1 and 3 and the constants of
    §3.5, and reproduced the index sieve
    ([t08_checks.py](../../../../research/shared/consolidation_2026-10-26/review_checks/rev1/t08_checks.py),
    [t08_out.txt](../../../../research/shared/consolidation_2026-10-26/review_checks/rev1/t08_out.txt)). It found one clause of
    the branch README OVERSTATED: that this "extends the repository's own Ljunggren norm-sign lemma to every real
    parameter a ≥ 1".
- **Provenance and novelty:**
  - Both statements are known: Bennett–Togbé–Walsh (BTW), Int. J. Number Theory 2 (2006) 195–206, Theorem 1.2 with
    m = 2, 3.
  - t=6: BTW used Baker-type linear forms in logarithms (KANT, ThueSolve), and Togbé–Voutier–Walsh (2005) wrote of
    this index class "We do not seem to be able to apply the hypergeometric method in this case". No analytic proof
    was found: the branch's searches, and ours
    ([novelty_web_log.md](../../../../research/shared/consolidation_2026-10-26/a_tracks/work/novelty_web_log.md),
    Q6). So this is a new proof of a known result, by an approximation method that the owner counts as classical
    in spirit ([approx.normsign](../../methods/approx.normsign.md), owner's note).
  - t=12: BTW's hypergeometric Theorem 3.7 (m ≥ 3), with a Runge-type computation, already proves the odd half.
    Proof 2 stays in the same Thue–Siegel family (its approximants agree numerically with the diagonal Padé
    approximants of ((a+i)/(a−i))^(1/4)), with the norm-sign bracket in place of a remainder estimate and no
    Runge step. It is a variant, not a new mechanism.

## Open questions
1. **The even Pell index of t=6** (where Proof 2 stops):
   - An even-index solution gives coprime (u,v), v ≥ 1, 0 < u < √6·v, d ∈ {1,2,3,6}, d | u, with
     u⁴ + 24u³v − 36u²v² − 144uv³ + 36v⁴ = d², and u/v is a convergent of the root β₃ ≈ 0.2381 with
     |u/v − β₃| < 1.35/v⁴ (Theorem 4.1, Prop. 4.2).
   - Index sieve (Prop. 4.3; 22 moduli; no reciprocity): if X_k is a square, then k ≡ 0, 1, −1, −2 (mod 720720).
     An even-index solution other than (1,1) therefore has k ≡ 0 or −2 (mod 720720), k ≥ 720718 and
     x > 10^509666.
   - The obstruction (Prop. 4.4 and §4.3): a base point for β₃ with E > 1, κ < 4 and y₀ ≥ 22 is a convergent of
     β₃. Among all base points with y₀ ≤ 2·10⁵, and all convergents up to denominators near 10^1057, only 5/21
     qualifies, with κ = 3.381 (3.045 with Voutier's constants, 3.026 with the best arithmetic), while κ < 3 is
     needed. This part is a measurement, not a theorem.
2. The even Pell index of t=12: no base point with κ < 23 for y₀ ≤ 5·10⁴ (a measurement).
3. t=2 (Bumby): κ ≈ 4.06 even for the odd half, so Proof 2's method does not apply.
4. Prove that the approximants are the diagonal Padé approximants of ((a+i)/(a−i))^(1/4) (proposed route:
   g₁ + ig₂ = H/(1 − iH)).
5. The whole BTW family by Proof 2's method, as planned in the branch's
   [02_exploration_log.md](source/08_analysis_pade_thue/notes/02_exploration_log.md) §10 (unverified).

## Literature
- Ljunggren (1942): X⁴−DY²=1 has at most two positive solutions, proved with quartic fields and Skolem's method.
- Bennett–Walsh and Akhtari: aX⁴−bY²=1.
- Togbé–Voutier–Walsh, Acta Arith. 120 (2005): the even-index class, analytically for t ≥ 204 and by computation
  (SIMATH, KANT) below that.
- Bennett–Togbé–Walsh (2006): the family t=m(m+1).
- Bumby (1967): t=2, by arithmetic in Q(√−2,√−3).
- Walsh (2026, arXiv:2603.05809; Funct. Approx. Comment. Math., 2026): an elementary proof for t=2 that uses Jacobi
  symbols.
