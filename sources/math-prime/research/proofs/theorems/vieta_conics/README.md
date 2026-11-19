# Vieta conics: the k for which x²−kxy+y²+lx+my=0 has positive solutions
Status: WRITTEN-AUDITED (Corollary B passed its review on condition of a repair that was then made, and the repaired
text was not re-reviewed)     Methods: [involution.vieta](../../methods/involution.vieta.md)
Generalizes: Cipu's theorem, the case m = 0 (not a card)     Generalized by: none.
Related: [vieta_gcd_invariant](../vieta_gcd_invariant/README.md) treats the case l = m = r.

## Statement
Let l, m ≥ 0 be integers, not both 0, let k ≥ 3, and put F = x² − kxy + y² + lx + my. Let S(l,m) be the set of
k ≥ 3 for which F = 0 has a solution in positive integers; by Proposition 3.4 of the source, this is the same as
infinitely many. Put

  D1 = {y + m + (l+1)/y : y ≥ 1, y | l+1, y² ≤ l+1},  D2 = {x + l + (m+1)/x : x ≥ 1, x | m+1, x² ≤ m+1}.

**Theorem A.**
(a) k ∈ S(l,m) if and only if some positive (x,y) with F(x,y) = 0 has x² ≤ y² + my and y² ≤ x² + lx. These are
    the lattice points of the arc between the points of horizontal and vertical tangency.
(b) Such a point with x = 1 has k ∈ D1, and one with y = 1 has k ∈ D2. One with x, y ≥ 2 has k ≤ 2 + (l+m)/2,
    with equality only at (2,2).
(c) S(l,m) = D1 ∪ D2 ∪ T with T ⊂ [3, 2 + (l+m)/2], and max S(l,m) = l + m + 2.
(d) Every k ≥ 3 with (k − 2) | (l + m) lies in S(l,m).

**Corollary B (the equation of Miller–Nikolakopoulos–Srinivasan).** For r ≥ 1, (a+r)/b + (b+r)/a = k has a solution
in positive integers if and only if k ∈ S(r,r). The elements above r + 2 are the numbers y + r + (r+1)/y with
y | r+1 and y² ≤ r+1; the elements up to r + 2 form T(r,r), found by the finite search of (a). For example
S(3,3) = {3,4,5,7,8} and S(10,10) = {3,4,5,6,7,10,12,22}.

**Corrected claim.** This determines S(l,m) by a finite search for each (l,m); it is not a closed form, because no
closed form for T is known. The branch's top README (item 5) calls it "the complete set of k", which our review of
2026-10-26 marked OVERSTATED. Effective solvability for every r already follows from Mills (1954, 1956).

**Erratum.** x² − 13xy + y² + 27x = 0 has no solution in positive integers, so 13 ∉ S(27,0). The table of
Feng–Yuan–Hu, INTEGERS 13 (2013) #A8, lists it.

## Proof 1: `involution.vieta`
- **Origin:** branch `origin/claude/math-prime-cross-domain-g1vwsn`, tip `26c718ac`, folder
  `cloud_agent_claude/02_conics_and_diophantine/` (track 02). The files are copied byte for byte into
  [source/](source/).
- **Files:**
  - [T1_fundamental_arc.md](source/02_conics_and_diophantine/theorems/T1_fundamental_arc.md): the arc theorem that
    (a) uses;
  - [T3_lm_family.md](source/02_conics_and_diophantine/theorems/T3_lm_family.md): Theorem 3 (= Theorem A),
    Corollaries 3.1 and 3.2, Proposition 3.4, Corollary 3.5 (= Corollary B), Remark 3.3;
  - [T2_lx_family.md](source/02_conics_and_diophantine/theorems/T2_lx_family.md): the case m = 0, with the erratum
    as Remark 2.4;
  - the scripts `t1_fundamental_arc_check.py`, `t3_lm_family.py`, `t3b_positive_solutions.py`, `t2_lx_family.py`,
    the PARI re-implementation `x1_independent_gp_checks.gp` and the 45 identities of `x2_symbolic_identities.py`,
    each with its output in [evidence/](source/02_conics_and_diophantine/evidence/).
- **Steps:**
  1. The jumps s₁(x,y) = (ky − l − x, y) and s₂(x,y) = (x, kx − m − y) are the affine reflections in the diameters
     F_x = 0 and F_y = 0, which are conjugate to the axes. With N = l² + klm + m² > 0, every orbit of ⟨s₁, s₂⟩
     meets the arc between the points of horizontal and vertical tangency exactly once (T1). This gives (a).
  2. On the arc, 2(k − 2)xy = 2(x − y)² + 2lx + 2my, and one inequality gives (b).
  3. Divisor points (1, y), (x, 1) and diagonal points (t, t) give D1, D2 and (d).
  4. Proposition 3.4: the tangent line at the origin separates the two branches, so every positive solution lies on
     the branch that carries the arc.
  5. Erratum: (2e − 13c)² − 165c² = −108/g ∈ {−108, −36, −12, −4} has no solution modulo 5, 11 or 3. A second proof
     checks the four candidate arc points.
- **Lean:** none.
- **Audit:**
  - The branch's reviews:
    - [02__T1_fundamental_arc.md](source/reviews/02__T1_fundamental_arc.md): PASS, MATCH.
    - [02__T2_lx_family.md](source/reviews/02__T2_lx_family.md): PASS; the erratum CONFIRMED four ways (a
      certificate modulo 165, PARI, brute force to y ≤ 10⁷, the arc check).
    - [02__T3_lm_family.md](source/reviews/02__T3_lm_family.md): PASS, MATCH.
    - [02__T3_T5_revision_items.md](source/reviews/02__T3_T5_revision_items.md): Proposition 3.4 PASS; Corollary B
      PASS once repair R1 was made. R1 was made in the second revision, which was not reviewed again.
  - Our check, 2026-10-26: rev1 judged the README's claim OVERSTATED and did not re-derive the theorem. Its brute
    force for track 13 lists the k for r ≤ 12; they agree with S(r,r) above
    ([t13_general_out.txt](../../../../research/shared/consolidation_2026-10-26/review_checks/rev1/t13_general_out.txt)).
    The consolidation re-checked the erratum by brute force for x ≤ 2·10⁵ and by the arc box
    ([errata_checks.py](../../../../research/shared/consolidation_2026-10-26/a_tracks/work/checks/errata_checks.py)).
- **Provenance and novelty:**
  - Known: Mills (Proc. Amer. Math. Soc. 5 (1954); Canad. J. Math. 8 (1956)) reduces every such conic, with any
    linear terms, to finitely many fundamental points, effectively. Lemmermeyer (arXiv 2601.15229, 2026) uses the
    tangency box in the homogeneous case. Cipu (INTEGERS 15 (2015) #A47) settles m = 0. Miller–Nikolakopoulos–
    Srinivasan (arXiv 2605.19083) settle l = m = r ≤ 2 and write that their descent fails for r ≥ 3.
  - No match for the explicit structure of S(l,m) with m ≥ 1, nor for the k-sets of the MNS equation for r ≥ 3: the
    branch's searches A1–A8, its referees' S1–S12, and ours of 2026-10-26
    ([novelty_web_log.md](../../../../research/shared/consolidation_2026-10-26/a_tracks/work/novelty_web_log.md),
    Q3).
  - The erratum: no published correction of the Feng–Yuan–Hu table was found. Cipu's Proposition 1 already gives
    the correct S(27) without remarking on the discrepancy.

## Open questions
- A closed form for T(l,m). Its elements need not come from divisor or diagonal points: 7 ∈ S(10,3) comes from the
  arc point (2,3).
- Mixed signs lm < 0: the structure (c) is proved only for l, m ≥ 0; for lm < 0 only the finite arc test of T1
  applies (Remark 3.3).
