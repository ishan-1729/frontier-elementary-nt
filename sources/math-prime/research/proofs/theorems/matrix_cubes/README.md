# Cubes in M₂(Z): every 2×2 integer matrix is a sum of three cubes
Status: Theorems A and B LEAN; Theorem C WRITTEN-AUDITED (computer-assisted at p = 2, 3); two-cube search over Z
computational; two conjectures (table below)
Methods: [matrix.cayley_hamilton](../../methods/matrix.cayley_hamilton.md)
Generalizes: none     Generalized by: none

## Statement
A matrix A ∈ M₂(Z) is a *sum of s cubes* if A = X₁³ + … + X_s³ with every Xᵢ ∈ M₂(Z).

**Theorem A (three cubes).** Every A ∈ M₂(Z) is X³ + Y³ + Z³ with X³ = X, Y³ = Y and tr Z ∈ {0, ±1, ±2}. The
matrix 3I is not a sum of two cubes, and it is not a cube (det 3I = 9). So the least s for which every element of
M₂(Z) is a sum of s cubes is 3.

**Theorem B (two-cube obstruction).** If A = X³ + Y³ with X, Y ∈ M₂(Z) and every entry of A is divisible by 3,
then 9 | tr A. The same holds in M₂(Z/9Z).

**Theorem C (local theorem).** Put O = {A : A ≡ 0 (mod 3), tr A ≢ 0 (mod 9)}. A matrix A ∈ M₂(Z) is a sum of two
cubes in M₂(Z_p) for every prime p, and in M₂(R), if and only if A ∉ O.

| Result | Status |
|---|---|
| Theorem A, including w₃(M₂(Z)) = 3 | LEAN: the headline in the author's file, the refinement X³ = X, Y³ = Y, tr Z ∈ {0, ±1, ±2} only in the referee's file |
| Theorem B | LEAN (over Z); over Z/9Z by an exhaustive computation |
| Theorem C | WRITTEN-AUDITED: p ≥ 5 and R by hand; p = 2, 3 by a finite exact computation with Hensel lemma H′ |
| (−6 −5; −1 −1) and (−6 −6; −5 −6), both outside O, are not X³ + Y³ with \|tr X\|, \|tr Y\| ≤ 500; census of [−9,9]⁴: of 10,000 representatives, 215 have no representation with traces ≤ 150, and 19 searches are incomplete | computational (two independent implementations) |
| Conjecture 3.6: (−6 −5; −1 −1) is not a sum of two cubes in M₂(Z), so the local–global principle fails for two cubes | CONJECTURE |
| Conjecture 3.7: the proportion of matrices of height ≤ N that are sums of two cubes tends to 0 | CONJECTURE (weak evidence) |

## Proof 1: `matrix.cayley_hamilton`
- **Origin:** branch `origin/claude/math-prime-cross-domain-g1vwsn`, tip `26c718ac`, folder
  `cloud_agent_claude/04_matrix_waring_cubes/` (track 04 of the cloud Claude programme, 18–19 October 2026).
  The files are copied byte for byte into [source/](source/), the Lean files into [lean/](lean/).
- **Files:**
  - [T2_three_cubes.md](source/04_matrix_waring_cubes/theorems/T2_three_cubes.md): Theorem A, Lemmas 1–4, Cases A–C;
  - [T1_two_cube_obstruction.md](source/04_matrix_waring_cubes/theorems/T1_two_cube_obstruction.md): Theorem B;
  - [T4_local_two_cubes.md](source/04_matrix_waring_cubes/theorems/T4_local_two_cubes.md): Theorem C, Lemmas D, H,
    H′, S, C;
  - [T3_two_cubes.md](source/04_matrix_waring_cubes/theorems/T3_two_cubes.md): Proposition 3.2, Algorithm 3.3,
    Theorem 3.4, Proposition 3.5, Conjectures 3.6 and 3.7;
  - [scripts/](source/04_matrix_waring_cubes/scripts/), with commands, exit codes and outputs in
    [evidence/](source/04_matrix_waring_cubes/evidence/); the search log
    [LITERATURE.md](source/04_matrix_waring_cubes/literature/LITERATURE.md); the failed approaches in
    [exploration_log.md](source/04_matrix_waring_cubes/notes/exploration_log.md).
- **Steps of Theorem A:**
  1. Two cubes with X³ = X (rank-one idempotents, involutions and their negatives) realise every off-diagonal pair,
     once the diagonal is one of (0,0), (1,1), (−1,−1), (1,0), (−1,0) (Lemma 1).
  2. Shifted cube (Lemma 2, Cayley–Hamilton): if tr W = 3c, then (W − cI)³ = (3c² − det W)·W − c³I.
  3. An SL₂(Z)-conjugate of A gets one diagonal entry with a prescribed residue modulo α = 3c² − det W, because the
     primitive binary form v ↦ det(v, Bv) takes values prime to any modulus (Lemma 3).
  4. A 2-adic valuation argument chooses c ∈ {±1, ±2} and p ∈ {−1, 0, 1} (Lemma 4). Cases: |tr A| ≤ 2, A scalar,
     and the rest.

  Theorem B follows from tr M³ = t³ − 3td (t = tr M, d = det M) and a congruence modulo 9.
- **Lean:**
  - [lean/ThreeCubes.lean](lean/ThreeCubes.lean) (781 lines, `import Mathlib` only) proves `three_cubes`,
    `two_cube_obstruction`, `three_I_not_two_cubes` and `waring_number_three`.
  - Build record: [lean_ThreeCubes.md](source/04_matrix_waring_cubes/evidence/lean_ThreeCubes.md). Exit 0 with
    Lean v4.19.0 and Mathlib `c44e0c8`, the pins of this repository's `lean/` project; the axioms are propext,
    Classical.choice and Quot.sound; there is no `sorry` and no `native_decide`. The branch's referee compiled the
    file again: [lean_compile_author_copy.log](source/reviews/evidence/04_T1T2/lean_compile_author_copy.log).
  - The author's file states only that every matrix is a sum of three cubes. The refinement X³ = X, Y³ = Y,
    tr Z ∈ {0, ±1, ±2} is kernel-checked only in the referee's [lean/Refined04.lean](lean/Refined04.lean), the
    author's file followed by `rev_three_cubes_refined`, which states the whole of Theorem A's first sentence
    ([log](source/reviews/evidence/04_T1T2/lean_compile_refined.log): exit 0, standard axioms). The plain form
    sInf {s | every A is a sum of s cubes} = 3 is `rev_waring_number_eq_three` in the referee's
    [lean/Fidelity04.lean](lean/Fidelity04.lean) ([log](source/reviews/evidence/04_T1T2/lean_compile_fidelity.log)).
  - Nothing was rebuilt in this repository. With its toolchain:
    `cd lean && lake env lean ../research/proofs/theorems/matrix_cubes/lean/ThreeCubes.lean`.
- **Audit:**
  - The branch's reviews, one referee per claim group, from the same session and model family:
    - [04__T1_T2_three_cubes.md](source/reviews/04__T1_T2_three_cubes.md): Theorems A and B PASS, fidelity MATCH.
    - [04__T4_local.md](source/reviews/04__T4_local.md): Theorem C PASS, MATCH.
    - [04__T3_two_cubes.md](source/reviews/04__T3_two_cubes.md): PASS for every item labelled proved; the bounded
      searches were redone without PARI's `qfbsolve` and `quadunit`.
    - [04__T5_L.md](source/reviews/04__T5_L.md): PASS for the erratum below.

    The track's revision after these reviews was never made, so the slips listed in §5 of the first review are
    still in the branch text.
  - Our check, 2026-10-26 (rev1): CONFIRMED. An independent re-implementation of Theorem A, written from the text,
    verified 2 × 16,162 matrices with 0 failures. It excluded 3I modulo 9 by exhausting the 513 cube classes, and
    swept Lemma 4 over 1,220,590 cases:
    [checks/rev1/t04/](../../../../research/shared/consolidation_2026-10-26/review_checks/rev1/t04/).
- **Provenance and novelty:**
  - Gadre–Katre (arXiv math/0702445, 2007), Remark 5.6: "For M2(Z) whether 3 cubes suffice is open". They attribute
    four cubes to an unpublished Katre–Sheth preprint.
  - No later resolution was found, in the branch's searches (LITERATURE.md; about 30 papers screened by its referee)
    or in ours of 2026-10-26 (27 queries:
    [novelty_web_log.md](../../../../research/shared/consolidation_2026-10-26/a_tracks/work/novelty_web_log.md),
    Q1). Not readable: the Katre–Sheth preprint and D. N. Sheth's ICM 2002 short communication "2×2 matrices as
    sums of cubes". The lower bound (3I) was plausibly known to Katre–Sheth.
  - Theorem B is the case A₁ = 0, k = 3, n = 2 over Z/9Z of Choudhary–Kishnani–Singh (CKS), arXiv 2607.07755,
    Prop. 4.1. CKS's proof skips the key step, and their general statement is false (Corrections below). Theorem B
    is the first complete proof of that special case.
  - Theorem C: part 1 for large p is CKS Theorem 1.2, with an explicit constant of about 3.9·10³, and Lemma D is
    CKS Lemma 3.1 in general form. The complete statement, including p = 2, 3 and the exact 3-adic set, was not
    found.
  - The branch's two-square theorem (T5) is Kaneyama (2013), result and method, and is not imported.

## Corrections to the branch's text
- The status of T2 and the track README say that the whole of Theorem 2 is Lean-checked "on the exact statement".
  This holds for the headline only (see Lean above).
- T1's novelty label "NEW (no match found)" should read "first complete proof of a special case claimed in
  arXiv 2607.07755".
- The track README's "no claim of this track has yet been reviewed", and the top README's listing of Theorem 4 as
  under review, are stale: all four reviews were filed with PASS.
- The note on CKS in LITERATURE.md says that the conclusion of Prop. 4.1 is correct for A₁ = 0. That is false in
  general: 9I = I³ + (2I)³ in M₂(Z/27Z), and 27E₁₁ = (3E₁₁)³ (mod 81). It holds only over Z/9Z with n ≤ 2
  (review 04__T5_L).

## Open questions
1. Conjecture 3.6. The missing piece is an explicit bound for |tr X| and |tr Y| in terms of A. For
   (−6 −5; −1 −1) the same 14 candidate tuples survive every bound from 45 to 500. By the diagnostic of T3,
   §3.6(d), four of them give conics without a rational point and ten fail integrality on every unit orbit. The
   branch's track 11 was set up for this question and stopped after 25 minutes, without a claim.
2. Conjecture 3.7.
3. Fourth powers in M₂(Z), and n × n matrices with n ≥ 3 (not started).
