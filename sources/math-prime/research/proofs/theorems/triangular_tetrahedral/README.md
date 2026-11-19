# Triangular = tetrahedral: C(n,2) = C(m,3) (1611.a1, rank 2)
Status: WRITTEN (in-session review only). Theorem A (twelve quartic Thue equations) and Theorem B (the integral
points) are both proved. The 2-descent hypothesis (H), first taken from rank E(Q) ≤ 2, is proved in
[hypothesis_H.md](hypothesis_H.md) (2026-10-15).
Methods: [ghost.unit.multi_unit](../../methods/ghost.unit.multi_unit.md) +
[lifting.binomial.multi](../../methods/lifting.binomial.multi.md) (Proposition 17 only); for (H):
[norm.triples.multiplier](../../methods/norm.triples.multiplier.md)
Generalizes: none     Generalized by: none

## Statement
**Theorem A.** For the three forms below and c ∈ {1, 4, 9, 36}, the integer solutions of F(m,n) = c are exactly
±(m,n) for the listed pairs (m,n):

| class | F(m,n) | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|---|
| (1,0) | 36m⁴ − 6mn³ + n⁴ | (0,1) | (1,2) | none | (1,0), (1,6) |
| (0,1) | −5m⁴ + 14m³n + 24m²n² + 2mn³ + n⁴ | (0,1) | (−1,1) | (4,1) | (1,1) |
| (1,1) | −8m⁴ + 4m³n + 24m²n² + 4mn³ + n⁴ | (0,1) | none | (−1,1), (2,1) | none |

**Theorem B.** Let E: Y² = X³ − 144X + 1296 and K = Q(θ), θ³ − 144θ + 1296 = 0. Then (H) holds: every integral
point P of E has X(P) − θ ∈ μK*² for some μ ∈ {1, 9−θ, −15−θ, −12−θ} (Theorem H of
[hypothesis_H.md](hypothesis_H.md)). These are the classes of the 2-descent in
[weierstrass_2descent.out.txt](../../targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt), section CURVE
1611.a1. Hence:
- the integral points of E are the 28 points with X ∈ {−15, −12, −8, 0, 4, 9, 12, 24, 40, 48, 108, 252, 420, 4953} and
  both signs of Y;
- the LMFDB model y² + y = x³ − 9x + 20 (1611.a1) has 22 integral points;
- the numbers that are both triangular and tetrahedral, C(n,2) = C(m,3) > 0, are exactly 1, 10, 120, 1540 and 7140.

## Proof 1: strategy R1 (`ghost.unit.multi_unit` + `lifting.binomial.multi`)
- **Files:**
  - [proof.md](proof.md): the written proof, §§0–7;
  - [scripts/](scripts/): each script with its output (`*.out.txt`). Shared exact arithmetic is in
    [core.py](scripts/core.py), and the exclusion helpers are in [lifting2.py](scripts/lifting2.py);
  - [literature_check.md](literature_check.md): Avanesov's paper read in full, other proofs, provenance;
  - [review.md](review.md): the in-session verification and review record.
- **Steps:**
  1. *Reduction* (prior work: [weierstrass_2descent.out.txt](../../targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt)).
     - The 2-descent over K = Q(θ), θ³ − 144θ + 1296 = 0, maps each non-trivial class to coprime (m,n) with
       F(m,n) = g², g ∈ {1, 2, 3, 6}. The conic parametrization this uses is proved in proof.md Lemma 1.3′ (written
       2026-10-15; before that it was only checked by script).
     - The trivial class is solved directly.
     - (H) is Theorem H of [hypothesis_H.md](hypothesis_H.md). It is proved by:
       - a division with multipliers {1, 2} in Z[λ], λ = θ/6, certified by 398 exact boxes;
       - Bézout, prime elements and square allocation;
       - the units ±εᵏ, by a one-unit descent;
       - a congruence modulo 4.

       Z[λ] is not norm-Euclidean: the point ξ* fixed by ε needs norm 121/109. The facts N1 and N2 of the descent
       script are replaced by the integral closedness of Z[λ] (hypothesis_H.md (3.6)).
  2. *Thue fields* ([00_fields](scripts/00_fields.out.txt)).
     - F(m,n) = N(n − mα) in the maximal order O of Q(α), α a root of F(1,t); the fields are S4 of signature (2,1).
     - A solution is an element of norm c whose two marker coordinates vanish.
  3. *Coverage: the two-unit descent* ([01_coverage](scripts/01_coverage.out.txt)).
     - A least-E element of a unit orbit satisfies eight linear inequalities.
     - The vertices of their polygon bound E ≤ 34√c, 27√c and 13√c in the three fields.
     - An exact enumeration then gives 16 orbit seeds.
     - By-product: O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩ in each field, without the unit theorem.
  4. *Exclusion* ([02_exclusion](scripts/02_exclusion.out.txt)).
     - The kernel lattices of the units are taken modulo q ∈ {3, 4, 7, 8, 9, 13, 16, 25, 27, 31, 41}.
     - Each known zero is isolated by the Jacobian criterion (Proposition 17).
     - Every other candidate class is removed, either by congruences or by Proposition 17's uniqueness in the
       Λ_q-class of a known zero (siblings).
     - The zero sets are exactly the known solutions.
  5. *Reconstruction* ([03_reconstruction](scripts/03_reconstruction.out.txt)): back to the 28 points and to
     C(n,2) = C(m,3).
- **What is elementary.**
  - Theorem A uses exact integer and rational computation, with real numbers only in rational interval enclosures
    (sign changes and Sturm's theorem). It also uses the convexity lemma 3.2, the finite-polynomial lifting lemma
    and congruences.
  - Theorem B adds the 2-descent and [hypothesis_H.md](hypothesis_H.md), which use:
    - exact arithmetic in Z[λ];
    - a box table checked in exact rationals;
    - interval enclosures of one real root, for the units;
    - congruences.

    No rank, Mordell–Weil theorem, ideal theory or class number is used.
  - A second proof of unique factorization, by Dedekind's box principle, is included as a cross-check. It is not part
    of the proof: its lattices P(p, g) are prime ideals in modern language, so only Sections 1–5 of hypothesis_H.md
    meet the no-ideal-theory standard.
- **Lean:** none.
- **Audit:** none outside the session. [review.md](review.md) records the in-session checks:
  - independent re-implementations of the coverage and of the exclusion;
  - a brute-force box;
  - planted negative tests of Proposition 17;
  - adversarial audits and a referee round on proof.md;
  - for hypothesis_H.md, the build checks and the verification and referee rounds recorded in review.md
    (rounds 3 and 4).
- **Provenance and novelty** ([literature_check.md](literature_check.md)):
  - The theorem is Avanesov's (Acta Arith. 12, 1966/67). It was re-solved by Stroeker–de Weger (1999, elliptic
    logarithms) and by Gallegos-Ruiz–Katsipis–Tengely–Ulas (2020, Magma).
  - The mechanism is Avanesov's: cubic-field descent, then the same three (2,1) quartic fields, then Skolem's method
    at p = 3 with Proposition 17's determinant.
  - Its congruence-only form has a precedent on the sister curve 8604.a1 (Uchiyama, Tsukuba J. Math. 8, 1984).
  - The coverage lemma is an effective form of Koecher's reduction theory.
  - **No novelty is claimed.**
  - Avanesov's remark that all three fields have discriminant −2⁴·3²·179 is false for one of them, which has
    −2²·3²·179.
  - *For (H)* ([literature_check.md](literature_check.md) §5):
    - Avanesov took class number 1 from the Delone–Faddeev table.
    - Cavallar–Lemmermeyer (1998) list the field of discriminant −716 as "N 121/109": class number 1, not
      norm-Euclidean. Their value rests on a floating-point cover, with exact checks at the exceptional points.
      hypothesis_H.md gives a written proof of the lower bound at ξ*; the upper bound is not proved there.
    - The division with multipliers is the Dedekind–Hasse criterion (Hasse 1928), and §6 uses Dedekind's box
      principle (1894).
    - No earlier elementary treatment of (H) for this curve was found. No novelty is claimed.

## Open questions
- *Done 2026-10-15:* (H) for integral points, in [hypothesis_H.md](hypothesis_H.md). The same method (multiplier
  division, one-unit descent, square allocation, a congruence) should give (H) for 8604.a1 and 10404.f2. It did, the
  same day (reviewed 2026-10-26): [triangular_square_pyramidal](../triangular_square_pyramidal/hypothesis_H.md), whose
  field (discriminant −239) is norm-Euclidean, so a plain Euclidean division suffices, and
  [mordell_k_minus17](../mordell_k_minus17/hypothesis_H.md), with multipliers {1, 2}. Both cards are WRITTEN.
- The same scripts, with new data, treat 8604.a1 and 10404.f2: see their cards.
- The exclusion sieve's termination is checked, not proved a priori.
