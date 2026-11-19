# Triangular = square-pyramidal: n(n+1)/2 = m(m+1)(2m+1)/6 (8604.a1, rank 2)
Status: WRITTEN (in-session review only). Theorem A (twelve quartic Thue equations) and Theorem B (the integral
points) are proved; the 2-descent class list (H) is proved in [hypothesis_H.md](hypothesis_H.md), without the rank.
Methods: [ghost.unit.multi_unit](../../methods/ghost.unit.multi_unit.md) +
[lifting.binomial.multi](../../methods/lifting.binomial.multi.md) (Proposition 17 only); (H):
[norm.triples](../../methods/norm.triples.md) (Euclidean division in Z[λ], one-unit descent)
Generalizes: none     Generalized by: none

## Statement
Let E: y² = x³ − 9x + 81 (LMFDB 8604.a1) and K = Q(θ), θ³ − 9θ + 81 = 0.

**Theorem A.** For the three class forms below and c ∈ {1, 4, 9, 36}, the integer solutions of F_cl(m,n) = c are
exactly ±(m,n) for the listed pairs (m,n):

| class | F_cl(m,n) | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|---|
| (1,0) | 9m⁴ − 72m³n − 18m²n² + n⁴ | (0,1) | none | (1,0) | none |
| (0,1) | 9m⁴ − 72m³n + 18m²n² + n⁴ | (0,1) | none | (1,0) | (1,3) |
| (1,1) | 36m⁴ − 72m³n + n⁴ | (0,1), (2,1) | (1,4) | none | (1,0) |

It is proved through the twelve equations G(m,n) = c of the monic *field forms* (Theorem A′ of
[proof.md](proof.md)). For class (1,0) the field form is G = 4m⁴ − 5m³n − 3m²n² − 2mn³ + n⁴ with
F_cl(m,n) = G(2m, m+n), because (1+α)/2 is integral for a root α of F_cl(1,t). For the other two classes G = F_cl.

**Theorem H** ([hypothesis_H.md](hypothesis_H.md)). Every integral point P of E has x(P) − θ ∈ μK*² for some
μ ∈ {1, 3−θ, −3−θ, −θ}. These are the classes of the 2-descent in
[weierstrass_2descent.out.txt](../../targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt), section CURVE
8604.a1; this is hypothesis (H), which rank E(Q) ≤ 2 would also imply.

**Theorem B.**
- the integral points of E are the 24 points with x ∈ {−5, −3, 0, 3, 7, 9, 24, 33, 39, 513, 1099, 5112} and both signs
  of y;
- the numbers that are both triangular and square-pyramidal, n(n+1)/2 = m(m+1)(2m+1)/6 > 0, are exactly 1, 55, 91
  and 208335 (m = 1, 5, 6, 85).

## Proof 1: strategy R1 (`ghost.unit.multi_unit` + `lifting.binomial.multi`)
- **Files:**
  - [proof.md](proof.md): the written proof, §§0–7;
  - [scripts/](scripts/): each script with its output (`*.out.txt`). Shared exact arithmetic is in
    [core.py](scripts/core.py), and the exclusion helpers are in [lifting2.py](scripts/lifting2.py), a verbatim copy
    of the template's;
  - [literature_check.md](literature_check.md): prior proofs and provenance;
  - [review.md](review.md): the in-session verification and review record.
- **Steps:**
  1. *Reduction* (prior work: [weierstrass_2descent.out.txt](../../targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt);
     re-derived in [03_reconstruction](scripts/03_reconstruction.out.txt)).
     - The 2-descent over K maps each non-trivial class to coprime (m,n) with F_cl(m,n) = g², g ∈ {1, 2, 3, 6},
       through the parametrization of a conic (Lemma 1.3′ of [proof.md](proof.md)).
     - The trivial class is solved directly.
     - The class list (H) is Theorem H of [hypothesis_H.md](hypothesis_H.md); integral closedness of Z[λ] (there,
       (3.6)) replaces the facts N1 and N2 of the descent script.
  2. *Thue fields* ([00_fields](scripts/00_fields.out.txt)).
     - G(m,n) = N(n − mβ) in the maximal order O of Q(β); the fields are S4 of signature (2,1).
     - A solution is an element of norm c whose two marker coordinates vanish.
  3. *Coverage: the two-unit descent* ([01_coverage](scripts/01_coverage.out.txt)).
     - A least-E element of a unit orbit satisfies eight linear inequalities.
     - The vertices of their polygon bound E ≤ 5.31√c, 101.6√c and 200.8√c in the three fields.
     - An exact enumeration then gives 12 orbit seeds, one per class and c.
     - By-product: O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩ in each field, without the unit theorem.
  4. *Exclusion* ([02_exclusion](scripts/02_exclusion.out.txt)).
     - The kernel lattices of the units are taken modulo q ∈ {4, 7, 8, 13, 16, 23, 43, 49}.
     - Each known zero is isolated by the Jacobian criterion (Proposition 17).
     - Every other candidate class is removed, either by congruences or as a sibling: a class inside the Λ_q-class
       of a known zero at which Proposition 17 holds for the same q.
     - The zero sets are exactly the known solutions.
  5. *Reconstruction* ([03_reconstruction](scripts/03_reconstruction.out.txt)): back to the 24 points and to
     T(n) = P(m).
- **What is elementary.**
  - Theorem A uses exact integer and rational computation, with real numbers only in rational interval enclosures
    (sign changes and Sturm's theorem). It also uses the convexity lemma 3.2, the finite-polynomial lifting lemma and
    congruences.
  - Theorem B adds the 2-descent and (H), proved in [hypothesis_H.md](hypothesis_H.md): Z[λ] (λ = θ/3) is
    norm-Euclidean, certified by 116 exact boxes, which gives the divisibility theory that class number 1 provided
    (Bézout, prime elements, square allocation); its units are ±εᵏ, by a one-unit descent; and a congruence modulo 4
    removes the unit's coset. Integral closedness of Z[λ] replaces the facts N1 and N2. No rank, ideals or class
    numbers are used.
- **Lean:** none.
- **Audit:** none outside the session. [review.md](review.md) records the in-session checks; round 3 covers
  hypothesis_H.md.
- **Provenance and novelty** ([literature_check.md](literature_check.md)):
  - The theorem is classical: Avanesov (Volzh. Mat. Sb. 8, 1971), Finkelstein–London (J. Number Theory 4, 1972) and
    Uchiyama (Tsukuba J. Math. 8, 1984, 131–157). The first two were not read.
  - Uchiyama's paper was read in full. His cubic ring is our O_K; his twelve sub-case equations are our twelve
    equations F_cl = c; his three quartic fields and his unit groups are ours; his exclusion is a two-index argument
    written as finite congruences modulo powers of 2.
  - His printed Thue step is incomplete: wrong integral bases in two fields, printed "units" that are not units in
    three propositions, and a congruence step that is false in most tested cases and whose conclusion fails near the
    known solutions. His Prop. 8 cannot be completed modulo powers of 2 at all: the markers of its orbit have a 2-adic
    zero that is not an integer zero, and odd moduli are needed.
  - The theorem is right: Theorem A here with Uchiyama's §1 (given h_K = 1 and the unit of K) proves it. This card
    proves it without that input: hypothesis_H.md replaces his class number (stated without a source) and his
    tabulated unit.
  - **No novelty is claimed.** This is a reconstruction of Uchiyama's proof, complete and checked, with a reduction
    that uses no ideals or class numbers.

## Open questions
- *Done:* an elementary proof of (H) for integral points ([hypothesis_H.md](hypothesis_H.md)). R1 is now complete
  and elementary on this curve.
- The exclusion sieve's termination is checked, not proved a priori.
