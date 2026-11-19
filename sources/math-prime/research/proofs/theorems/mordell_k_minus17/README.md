# Mordell's equation y² = x³ + 17 (x² − 17 = y³; 10404.f2, rank 2)
Status: WRITTEN (in-session review only). Theorem A (eight quartic Thue equations) and Theorem B (the integral
points) are proved; the 2-descent class list (H) is proved in [hypothesis_H.md](hypothesis_H.md), without the rank and
without Cassels's ideal-theoretic proof.
Methods: [ghost.unit.multi_unit](../../methods/ghost.unit.multi_unit.md) +
[lifting.binomial.multi](../../methods/lifting.binomial.multi.md) (Proposition 17 only); (H):
[norm.triples.multiplier](../../methods/norm.triples.multiplier.md) (division with multipliers {1, 2}, one-unit
descent)
Generalizes: none     Generalized by: none

## Statement
Let E: y² = x³ + 17 (LMFDB 10404.f2) and K = Q(θ), θ³ + 17 = 0.

**Theorem A.** For the three class forms below and c among the listed right sides, the integer solutions of
F_cl(m,n) = c are exactly ±(m,n) for the listed pairs (m,n):

| class | F_cl(m,n) | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|---|
| (1,0) | −12m⁴ − 24m³n + 12m²n² + n⁴ | (0,1) | (1,2) | — | — |
| (0,1) | 16m⁴ − 8m³n − 8mn³ + n⁴ | (0,1), (1,1) | none | — | — |
| (1,1) | −3m⁴ + 20m³n − 12m²n² + 4n⁴ | none | (0,1), (6,1) | (1,1) | none |

It is proved through the eight equations G(m,n) = c of the monic *field forms* (Theorem A′ of [proof.md](proof.md)):
- class (1,0): G = F_cl;
- class (0,1): G = m⁴ − m³n − 4mn³ + n⁴ with F_cl(m,n) = G(2m, n), because α/2 is integral for a root α of
  F_cl(1,t);
- class (1,1): G = −12m⁴ + 40m³n − 12m²n² + n⁴ with 4F_cl(m,n) = G(m, 2n) and F_cl(2m, n) = 4G(m,n), because F_cl
  is not monic.

**Theorem H** ([hypothesis_H.md](hypothesis_H.md)). Every integral point P of E has x(P) − θ ∈ μK*² for some
μ ∈ {1, −2−θ, 4−θ, −1−θ}. These are the classes of the 2-descent in
[weierstrass_2descent.out.txt](../../targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt), section CURVE
10404.f2; this is hypothesis (H), which rank E(Q) ≤ 2 would also imply.

**Theorem B.** The integral points of E are the 16 points with x ∈ {−2, −1, 2, 4, 8, 43, 52, 5234} and both signs of y.

(H) is also a published theorem: Cassels (Acta Math. 82, 1950, §22) proves the four-class list for this curve, with
ideals, the class number and a unit of Q(∛17), and no class field theory. hypothesis_H.md replaces that input. The
theorem itself is classical (Nagell 1930).

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
     - The 2-descent over K maps each non-trivial class to coprime (m,n) with F_cl(m,n) = g², g ∈ {1, 2} for classes
       (1,0) and (0,1) and g ∈ {1, 2, 3, 6} for class (1,1), through the parametrization of a conic (Lemma 1.3′ of
       [proof.md](proof.md)).
     - The trivial class is solved directly.
     - The class list (H) is Theorem H of [hypothesis_H.md](hypothesis_H.md); integral closedness of O_K (there,
       (3.6)) replaces the facts N1 and N2 of the descent script.
  2. *Thue fields* ([00_fields](scripts/00_fields.out.txt)).
     - G(m,n) = N(n − mβ) in the maximal order O of Q(β); the fields are S4 of signature (2,1).
     - A solution is an element of norm c whose two marker coordinates vanish.
  3. *Coverage: the two-unit descent* ([01_coverage](scripts/01_coverage.out.txt)).
     - A least-E element of a unit orbit satisfies eight linear inequalities.
     - The vertices of their polygon bound E ≤ 369.1√c, 7.12√c and 376.3√c in the three fields.
     - An exact enumeration then gives 17 orbit seeds. Class (0,1) has no element of norm 4; class (1,1) has six
       orbits each for c = 9 and c = 36.
     - By-product: O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩ in each field, without the unit theorem.
  4. *Exclusion* ([02_exclusion](scripts/02_exclusion.out.txt)).
     - The kernel lattices of the units are taken modulo q ∈ {3, 4, 5, 7, 8, 9, 16, 19, 27}.
     - Each known zero is isolated by the Jacobian criterion (Proposition 17).
     - Every other candidate class is removed, either by congruences or as a sibling: a class inside the Λ_q-class
       of a known zero at which Proposition 17 holds for the same q.
     - The zero sets are exactly the known solutions.
  5. *Reconstruction* ([03_reconstruction](scripts/03_reconstruction.out.txt)): back to the 16 points.
- **What is elementary.**
  - Theorem A uses exact integer and rational computation, with real numbers only in rational interval enclosures
    (sign changes and Sturm's theorem). It also uses the convexity lemma 3.2, the finite-polynomial lifting lemma and
    congruences.
  - Theorem B adds the 2-descent and (H), proved in [hypothesis_H.md](hypothesis_H.md): the ring Z + Zθ + Zω is not
    norm-Euclidean, but a division with multipliers {1, 2} (Dedekind–Hasse), certified by 795 exact boxes, gives the
    divisibility theory that class number 1 provided (Bézout, prime elements, square allocation); its units are ±ηᵏ,
    η = 18 + 7θ, by a one-unit descent; and the norm condition alone leaves the four classes. Integral closedness replaces the facts N1 and N2. No rank, ideals or class numbers are used
    (rank 2 is also proved by Cassels (1950) and, unconditionally, by Gebel–Pethő–Zimmer (1998)).
- **Lean:** none.
- **Audit:** none outside the session. [review.md](review.md) records the in-session checks; round 3 covers
  hypothesis_H.md.
- **Provenance and novelty** ([literature_check.md](literature_check.md)):
  - The theorem is classical: Nagell (1930), then Hemer (1952, 1954), London–Finkelstein (1973), Gebel–Pethő–Zimmer
    (1998, elliptic logarithms) and Bennett–Ghadermarzi (2015, cubic Thue equations with Baker's method).
  - The algebraic proofs whose method is known (Hemer by its review; Nagell by inference) take the cubic route.
    Under Mordell's correspondence the 16 points fall into four cubic Thue equations, in four cubic fields with one
    unit each (computed here). R1 takes the quartic route, through Q(∛17): eight quartic Thue equations in three S4
    fields with two units each. The routes share no field and no equation.
  - Related: [mordell_negative_k](../mordell_negative_k/README.md) (LEAN) proves finiteness of x² + k = y³ for every
    k < 0, so for this equation, by the one-unit cubic route; it does not give the list.
  - No proof by the quartic route, or by any two-exponent Skolem-type (unit) argument, was found for this theorem, for
    x² − 17 = y³, or for a family containing it. **No novelty is claimed yet:** Nagell 1930, Hemer 1952/1954,
    London–Finkelstein 1973 and Ljunggren 1963 were not read, and whether the outline shared with the cubic route makes
    R1 a rewording is the owner's decision.
  - Not new: the quartic route itself (named by Bennett–Ghadermarzi; implemented in KANT with Baker-type solving),
    the two-index mechanism (Skolem 1933–35, Ljunggren 1942; on the sister curves, Avanesov and Uchiyama), the
    one-exponent congruence lemma (Delone; Hemer), the descent's class list (Cassels 1950) and the coverage lemma
    (Koecher's reduction theory).

## Open questions
- *Done:* an elementary proof of (H) for integral points ([hypothesis_H.md](hypothesis_H.md)). R1 is now complete and
  elementary on this curve.
- A complete elementary proof by the cubic route (four one-unit cubic Thue equations), as an independent second
  proof.
- The exclusion sieve's termination is checked, not proved a priori.
