# ghost.unit.multi_unit: unit ghosts with two or more independent units
Parent: [ghost.unit](ghost.unit.md)     Children: none yet

**Status: WRITTEN (2026-10-15)**, as part of [triangular_tetrahedral](../theorems/triangular_tetrahedral/README.md),
its first use. It was refereed in-session only ([review.md](../theorems/triangular_tetrahedral/review.md)).

## What it does
It extends `ghost.unit` from one unit to several, in two halves.
- **Coverage.** Every element of given norm c is a unit times one of finitely many seeds. The seeds are found by an
  explicit linear certificate and an exact enumeration. No unit theorem, fundamental unit, ideal or class number is
  used.
- **Exclusion.** On each orbit, the marker coordinates vanish only at the known solutions. This uses:
  - kernel lattices of the units modulo prime powers;
  - the Jacobian criterion of [lifting.binomial.multi](lifting.binomial.multi.md) (Proposition 17);
  - a congruence sieve.

## Coverage lemma (reduced representatives)
**Setting.** Let L be a number field with places v and O an order of L. For x ∈ O ⊗ R put s_v(x) = |σ_v(x)|² and
E(x) = Σ_v s_v(x). Let G be a group of units **of norm +1**, U ⊂ G finite, and f_v(u) = s_v(u).
1. **Reduction.** For N(z) = c, choose w of least E in the finite set {gz : g ∈ G, E(gz) ≤ E(z)}. Then N(w) = c, and
   E(uw) − E(w) = Σ_v s_v(w)(f_v(u) − 1) ≥ 0 for every u ∈ U. So s(w) lies in the cone cut out by these linear
   inequalities.
2. **Certificate.** Replace the f_v(u) by rational upper bounds, and normalize by Σ_v r_v = 1. If every vertex of
   the resulting polytope has all coordinates positive, then the product Π r_v^{d_v} (d_v = 1 real, 2 complex) is at
   least its least value b_min at a vertex, by concavity of its logarithm. Hence E(w)^{deg L} ≤ c²/b_min.
3. **Enumeration.** The seeds are the elements of norm c with E ≤ T, found by an exact Fincke–Pohst enumeration with a
   rational lower Gram matrix, and merged into orbits only through exactly verified unit quotients.
   - If G contains units of norm −1, enumerate norm ±c instead.
4. **A purely linear variant.** s_v ≥ κ_v E, with κ_v the least v-coordinate over the vertices. The box is larger.
5. **Independence.** The certificate forces the log size vectors of U to span the hyperplane Σ d_v y_v = 0. So the
   units are independent, and G has full rank.
6. **The unit group.** With c = ±1, the same enumeration decides whether the given units are fundamental.

**Where the proof is.** [proof.md](../theorems/triangular_tetrahedral/proof.md) §3 proves the case used: degree 4,
signature (2,1), weights 1 (Lemma 3.2, Remarks (i)–(iii)). The general case is the same argument with the weights
d_v. Code: [01_coverage.py](../theorems/triangular_tetrahedral/scripts/01_coverage.py).

## Exclusion on an orbit
For a seed w and a prime power q ≠ 2, put Λ_q = {(a,b) : ε₁^aε₂^b ≡ 1 mod q}.
1. Each class mod Λ_q is ≡ its seed mod q. So the candidates are the classes whose residue has vanishing markers.
2. A class containing a known zero z is certified by Proposition 17 at z. Its coset system is multiplication by a
   basis of units ≡ 1 mod q.
3. **Sieve.** Classes are refined modulo intersections of the Λ_q. A class that contains no known zero, but lies in the
   Λ_q-class of a known zero at which Proposition 17 holds **for that same q**, is empty.
4. The orbit is done when every candidate class is certified or empty.

Two points about signs and decomposition:
- **Phantom candidates.** When −1 ≡ ε^{(a,b)} (mod q), the negative of a solution gives a phantom candidate. A
  modulus at which −1 is not a residue of the unit group removes it, or the sibling rule at a second modulus does.
- **Use Λ_q, not h₁Z×h₂Z.** The cosets of h₁Z×h₂Z overcount by [Λ_q : h₁Z×h₂Z]. When p ∤ h₁h₂, those cosets share
  one p-adic residue disc.

Code: [lifting2.py](../theorems/triangular_tetrahedral/scripts/lifting2.py),
[02_exclusion.py](../theorems/triangular_tetrahedral/scripts/02_exclusion.py).

## Requirements
- **Coverage.** Units whose log size vectors span the trace-zero hyperplane, which needs the full unit rank. Any
  signature works.
- **Exclusion.**
  - As many markers as units. For a quartic Thue equation this means signature (2,1): 2 markers against 2 units.
  - For each known zero, some modulus at which it passes Proposition 17.
  - A sieve that terminates. Termination is not guaranteed a priori.

## Limits
- The exclusion sieve has no a priori termination guarantee. A profinite zero of the markers that is not an integer
  zero would block every sieve; a p-adic one blocks only the moduli p^k. One 2-adic zero that is not an integer zero
  was met, in 8604.a1 (orbit (0,1), c = 36; Hensel certificate in
  [check_2adic_zero.py](../theorems/triangular_square_pyramidal/scripts/check_2adic_zero.py)), and odd moduli passed
  it.
- More units than markers is Skolem's failure case, for example totally real quartic fields. The coverage still
  works there; the exclusion does not.
- It solves the Thue equations it is given. Reaching those equations from integral points is a separate step: for
  1611.a1, 8604.a1 and y² = x³ + 17 it is the 2-descent, which needs hypothesis (H), implied by rank ≤ 2. For
  y² = x³ + 17, (H) is also Cassels's theorem (1950, §22), which uses ideals. For all three, (H) is proved
  elementarily: [1611.a1](../theorems/triangular_tetrahedral/hypothesis_H.md),
  [8604.a1](../theorems/triangular_square_pyramidal/hypothesis_H.md),
  [y² = x³ + 17](../theorems/mordell_k_minus17/hypothesis_H.md).

## Theorems using it
- [triangular_tetrahedral](../theorems/triangular_tetrahedral/README.md) (Theorem A: twelve quartic Thue equations of
  signature (2,1); 1611.a1).
- [triangular_square_pyramidal](../theorems/triangular_square_pyramidal/README.md) (Theorem A: twelve quartic Thue
  equations of signature (2,1); 8604.a1).
- [mordell_k_minus17](../theorems/mordell_k_minus17/README.md) (Theorem A: eight quartic Thue equations of signature
  (2,1); 10404.f2, y² = x³ + 17).

## Provenance
- **Exclusion.** It is Skolem's method (1933–35) in finite-polynomial form.
  - Avanesov used the same mechanism on the same fields in 1966/67, at p = 3, with p-adic series.
  - Uchiyama (Tsukuba J. Math. 8, 1984) wrote a two-index argument of this kind entirely as congruences mod 2^{h+c},
    for the sister curve 8604.a1. Read in full on 2026-10-15, his printed argument has gaps near the known solutions
    ([triangular_square_pyramidal literature check](../theorems/triangular_square_pyramidal/literature_check.md)
    §1.7). Proposition 17, with a kernel-lattice sieve that includes odd moduli, replaces his exclusion; for his
    Prop. 8 no argument modulo powers of 2 suffices.
- **Coverage.** It is an effective, certificate form of Koecher's reduction theory of unary forms: the minimizers of
  the trace form over a unit orbit form a polyhedral cone cut out by finitely many units.
  - Koecher 1960 was not accessed; the attribution follows Porter–Bali–Leibak, arXiv:2302.00505.
  - Bounds for a reduced associate are classical: Levesque–Waldschmidt, arXiv:1505.06653, Lemma 6; Shorey–Tijdeman,
    Lemma A.15; Fieker's thesis, 1997.
  - The weights-1 vertex-product certificate was not found in the sources read, and nothing is claimed for it.
- **By-product.** For 1611.a1, 8604.a1 and 10404.f2, O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩ in each of the nine fields, proved
  without the unit theorem.
- Sources: [literature check](../theorems/triangular_tetrahedral/literature_check.md).
- No novelty is claimed.
