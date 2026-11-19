# norm.triples: Euclidean arithmetic of integer triples
Parent: norm     Sibling: [norm.pairs](norm.pairs.md)     Children: [norm.triples.multiplier](norm.triples.multiplier.md)

## What it does
Triples (a,b,c) with an explicit multiplication rule stand for the elements of a cubic ring in a fixed basis, usually
a+bθ+cθ² (for Q(∛17) the basis is 1, θ, ω with ω = (1+θ+θ²)/3). It was first done for θ³=4θ−2, the ring of
Mordell's equation y(y+1)=x(x+1)(x+2); then for λ³ = λ − 3 and, through the child, for λ³ = 4λ − 6 and the ring of
integers of Q(∛17) (hypothesis (H) of three R1 cards). Every property used is derived directly:
- multiplication by reduction modulo the cubic;
- a norm without conjugates, the determinant of the multiplication matrix, and an adjugate triple α* with
  αα*=N(α);
- a proof that the norm vanishes only at 0 (here a 2-adic descent, since t³−4t+2 is Eisenstein at 2);
- Euclidean division, from a finite table of rational boxes that cover a fundamental domain, each with a
  translation and an exact Taylor bound below 1;
- gcd and Bézout identities, prime elements of prime norm, exact powers, and coprime square allocation;
- the units, by descent on a positive definite height, and from them the units modulo squares. For a totally real
  ring the height is T(u)=Tr(u²), a positive integer; for a complex cubic, with one unit, it is T₂ = σ₁² + 2Q, a
  real-valued form handled through rational enclosures of the real root.

No class-number theorem, unit theorem or table of field invariants is assumed. The integers of the field are not
needed either, because all arithmetic stays among integer triples.

Proof 2 of [mordell_1963](../theorems/mordell_1963/README.md) runs the same arithmetic on all rational points of the
curve, not only on integral ones. Mordell's factorization X − θ = (unit)·θ·(square) is extended to reduced rational
points, the unit class modulo squares is read off from the real component of the point, and the square root then
gives the half of a point (halving). This is the descent step of [ghost.rational](ghost.rational.md) for a curve with
no isogeny. It is the classical complete 2-descent through the cubic field, whose halving step is Mordell's
(Diophantine Equations, Ch. 16, Thm 6; also Silverman, AEC Ex. 10.9(a)), here with the field facts proved rather than
taken from tables.

## Ingredients (with the file that proves each)
All in [theorems/mordell_1963/](../theorems/mordell_1963/README.md):

| Ingredient | Proof | Check |
|---|---|---|
| Multiplication, determinant norm, adjugate triple; the norm vanishes only at 0 | [proof.md](../theorems/mordell_1963/proof.md) §1 | [01_triples.py](../theorems/mordell_1963/scripts/01_triples.py) |
| Euclidean division: 24 boxes cover [−½,½]²×[0,½], every bound below 1 (largest 2023/2048) | §2 | [02_euclid_cover.py](../theorems/mordell_1963/scripts/02_euclid_cover.py) |
| Gcd, Bézout, Euclid's lemma, primes of prime norm, exact powers, square allocation | §3 | (logic only) |
| Units ±εʳηˢ with ε=θ−1, η=2θ−1: three multipliers lower T above 52, and 20 units have T≤52 | §4 | [03_units.py](../theorems/mordell_1963/scripts/03_units.py) |
| The class of a rational point, read off from its real component, and the halving lemma | [ghost/proof.md](../theorems/mordell_1963/ghost/proof.md) (3.1)–(3.8) | [02_class.py](../theorems/mordell_1963/ghost/scripts/02_class.py), [03_halving.py](../theorems/mordell_1963/ghost/scripts/03_halving.py) |

## Requirements: when it applies
- **Division.** Every point of the closed fundamental cube must have an integer translate with |N|<1. A finite
  table of boxes then exists by compactness, and a search finds it (here
  [cover_search.py](../theorems/mordell_1963/scripts/cover_search.py), a k-d split to depth 7). For θ³=4θ−2 a grid
  search suggests that the worst point is θ²/2, where the best translate has |N|=½ exactly (Remark 2.3 of the
  proof proves the lower bound ½ there). The gap between ½ and 1 is what keeps the table small.
- **Units.** The descent as written is for a totally real cubic, with two independent units. It needs rigorous
  intervals for the three real roots and, for each real place k, a unit that is smaller than 1 at the two
  other places. A complex cubic needs one unit and a one-dimensional version of the same descent. This is done
  for λ³ = 4λ − 6 in the child [norm.triples.multiplier](norm.triples.multiplier.md)
  ([hypothesis_H.md](../theorems/triangular_tetrahedral/hypothesis_H.md) §4) and for the ring of integers of Q(∛17),
  basis (1, θ, ω) ([mordell_k_minus17/hypothesis_H.md](../theorems/mordell_k_minus17/hypothesis_H.md) §4); and for
  λ³ = λ − 3, with a plain Euclidean division, in
  [triangular_square_pyramidal/hypothesis_H.md](../theorems/triangular_square_pyramidal/hypothesis_H.md) §4.
- **Allocation.** Coprime square allocation is the cube allocation of [norm.pairs](norm.pairs.md) with exponent 2. A common divisor of the two
  factors that is a prime element is separated by exact powers, as with 4θ−3 in Mordell's equation.

## Limits: where it stops, with examples
- It needs a norm-Euclidean ring. Without one, division fails and the allocation lemma has no proof. The
  analogous failure for pairs is c=−23 (see [norm.pairs](norm.pairs.md)).
  - Failing cubic examples: λ³ = 4λ − 6 (discriminant −716), where the point ξ* fixed by the unit needs norm
    121/109 ([hypothesis_H.md](../theorems/triangular_tetrahedral/hypothesis_H.md), Remark 2.5); and the ring of
    integers of Q(∛17) (discriminant −867), where a point fixed by −η needs 1115/1028
    ([mordell_k_minus17/hypothesis_H.md](../theorems/mordell_k_minus17/hypothesis_H.md), Remark 2.5).
  - A complex cubic that works: λ³ = λ − 3 (discriminant −239), with 116 boxes; the point fixed by the unit needs
    8/9 ([triangular_square_pyramidal/hypothesis_H.md](../theorems/triangular_square_pyramidal/hypothesis_H.md),
    Remark 2.5).
  - For λ³ = 4λ − 6 and for Q(∛17) the child [norm.triples.multiplier](norm.triples.multiplier.md) replaces the division by a
    Dedekind–Hasse division with multipliers {1, 2}. Another ring needs its own certificate, and integer
    multipliers need not suffice even when the class number is 1.
- The box table is a certificate found by computer search and checked by exact rational arithmetic. For θ³=4θ−2
  every bound is at most 2023/2048 (all below 1), and 24 boxes suffice; the worst point needs ½. One more
  measurement: for λ³ = λ − 3, whose worst point needs 8/9, the table has 116 boxes (largest Taylor bound
  130993/131072). Two rings do not show how the length grows as the worst point approaches 1.
- The unit descent needs a threshold, the largest (d_k−c_k)/(1−c_k) over the three multipliers (here below 53),
  and an enumeration of the units below it. The threshold grows with the largest real values of the
  multipliers.

## Theorems using it
- [mordell_1963](../theorems/mordell_1963/README.md): y(y+1)=x(x+1)(x+2). Proof 1 uses it on integral points (Mordell's
  case analysis); proof 2 on rational points (halving).
- [triangular_tetrahedral](../theorems/triangular_tetrahedral/README.md), hypothesis (H), through the child
  [norm.triples.multiplier](norm.triples.multiplier.md): §3 there is this card's divisibility theory, from Bézout on.
- [triangular_square_pyramidal](../theorems/triangular_square_pyramidal/README.md), hypothesis (H): plain Euclidean
  division in Z[λ], λ³ = λ − 3 (116 boxes), with the one-unit descent of the child.
- [mordell_k_minus17](../theorems/mordell_k_minus17/README.md), hypothesis (H), through the child (multipliers
  {1, 2} in the ring of integers of Q(∛17)).

## Provenance
- *Proved here,* for θ³=4θ−2: the ring of integer triples is integrally closed, its units are ±(θ−1)ʳ(2θ−1)ˢ, and it
  is norm-Euclidean, with every rational point within 395/512 of the lattice in norm (from the box table) and the
  point θ²/2 needing ½.
- *Known in the literature* (see the [literature check](../theorems/mordell_1963/literature_check.md)). These facts
  are classical: Mordell took the integral basis, the class number and the units from the table of Delaunay and
  Faddeev, and Billevich's 1956 table also lists a fundamental pair of units. Lemmermeyer's tables give the Euclidean minimum ½,
  which Lemmermeyer credits to Clarke (1951, not read).
- *Box covers.* Covering a fundamental domain by small boxes is the standard computational method for Euclidean
  minima: Cavallar and Lemmermeyer (1998) work in floating point with a safety margin, and Lezowski (2014) with
  proven error bounds. Here each box is checked in exact rational arithmetic, by a Taylor bound and, as a sharp
  second check, by Bernstein coefficients.
- *Trace height.* Tr(u²) appears in the unit literature as a successive-minima criterion for fundamental units
  (Godwin 1984, Cusick 1984) and as an enumeration bound (Pohst and Zassenhaus). Up to sign, the three descent
  multipliers fit the description of Berwick's (1932) extremal units recorded in the literature check (Berwick's
  paper was not read; μ₁, −μ₂, −μ₃ have norm +1, as his units do). The closest published step is in Ennola
  (1987), inside a proof by contradiction that identifies the unit minimising Tr(u²): multiplying a supposed
  other minimiser by Berwick's unit would lower Tr(u²). It is not used as a reduction. The literature check
  found no descent on Tr(u²) of the form used here.
- Louboutin and Lee (2017; only the abstract and a review were read) determine the units of this ring without
  tables, by regulator bounds.

Several sources were seen only as abstracts or citations. No novelty is claimed.
