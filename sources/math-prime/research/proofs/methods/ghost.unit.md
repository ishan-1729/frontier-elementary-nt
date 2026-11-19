# ghost.unit: integer ghosts on a norm form with one unit
Parent: [ghost](ghost.md)     Children: [ghost.unit.multi_unit](ghost.unit.multi_unit.md)

## What it does
Take C²=A³+d with d>0, that is x²+k=y³ with k=−d<0.
1. **Finitely many cubic equations.** A pigeonhole search for short pairs, plus Bézout, puts every
   solution under one of finitely many labels. Solutions with the same label differ by a cube, which gives
   finitely many binary cubic equations F(v,t)=q³ with q²∣4d.
2. **Enlargement.** Each irreducible cubic is enlarged to integer triples z with N(z)=m in the system
   ϑ³=3Aϑ+2C. The original solutions are the triples with w=0.
3. **Multiplier.** A norm-one multiplier ε is built by pigeonhole, with M(ε)≡I (mod 4).
4. **Descent.** Outside a box of size O(|m|^{2/3}), the height H=u²+v²+w² at least halves at a parent. So every
   ghost is Tⁿβ for one of finitely many seeds β.
5. **Count.** [lifting.binomial](lifting.binomial.md) bounds how many n give w=0 on each orbit.

The output is finiteness with a count per orbit, not the positions of the solutions, and no algorithm.

## Requirements
Exactly one independent unit. In the Mordell reduction the cubic X³−3AX−2C has discriminant
−108d<0, hence one real root and one unit.

## Limits
With two independent units, w=0 is one condition in two exponents, and the 2-adic count no longer isolates
zeros. This happens for a totally real cubic, as in x²+k=y³ with k>0. The general ghost route for k>0 stayed open. That theorem was proved by
[approx.binomial](approx.binomial.md) instead.

## Generalization: [ghost.unit.multi_unit](ghost.unit.multi_unit.md) (WRITTEN, 2026-10-15)
Two or more units, with as many markers as units.
- **Coverage: done.** A least-E element of a unit orbit lies in a polyhedral cone cut out by finitely many units, and
  the vertices of that cone bound its size. This is an explicit certificate, with no unit theorem. With c = ±1 it
  also proves that two given units are fundamental.
- **Exclusion: done.** It uses the Jacobian criterion of [lifting.binomial.multi](lifting.binomial.multi.md)
  (Proposition 17), applied on the classes of the kernel lattices of the units modulo prime powers, with a
  congruence sieve.
- **First use:** [triangular_tetrahedral](../theorems/triangular_tetrahedral/README.md), twelve quartic Thue
  equations of 1611.a1. Later uses: [triangular_square_pyramidal](../theorems/triangular_square_pyramidal/README.md)
  (8604.a1) and [mordell_k_minus17](../theorems/mordell_k_minus17/README.md) (y² = x³ + 17).
- **Where it applies.** Only to Thue fields with a complex place. The four rank-2 targets have none; the curves
  1611.a1, 8604.a1 and y²=x³+17 do ([rank2_signatures_2026-10](../targets/rank2_signatures_2026-10.md)).

## Theorems using it
[mordell_negative_k](../theorems/mordell_negative_k/README.md).

## Provenance
The team's own audit calls the written argument an elementary reconstruction of Skolem's classical
p-adic mechanism, done with finite polynomials. See
[ROOT_SPECIFIC_SKOLEM_ELEMENTARITY_AUDIT_20260924.md](../../shared/prior/ghost_parametrization/team_work/team2/ROOT_SPECIFIC_SKOLEM_ELEMENTARITY_AUDIT_20260924.md)
and [ROOT_POSITIVE_MORDELL_FINITE_WRITTEN.md](../../shared/prior/ghost_parametrization/team_work/team2/ROOT_POSITIVE_MORDELL_FINITE_WRITTEN.md).
In those files "positive c" means Y²=X³+c.
