# ghost.rational: rational ghosts on an enlarged genus-1 equation
Parent: [ghost](ghost.md)     Children: `ghost.rational.multi_index` (open, below)

## What it does
The original equation is the slice "a coordinate equals 1" of a homogeneous equation. The
primitive integer points of that equation, the nodes, are the rational points of a genus-1 curve.
- **Ljunggren:** the nodes are r⁴+q²=2s⁴ with gcd(r,s)=1, and the original solutions are the nodes with r=1.
- **x²+7=y³:** the ghosts are the reduced rational points (A/B², C/B³) of Y²=X³−7, and the original
  solutions are those with B=1.
- **x²+x+1=y³:** the same, on Y²+Y+1=X³.
- **y(y+1)=x(x+1)(x+2)** (proof 2): the reduced rational points (A/B², C/B³) of y²+y=u³−u with u=x+1, and the
  original solutions are those with B=1.

The proof then runs in four steps:
1. **Descent.** Every node descends, by a strictly decreasing integer height, to finitely many seeds.
   - Ljunggren uses Pythagorean parametrizations and the Vierzahlensatz, with height s.
   - x²+7=y³ uses two 3-isogeny maps, whose parents exist by `norm.pairs`, with height A.
   - x²+x+1=y³ uses two 3-isogeny maps through y²+y=x³+20, with height A. The parents come from Eisenstein
     pairs (`norm.pairs`) and from integer cube allocation, because the auxiliary curve is (y−4)(y+5)=x³.
   - y(y+1)=x(x+1)(x+2) has no isogeny to use (37a1 is alone in its isogeny class), so it halves: a point of the
     unbounded real component is D(R′) for a rational R′, by a square root in the ring of triples (`norm.triples`),
     and a point of the bounded component is first moved across by the chord S. The height is max(|A|, B²), and
     doubling cancels at most 37².
2. **One signed index.** Explicit rational maps on the cubic turn the descent tree into one
   signed index l. These are chord and tangent formulas, checked as polynomial identities with no group law,
   together with an intertwining identity: DA=A²D for Ljunggren, TS=S³T for k=7 and for x²+x+1=y³, DS=S²D for
   y(y+1)=x(x+1)(x+2).
3. **Fixed-value marker.** An integer Somos-4 recurrence carries the reduced numerators and
   denominators along the index. A primitive-numerator lemma turns the marker into a fixed value:
   f_l=−1 for Ljunggren, |h_n|=1 for k=7 and for x²+x+1=y³, a_l²=C² for the quartic family. For
   y(y+1)=x(x+1)(x+2) the marker B=1 is used directly, without a recurrence.
4. **Exclusion.** [lifting.somos](lifting.somos.md) excludes the marker at every index except the seeds. When the
   real curve has a bounded component that carries the odd indices, [lifting.halving](lifting.halving.md) does it
   without a prime (y(y+1)=x(x+1)(x+2)).

## Requirements
The rational points must collapse onto one signed index. Every completed case has
Mordell–Weil rank 1; both stuck cases have rank 2 (LMFDB, read 2026-10-14):

| Case | Curve | LMFDB | Mordell–Weil group | Status |
|---|---|---|---|---|
| ljunggren | y²=x³−2x | 256.b1 | Z ⊕ Z/2 | LEAN |
| mordell_k7 | y²=x³−7 | 5292.i1 | Z | LEAN |
| x2_x_1_y3 | y²+y=x³−1 | 243.a1 | Z | WRITTEN |
| mordell_1963 (proof 2) | y²+y=x³−x | 37.a1 | Z | WRITTEN |
| quartic t=2 (Bumby) | y²=x³−12x | 288.a1 | Z ⊕ Z/2 | LEAN |
| quartic t=3 | y²=x³−36x | 576.c3 | Z ⊕ (Z/2)² | WRITTEN-AUDITED |
| quartic t=7 | y²=x³−392x | 12544.f1 | Z ⊕ Z/2 | WRITTEN-AUDITED |
| quartic t=11 | y²=x³−1452x | 34848.c1 | Z ⊕ Z/2 | WRITTEN-AUDITED |
| quartic t=26 | y²=x³−18252x | 48672.c1 | Z ⊕ Z/2 | WRITTEN-AUDITED |
| quartic t=6 (7x⁴−6y²=1) | y²=x³−252x | 14112.c1 | Z² ⊕ Z/2 | OPEN |
| quartic t=12 (13x⁴−12y²=1) | y²=x³−117x | 48672.b1 | Z² ⊕ Z/2 | OPEN |

How the quartic-family curves are obtained:
- For (t+1)x⁴−t·y²=1, the curve listed is y²=x³−t²(t+1)x.
- Put v=ty and u=x. Then v²=au⁴+b with a=t(t+1), b=−t, and the map (X,Y)=(au², auv) sends this
  quartic onto y²=x³+abx.
- Curves related by scaling or isogeny have equal rank. For t=12, y²=x³−1872x is
  isomorphic to y²=x³−117x.

## Limits
At rank 2 the ghosts carry two independent indices. The marker becomes one condition on a
two-index lattice, and no all-index exclusion is known for that. The t=6 handover found the
primitive ghost node (1189,1475,2277967) outside both analysed trajectories
([handover](../../shared/prior/ghost_parametrization/ghost_parametrization_handover/HANDOVER_PROMPT.md)),
which is what a second index looks like.

## Generalization slot: `ghost.rational.multi_index` (OPEN)
The goal is an all-index exclusion of a fixed-value marker on a rank-2 lattice of rational ghosts. It would
reach the t=6 and t=12 members of [quartic_t_family](../theorems/quartic_t_family/README.md) and
several figurate-number problems (see [targets/](../targets/README.md)).

## Theorems using it
[ljunggren](../theorems/ljunggren/README.md) (proof 1), [mordell_1963](../theorems/mordell_1963/README.md) (proof 2),
[mordell_k7](../theorems/mordell_k7/README.md), [quartic_t_family](../theorems/quartic_t_family/README.md),
[x2_x_1_y3](../theorems/x2_x_1_y3/README.md).

## Provenance
Appendix B of the [ghost-descent Ljunggren proof](../theorems/ljunggren/ghost_descent/ljunggren_elementary_proof_descent.pdf) records:
- the signed quartic descent, attributed to Conrad in inherited working notes;
- Elkies (2014, online), who suggested pairing the cubic with an obstruction modulo powers of 3;
- Hone–Swart (2008), whose invariant specializes to the Somos recurrence used.

No priority has been established for the complete argument.

For [mordell_1963](../theorems/mordell_1963/README.md) (proof 2), the halving descent and the exclusion follow the route
of Silverman's AEC Ex. 10.9 and 9.13 (see its
[literature check](../theorems/mordell_1963/ghost/literature_check_ghost.md)); no novelty is claimed there.
