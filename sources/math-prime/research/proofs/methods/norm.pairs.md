# norm.pairs: Euclidean arithmetic of integer pairs
Parent: norm     Sibling generalization: [norm.triples](norm.triples.md)

## What it does
Pairs (a,b) with an explicit multiplication rule stand in for a+bω, for example with ω²=−2, or with
norm a²+ab+kb². Every property used is derived directly:
- conjugation and a multiplicative norm;
- Euclidean division, by explicit rounding regions whose error norm stays below 1;
- gcd and Bézout identities;
- the units, found by descent;
- coprime cube allocation: UV=W³ with U and V coprime gives U = unit · cube.

No unique-factorization or class-number theorem is assumed.

## Theorems using it
- [mordell_small_k](../theorems/mordell_small_k/README.md).
- [mordell_k7](../theorems/mordell_k7/README.md), for the existence of parents under both isogeny maps.
- [x2_x_1_y3](../theorems/x2_x_1_y3/README.md), for the Ψ-parents. These use Eisenstein pairs (ω²=−1−ω), where
  division has error norm at most 7/16 and the cube rule leaves the unit ω modulo cubes. The Φ-parents need only
  integer cube allocation, because the auxiliary curve factors over Z as (y−4)(y+5)=x³.

## Limits
It needs a norm-Euclidean pair ring. At c=−23 the cube extraction fails: 2+√−23 has norm 27 and an
explicit pair Bézout relation with its conjugate, yet it is not a pair cube. See the 2026-09-24 entries in
[CURRENT_STATE.md](../../shared/prior/ghost_parametrization/team_work/CURRENT_STATE.md).

## Generalization: [norm.triples](norm.triples.md) (done for θ³=4θ−2)
The same toolkit for triples a+bθ+cθ² in a cubic ring:
- explicit Euclidean division with error norm below 1;
- gcd and Bézout identities;
- units modulo squares, found by descent;
- square allocation.

It replaces the tables of cubic-field invariants that many classical solutions cite. The first case done is Mordell's
y(y+1)=x(x+1)(x+2), with θ³=4θ−2 and discriminant 148: [mordell_1963](../theorems/mordell_1963/README.md) (target in
[targets/](../targets/ghost_screen_2026-10-14.md)).

## Sources
- [mordell_descent_example_minus2.md](../../shared/prior/ghost_parametrization/team_work/team2/x2plus7/mordell_descent_example_minus2.md) (x²+2=y³).
- [main.tex §seven](../../shared/prior/elementary_diophantine_paper/main.tex), the Euclidean lemma with three regions.
