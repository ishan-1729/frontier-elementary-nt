# involution.parity: parity of fixed points of side involutions on hyperbolic polygons
Parent: involution     Children: none

## What it does
1. The binary forms of a given discriminant in a level lattice are written as integer triples, for example
   N x² + 4yz = p. A form of positive discriminant is a geodesic of the hyperbolic plane.
2. Take a convex hyperbolic polygon whose sides lie on mirrors of integral automorphs, or are folded by an integral
   Atkin–Lehner half-turn. Each side gives an explicit linear involution of the finite set of triples whose geodesic
   crosses that side.
3. The incidence-parity lemma (the sign changes of a cyclic sequence are even) gives
   Σ #Fix(ιⱼ) + Σ |A_{Vᵢ}| ≡ 0 (mod 2). The fixed points of the side involutions are representations by an
   indefinite form, or trivial factorisations. The corner terms are representations by the definite form.
4. Congruences on p make every term but one vanish. That term is odd, so the representation exists.

Using every side of a quadrilateral or a hexagon, instead of a first-return map on a triangle, keeps every
involution linear. On triangles the first-return maps are the windmills of Zagier, Jackson and Elsholtz.

## Ingredients (with the file that proves each)
All in [parity_genus_forms/source/06_involutions_and_parity/theorems/](../theorems/parity_genus_forms/source/06_involutions_and_parity/theorems/):
- the incidence-parity lemma, known in substance: `T0_incidence_parity_lemma.md`;
- the polygon set-ups for 7x² + 4yz = p and 15x² + 4yz = p, with their symbolic identities: `T1_x2_plus_7y2.md`,
  `T2_x2_plus_15y2_genus.md`;
- the transfer identity that measures the limit: `T3_obstruction_p1mod4.md`;
- the windmills as first returns: `T7_windmills_are_first_returns.md`.

## Requirements: when it applies
- The real quadratic field of the irrational mirror must not split p: then its fixed-point term vanishes by
  congruences. This is why only the p ≡ 3 (mod 4) halves of the classical statements are reached.

## Limits: where it stops, with examples
- p ≡ 1 (mod 4): for p ≠ 5, N(x² + 5y²) ≡ #{X² − 20Y² = p, 5Y² < p} (mod 2), and the real quadratic term cannot be
  removed. So p ≡ 1, 9 (mod 20) ⇒ p = x² + 5y² is not reached. The same happens for x² + 2y² with p ≡ 1 (mod 8),
  the case Elsholtz left open.
- Parity proves existence only where a count is odd. For three triangular numbers the ordered count t(k) is
  ≡ d(8k+3)/2 (mod 2), odd only on a set of density zero (Cooper–Riasanovsky 2013, Prop. 26; the branch's track 06,
  T4). So plain parity cannot prove Eureka; compare register entry R9 in
  [research/eureka/MECHANISMS.md](../../eureka/MECHANISMS.md).

## Theorems using it
- [parity_genus_forms](../theorems/parity_genus_forms/README.md).

## Provenance
- Imported from branch A (`origin/claude/math-prime-cross-domain-g1vwsn`, tip `26c718ac`, track 06).
- The parity principle and the windmill involutions are due to Zagier (1990), Heath-Brown (1984), Jackson (2000,
  2001) and Elsholtz (2010). Every involution proof found treats a form over a field of class number 1:
  x² + y², x² + 2y², x² + 3y² and 3x² + 4y².
- New here: the polygons with mirror and Atkin–Lehner sides, which give the first parity proofs found for x² + 7y²
  and over an imaginary quadratic field of class number 2 (Q(√−15), and Q(√−5)). Searches: the branch's log
  Q1–Q31, and the consolidation's
  [novelty_web_log.md](../../../research/shared/consolidation_2026-10-26/a_tracks/work/novelty_web_log.md), Q5.
  Jackson's 2000 Tatra paper and Heath-Brown's 1984 note were not read.
