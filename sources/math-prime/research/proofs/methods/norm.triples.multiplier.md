# norm.triples.multiplier: division with a multiplier, for a ring that is not norm-Euclidean
Parent: [norm.triples](norm.triples.md)     Children: none yet

**Status: WRITTEN (2026-10-15)**, as part of [triangular_tetrahedral](../theorems/triangular_tetrahedral/README.md),
in [hypothesis_H.md](../theorems/triangular_tetrahedral/hypothesis_H.md), its first use. Refereed in-session only.

## What it does
It replaces the Euclidean division of `norm.triples` when the ring of triples is not norm-Euclidean.
- **Division with a multiplier.** For every rational triple ξ outside the ring there are a multiplier s from a finite
  set and κ in the ring with 0 < |N(sξ − κ)| < 1. For α, β with β ∤ α this gives 0 < |N(sα − κβ)| < |N(β)|; this is
  the Dedekind–Hasse criterion.
- **Bézout by least norm.** Among the nonzero combinations xα + yβ, one of least |N| divides α and β.
- **The rest of `norm.triples`.** From Bézout on, the proofs are unchanged: Euclid's lemma, prime elements of prime
  norm, exact powers, coprime square allocation, integral closedness.
- **Units of a complex cubic.** One unit ε, by descent on T₂ = σ₁² + 2|σ₂|². If neither εu nor ε⁻¹u has smaller T₂,
  then T₂(u) is bounded; an enumeration of the units below the bound finishes. The parent card was written for the
  totally real case; this one-unit descent is also used with a plain division, for λ³ = λ − 3
  ([triangular_square_pyramidal](../theorems/triangular_square_pyramidal/hypothesis_H.md) §4).

## Ingredients (with the file that proves each)
All in [theorems/triangular_tetrahedral/](../theorems/triangular_tetrahedral/README.md), for λ³ = 4λ − 6
(discriminant −716):

| Ingredient | Proof | Check |
|---|---|---|
| Ring, norm, adjugate; the norm vanishes only at 0 (2-adic descent) | hypothesis_H.md §1 | [04_ring.py](../theorems/triangular_tetrahedral/scripts/04_ring.py) |
| Division with multipliers {1, 2}: 398 boxes, Taylor bounds ≤ 8191/8192, sharp maximum 4093/4096 | §2 | [05_dh_cover.py](../theorems/triangular_tetrahedral/scripts/05_dh_cover.py) |
| The ring is not norm-Euclidean: ξ* is fixed by ε, and every translate has \|N\| ≥ 121/109 | Remark 2.5 | 05_dh_cover.py |
| Bézout, primes, square allocation, integral closedness | §3 | (logic only) |
| Units ±εᵏ, by the T₂ descent (threshold 114) | §4 | [06_units.py](../theorems/triangular_tetrahedral/scripts/06_units.py) |
| Second proof of unique factorization, by Dedekind's box principle (constant 535/24) | §6 | [08_ufd_box.py](../theorems/triangular_tetrahedral/scripts/08_ufd_box.py) |

## Requirements: when it applies
- A finite multiplier set S and a box table. A table exists when every point of the closed fundamental domain has
  some s ∈ S and a translate with |N(sξ − κ)| < 1, by compactness. A search finds it
  ([cover_search_K.py](../theorems/triangular_tetrahedral/scripts/cover_search_K.py)).
- The points ξ ∉ R with sξ ∈ R need a translate with another multiplier. For S = {1, 2} these are the seven
  half-lattice classes, each with an explicit κ.
- A ring that satisfies the division is a principal ideal domain. So the method can work only when the class number
  is 1. It certifies that case by case, and needs no class-number theory.

## Limits: where it stops, with examples
- The table is a computer-found certificate, checked exactly. For λ³ = 4λ − 6 the margins are tight: the largest
  Taylor bound is 8191/8192.
- Whether a given set of multipliers works is decided only by the certificate, including the points ξ ∉ R with
  sξ ∈ R. Integer multipliers need not give a Dedekind–Hasse division even when the class number is 1: in
  Z[(1+√−19)/2] no integer s works at ω/2. An odd s fixes the coset of a half-lattice point, and an even s lands in R.
- Hurwitz (1919; not accessed, via Lemmermeyer) gives, for every field, an integer j < m with |N(jx − y)| < 1, where
  jx = y is allowed. That is weaker.
- The unit descent, as written, is for one unit (signature (1,1)). Totally real cubics have the three-multiplier
  descent of `norm.triples`.

## Theorems using it
- [triangular_tetrahedral](../theorems/triangular_tetrahedral/README.md): hypothesis (H), which makes Theorem B
  unconditional.
- [mordell_k_minus17](../theorems/mordell_k_minus17/README.md): hypothesis (H) for y² = x³ + 17, in the ring of
  integers Z + Zθ + Zω of Q(∛17) (θ³ = −17, ω = (1 + θ + θ²)/3; discriminant −867). The ring is not norm-Euclidean (a
  point fixed by −η needs 1115/1028); 795 boxes with multipliers {1, 2}, the seven half-lattice classes, and units
  ±ηᵏ with η = 18 + 7θ ([hypothesis_H.md](../theorems/mordell_k_minus17/hypothesis_H.md)).
- [triangular_square_pyramidal](../theorems/triangular_square_pyramidal/README.md): the one-unit descent only, for
  λ³ = λ − 3; its ring is norm-Euclidean, so the division there is that of `norm.triples`.

## Provenance
See [literature_check.md](../theorems/triangular_tetrahedral/literature_check.md) §5.
- The Dedekind–Hasse criterion: Dedekind (Werke II 373–375, from the Nachlass; not accessed, cited by Lemmermeyer),
  and Hasse (J. reine angew. Math. 159, 1928). Hurwitz (Math. Z. 3, 1919; not accessed) allows jx = y; see Limits.
- Hand coverings for Z[(1+√−19)/2] use the multipliers 1 and 2, and θ or θ̄ at the 2-torsion points: Cámpoli (Amer.
  Math. Monthly 95, 1988); R. A. Wilson (2016).
- Computer coverings for Euclidean minima: Cavallar–Lemmermeyer (1998), Lezowski (2014). Cavallar–Lemmermeyer list
  −716 as "N 121/109" and −867 as "N 1115/1028", which match Remark 2.5 of the two cards.
- Two-stage divisions are known for other complex cubic fields (−199, −307, −335, −351).
- Not found: a published multiplier certificate for −716 or −867, or for any cubic field.
- No novelty is claimed.
