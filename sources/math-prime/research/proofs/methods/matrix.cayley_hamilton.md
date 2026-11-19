# matrix.cayley_hamilton: Cayley–Hamilton linearisation of powers in M₂(Z)
Parent: matrix     Children: none

## What it does
For a 2×2 matrix X with trace t and determinant d, Cayley–Hamilton gives X² = tX − dI and
X³ = (t² − d)X − td·I. A power of X is thus linear in X, with coefficients in Z[t, d], and a Waring-type problem in
M₂(Z) becomes linear algebra with arithmetic coefficients.
1. *Linear cubes.* X³ = X holds for rank-one idempotents (t = 1, d = 0), their negatives, and involutions
   (t = 0, d = −1). Sums of two of them realise every off-diagonal pair once the diagonal is one of five pairs.
2. *Shifted cube.* On the plane tr W = 3c, (W − cI)³ = (3c² − det W)·W − c³I, so one cube is a scalar multiple of W
   up to a scalar matrix.
3. *Arithmetic.* What remains is a congruence for one diagonal entry of an SL₂(Z)-conjugate of the target. It is
   solved with the binary form v ↦ det(v, Bv) of the conjugacy class, which takes values prime to any modulus when
   B is primitive, and a 2-adic valuation argument.
4. *Two cubes.* A = X³ + Y³ becomes a linear relation in the lattice (sl₂(Z), −det). Locally it is linear algebra
   over finite fields, through the differential H ↦ Y²H + YHY + HY² of Y ↦ Y³, with Hensel lifting.

## Ingredients (with the file that proves each)
In [matrix_cubes/source/04_matrix_waring_cubes/theorems/](../theorems/matrix_cubes/source/04_matrix_waring_cubes/theorems/):
- Lemmas 1–4 of `T2_three_cubes.md`: the linear cubes, the shifted cube, conjugation, the 2-adic choice;
- Lemmas D, H and H′ of `T4_local_two_cubes.md`: invertibility of the differential, and Hensel lemmas;
- Proposition 3.2 of `T3_two_cubes.md`: the necessary conditions for two cubes.

## Requirements: when it applies
- 2×2 matrices over Z, where conjugacy classes correspond to binary quadratic forms.
- A usable shift: for cubes, c with tr W = 3c.

## Limits: where it stops, with examples
- Two cubes over Z: the method gives necessary conditions and a complete search for bounded traces, but no bound on
  the traces of a representation. So Conjecture 3.6, that (−6 −5; −1 −1) is not a sum of two cubes, stays open
  ([matrix_cubes](../theorems/matrix_cubes/README.md)).
- Fourth powers in M₂(Z), and n × n matrices with n ≥ 3, were not tried.

## Theorems using it
- [matrix_cubes](../theorems/matrix_cubes/README.md).

## Provenance
- Imported from branch A (`origin/claude/math-prime-cross-domain-g1vwsn`, tip `26c718ac`, track 04).
- For two squares in M₂(Z), the same linearisation with the same trace choices is Kaneyama's (2013, Lemme 7). The
  branch's two-square theorem is a condensed version of his and is not imported.
- For cubes, the shifted cube on tr W = 3c, combined with idempotent and involution cubes, was not found. Searches:
  the branch's LITERATURE.md, and the consolidation's
  [novelty_web_log.md](../../../research/shared/consolidation_2026-10-26/a_tracks/work/novelty_web_log.md), Q1.
- The correspondence between GL₂(Z)-conjugacy classes of integer matrices and binary quadratic forms is classical
  (Latimer–MacDuffee). The value lemma is, for example, Lemma 2.25 of Cox's *Primes of the form x² + ny²*; the
  branch re-proves it.
