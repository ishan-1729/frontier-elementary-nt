# lifting.somos: exact lifting on Somos-4 sequences
Parent: lifting     Children: none yet

## What it does
It excludes a fixed value of an integer observable of a Somos-4 sequence, at every index.
1. Scale the sequence so that its transition becomes a polynomial map with integer coefficients, some
   power of which is the identity modulo p.
2. Finite differences of the observable gain one factor of p per order, and the finite binomial theorem
   controls long jumps.
3. A nonzero first difference at the seed residue classes then rules out every other index, at a
   precision chosen per index.

| Theorem | Prime | Exact statement used |
|---|---|---|
| ljunggren (proof 1) | 3 | G(j)=b_{2j}b_{2j−2}−9b_{2j−1}² equals 4 only for j=0,1 (long-jump congruence modulo 3^{e+3}) |
| mordell_k7 | 5 | v₅(s_{c+5q}−s_c)=2+v₅(q) for c=0,−1, so \|s_k\|=1 only at k=0 |
| quartic t=7 | 2 | v₂(a_n−a_k)=v₂(n−k)+3 for n≡k (mod 8) |
| x2_x_1_y3 | 2 | v₂(s_{c+8q}−s_c)=4+v₂(q) for c=0,−1, so \|s_k\|=1 only at k=0 |

## Requirements
The method needs a prime for which some power of the scaled transition is congruent to the identity, and a
nonzero leading increment at the seed classes. The prime is chosen case by case.
- **At p=2** the doubling step J²(z)=z+2^{a+1}Λ(z)+2^{2a}Θ(z), for J(z)=z+2^aΛ(z), needs a ≥ 2. So the period map
  must be ≡ id modulo 4, not just modulo 2 ([x2_x_1_y3](../theorems/x2_x_1_y3/proof.md), Lemma 7.3).
- **A prime dividing h₃ can still fail.** For x²+x+1=y³ the prime 3 gives every s_k ≡ 1 (mod 9). But the class
  k ≡ 5 (mod 9) has s_k ≡ 1 (mod 3⁵) at every k < 3600 checked, which looks like a 3-adic zero away from the
  known indices. The prime 2 works ([08_crosscheck.py](../theorems/x2_x_1_y3/scripts/08_crosscheck.py), 8.3).

## Theorems using it
[ljunggren](../theorems/ljunggren/README.md) (proof 1), [mordell_k7](../theorems/mordell_k7/README.md),
[quartic_t_family](../theorems/quartic_t_family/README.md), [x2_x_1_y3](../theorems/x2_x_1_y3/README.md).

## Sources
- §5 of the [ghost-descent Ljunggren proof](../theorems/ljunggren/ghost_descent/ljunggren_elementary_proof_descent.pdf).
- [main.tex §seven](../../shared/prior/elementary_diophantine_paper/main.tex), "Finite lifting".
- [ROOT_T7_UNIFORM_ASSEMBLY.md](../../shared/prior/ghost_parametrization/team_work/team3/ROOT_T7_UNIFORM_ASSEMBLY.md).
- Section 7 of the [x2_x_1_y3 proof](../theorems/x2_x_1_y3/proof.md) (the prime 2, period 8).

## Provenance
- **Poonen (2014).** His p-adic interpolation of iterates uses the same substitution-difference mechanism
  analytically. The finite statements are proved directly here.
- **Ayad (1991–93), the closest prior use, on other curves.** He writes integrality of mP as ψ_m = ±1 and lifts
  congruences of the division-polynomial sequence modulo prime powers:
  - J. Number Theory 38 (1991): at p = 2, for ψ₃ odd;
  - Manuscripta Math. 76 (1992) and Ann. Inst. Fourier 43 (1993): at odd p.

  [x2_x_1_y3](../theorems/x2_x_1_y3/literature_check.md) treats the case 2 | ψ₃ that he left open.
