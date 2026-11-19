# Catalogue of theorems

There is one row per theorem; follow the slug to its card, which holds the files, the Lean declarations and
the audit records. Status words are defined in [README.md](README.md). Method IDs are
explained in [methods/README.md](methods/README.md).

Rows marked "(branch review)" were consolidated on 2026-10-26 from the branches made after 2026-10-15; their
review was recorded in the originating branch (see [README.md](README.md#status-words) and
[the consolidation record](../shared/consolidation_2026-10-26/README.md)).

## Integral points: Ljunggren and Mordell equations

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [ljunggren](theorems/ljunggren/README.md) | x²+1=2y⁴, x,y>0 ⇔ (1,1),(239,13) | LEAN, two proofs | proof 1: `ghost.rational` + `lifting.somos`; proof 2: `approx.normsign` |
| [mordell_k7](theorems/mordell_k7/README.md) | x²+7=y³ ⇔ (±1,2),(±181,32) | LEAN | `ghost.rational` + `norm.pairs` + `lifting.somos` |
| [x2_x_1_y3](theorems/x2_x_1_y3/README.md) | x²+x+1=y³ ⇔ (0,1),(−1,1),(18,7),(−19,7) | WRITTEN | `ghost.rational` + `norm.pairs` + `lifting.somos` (prime 2) |
| [mordell_negative_k](theorems/mordell_negative_k/README.md) | x²+k=y³ has finitely many integer solutions for each fixed k<0 | LEAN | `ghost.unit` + `lifting.binomial` |
| [mordell_positive_k](theorems/mordell_positive_k/README.md) | x²+d=y³ has finitely many integer solutions for each fixed d>0 | LEAN | `approx.binomial` |
| [mordell_small_k](theorems/mordell_small_k/README.md) | x²+2=y³: (5,3); x²+11=y³: (4,3),(58,15) (positive solutions); procedure for x²+2h²=y³ | WRITTEN-AUDITED | `norm.pairs` |
| [mordell_1963](theorems/mordell_1963/README.md) | y(y+1)=x(x+1)(x+2) ⇔ x∈{0,−1,−2}, y∈{0,−1}; x=1, y∈{2,−3}; x=5, y∈{14,−15} (Mordell 1963, without his tables) | proof 1: WRITTEN-AUDITED; proof 2: WRITTEN | proof 1: `norm.triples`; proof 2: `ghost.rational` + `norm.triples` + `lifting.halving` |
| [lebesgue_nagell_d7](theorems/lebesgue_nagell_d7/README.md) | x²+7=yⁿ | PARTIAL | not yet classified |
| [mordell_k_minus17](theorems/mordell_k_minus17/README.md) | y²=x³+17 (x²−17=y³) ⇔ x∈{−2,−1,2,4,8,43,52,5234} (10404.f2, rank 2); eight quartic Thue equations of signature (2,1) solved completely; the 2-descent class list (H) proved without the rank or Cassels's ideals | WRITTEN (in-session review) | `ghost.unit.multi_unit` + `lifting.binomial.multi` (Prop. 17); (H): `norm.triples.multiplier` |

Together, `mordell_negative_k` and `mordell_positive_k` give finitely many solutions of x²+k=y³ for every k≠0.

## Quartic families

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [quartic_t_family](theorems/quartic_t_family/README.md) | (t+1)x⁴−t·y²=1: t=2 (Bumby) (1,1),(3,11); t=3,7,11,26 only (1,1); t=6,12: a solution of odd Pell index is (5,27), resp. (7,51); the even index and t=19 open | t=2 LEAN; t=3,7,11,26 WRITTEN-AUDITED; t=6,12 PARTIAL (odd index WRITTEN-AUDITED (branch review)); t=19 OPEN | proof 1: `ghost.rational` + `lifting.somos`; proof 2 (t=6,12, odd index): `approx.normsign` |

## Figurate numbers on rank-2 curves

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [triangular_tetrahedral](theorems/triangular_tetrahedral/README.md) | C(n,2)=C(m,3)>0 ⇔ 1, 10, 120, 1540, 7140 (1611.a1, rank 2); twelve quartic Thue equations of signature (2,1) solved completely; the 2-descent class list (H) proved without the rank | WRITTEN (in-session review) | `ghost.unit.multi_unit` + `lifting.binomial.multi` (Prop. 17); (H): `norm.triples.multiplier` |
| [triangular_square_pyramidal](theorems/triangular_square_pyramidal/README.md) | n(n+1)/2=m(m+1)(2m+1)/6>0 ⇔ 1, 55, 91, 208335 (8604.a1, rank 2); twelve quartic Thue equations of signature (2,1) solved completely; the 2-descent class list (H) proved without the rank | WRITTEN (in-session review) | `ghost.unit.multi_unit` + `lifting.binomial.multi` (Prop. 17); (H): `norm.triples` |

## Other Diophantine results

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [mixed_cubic](theorems/mixed_cubic/README.md) | x³+4y³−6hxy=1−2h³ has at most 1, 2, 0 integer pairs as h≡0,1,2 (mod 3) | WRITTEN | not yet classified |
| [x2_1250_11n](theorems/x2_1250_11n/README.md) | x²+1250=11ⁿ: only positive solution (9,3) | WRITTEN | not yet classified |
| [catalan](theorems/catalan/README.md) | a^m−b^n=1 with a,b,m,n>1 ⇒ (3,2,2,3) | OPEN (programme paused) | `ghost` (attempted) |
| [sum_and_square_sum](theorems/sum_and_square_sum/README.md) | d ≥ 8: integers y₁…y_d with Σyᵢ=s, Σyᵢ²=N exist ⇔ N ≥ M_d(s) (the balanced minimum) and N ≡ s (mod 2); answers Hoffmann's Problem 4.15 (Acta Arith. 2020) for m ≥ 8; sharp box threshold | WRITTEN-AUDITED (branch review) | not yet classified (four-square parity splitting) |
| [lucas_power_triples](theorems/lucas_power_triples/README.md) | U_n(A,±1): U_n, U_{n+rk}, U_{n+2rk} are never all r-th powers (r ≥ 2); every fixed modular sieve keeps infinitely many false positives | WRITTEN-AUDITED (branch review) | not yet classified (gap principle on r-th roots) |
| [cubic_angle_sector](theorems/cubic_angle_sector/README.md) | a ≥ 2,604,167: min \|U³−3aU²V−3UV²+aV³\| over U ≥ 2V > 0 is 8a, attained only at (3a,1); no values 0, ±1, ±2, ±3 there for a ≥ 1000 | WRITTEN-AUDITED (branch review) | `approx.normsign.all_degree` |

## Geometry and classical descents

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [integer_triangles](theorems/integer_triangles/README.md) | acute integer triangles of inradius 2p; prime-power inradius; bicentric quadrilaterals; 9(X³+Y³+Z³)=(X+Y+Z)³ | mixed: WRITTEN-AUDITED to PARTIAL | not yet classified |
| [radical_descent](theorems/radical_descent/README.md) | Fermat's square-area theorem, x⁴+y⁴=z², n²∣2ⁿ+1 ⇔ n∈{1,3}, related results | n²∣2ⁿ+1: LEAN; others written | `descent.classical` |

## Vieta jumping and conics

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [vieta_conics](theorems/vieta_conics/README.md) | x²−kxy+y²+lx+my=0 (l,m≥0): S(l,m)=D1∪D2∪T, with D1, D2 explicit divisor sets and T⊂[3, 2+(l+m)/2] found by a finite search; max S=l+m+2; the k for (a+r)/b+(b+r)/a=k, every r; 13∉S(27,0) (a published table corrected) | WRITTEN-AUDITED (branch review) | `involution.vieta` |
| [vieta_gcd_invariant](theorems/vieta_gcd_invariant/README.md) | (a+r)/b+(b+r)/a=k: r(a+b)/gcd(a,b)² = ka′b′−a′²−b′², multiplied by (k+2)/e² at each Vieta jump; r=2: 2(a+b)/gcd(a,b)²∈{1,2,3,4,5,8} (Conjecture 4.2 of Miller–Nikolakopoulos–Srinivasan) | WRITTEN-AUDITED (branch review) | `involution.vieta` |

## Quadratic forms and class numbers

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [barning_hall_class_numbers](theorems/barning_hall_class_numbers/README.md) | x²+y²=t²−m: Σ 1/\|Stab\| over the Barning–Hall chamber = (H(4m)+2H(m))/4; the number of SO⁺(2,1;Z)-orbits is (H(4m)+2H(m))/2+[m=2a²]/2+3[m=a²]/4 (OEIS A005794; the printed form conditional on two named imports) | WRITTEN-AUDITED (branch review) | reflection-group fundamental domain and Γ(2) reduction (classical; no card) |
| [parity_genus_forms](theorems/parity_genus_forms/README.md) | p≡3 (mod 4) prime: p≡11,15,23 (mod 28) ⇒ p=x²+7y²; p≡19,31 (mod 60) ⇒ p=x²+15y²; the other genera of discriminants −60 and −20; the windmills of Zagier, Jackson and Elsholtz as first returns | T1(a), T2, T3a WRITTEN-AUDITED (branch review); T1(b) the same after a repair; T7 WRITTEN | `involution.parity` |

## Matrix rings

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [matrix_cubes](theorems/matrix_cubes/README.md) | every A∈M₂(Z) is X³+Y³+Z³, and 3I is not X³+Y³, so the Waring number for cubes in M₂(Z) is 3; A is locally everywhere a sum of two cubes iff A∉O (a 3-adic class); two cubes over Z open (Conjecture 3.6) | LEAN (three cubes, w₃=3, the 3-adic obstruction; branch build record); local theorem WRITTEN-AUDITED (branch review) | `matrix.cayley_hamilton` |

## Additive problems: octahedral numbers

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [octahedral_bounded_sums](theorems/octahedral_bounded_sums/README.md) | Every N≥0 is a sum of eight octahedral numbers O(n)=(2n³+n)/3, given the three-square theorem for m≡1 (mod 4) (Dirichlet 1850); for N≤2.06·10⁴⁰ with no hypothesis (from ledger W011). Nine terms with a Lean-checked large range (N≥4,395,799,615) | WRITTEN-AUDITED (branch review) | not yet classified |

The lane objective, seven octahedral numbers for every N, remains open: see [research/waring/](../waring/problem.md).

## Determinants of Legendre symbols, and congruences from finite identities (cloud Claude branch, unreviewed)

Imported on 2026-10-26 from `origin/claude/math-prime-cross-domain-g1vwsn` (tip 26c718ac), tracks 12, 16 and 17. None of
these concerns the Eureka or Waring objectives. No independent review is recorded for any of them.

| Theorem | Statement | Status | Methods |
|---|---|---|---|
| [legendre_hankel_determinants](theorems/legendre_hankel_determinants/README.md) | Sun's Conjectures 3.2(i), 3.5(i), 3.6(i), 4.2, 4.3, 4.4, 4.5, 4.6 (arXiv:2405.03626) on det[x + ((j+k+c)/p) + (j/p)y + (k/p)z + (jk/p)w] over Hankel windows, and track 03's Conjecture 4.7; p ≡ 3 (mod 4) elementary (the h(−p)-signs via Mordell's theorem), p ≡ 1 (mod 4) via Dirichlet's class number formula and Gauss's sign | WRITTEN (unreviewed; checked numerically by the 2026-10-26 consolidation) | not yet classified |
| [truncated_binomial_valuation](theorems/truncated_binomial_valuation/README.md) | ν_p((1+y)^N − Σ_{m<k} C(N,m)y^m) = kν_p(y) + ν_p(C(N,k)) when ν_p(C(k+i,i)) < iν_p(y) for all i ≥ 1; for k = p^s it is p^s ν_p(y) + ν_p(⌊N/p^s⌋); the 2-adic valuation of Σ_{r=j}^{n} (−1)^r C(r,j) | WRITTEN (unreviewed; checked numerically by the 2026-10-26 consolidation) | not yet classified |
| [wolstenholme](theorems/wolstenholme/README.md) | H_{p−1} ≡ 0 (mod p²), H^{(2)}_{p−1} ≡ 0 (mod p), C(2p−1,p−1) ≡ 1 (mod p³) for p ≥ 5, from one interpolation identity, without Fermat's little theorem | WRITTEN (unreviewed; checked numerically by the 2026-10-26 consolidation) | not yet classified |
| [jacobsthal_kazandzidis](theorems/jacobsthal_kazandzidis/README.md) | C(ap,bp) ≡ C(a,b) (mod p^{3+ν_p(ab(a−b))}) for p ≥ 5, without the p-adic Gamma function | WRITTEN (unreviewed; checked numerically by the 2026-10-26 consolidation) | not yet classified |
| [martin_ruiz_congruences](theorems/martin_ruiz_congruences/README.md) | Σ_{a<p} H_a q_p(a) ≡ (Q₁² + Q₂)/2 − Q₁ (mod p); ((p−1)/2)! ≡ (−1)^{(p−1)/2} Σ_{i≤(p−1)/2} C(2i,i)4^{−i}(i/p) (mod p) | WRITTEN (unreviewed; checked numerically by the 2026-10-26 consolidation) | not yet classified |

## Lane objectives (open)

| Objective | Where |
|---|---|
| Eureka: every n≥0 is a sum of three triangular numbers, by a new elementary descent | [research/eureka/](../eureka/problem.md) |
| Waring: every N≥0 is a sum of seven octahedral numbers | [research/waring/](../waring/problem.md) |
