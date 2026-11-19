# Exact valuation of a truncated binomial expansion ("higher-order LTE"), and the 2-adic valuation of alternating column sums of Pascal's triangle
Status: WRITTEN (unreviewed; checked numerically by the 2026-10-26 consolidation)     Methods: not yet classified
Generalizes: the lifting-the-exponent lemma (the case k = 1)     Generalized by: none

Origin: branch `origin/claude/math-prime-cross-domain-g1vwsn`, tip 26c718ac, folder
`cloud_agent_claude/16_factorial_splitting_generalized/` (cloud Claude programme, track 16, 18–19 Oct 2026), files
T2, T3 and T4 of that folder, in [source/](source/). The track started from Ishan's "AB identity"; his own note and
image are private and are not reproduced here. No independent review is recorded.

## Statement
Notation: p is a prime, ν_p the p-adic valuation, and for integers N ≥ k ≥ 1 and y in a ring
W_{N,k}(y) := Σ_{i=0}^{N−k} C(N, k+i) y^i, so that (1+y)^N − Σ_{m<k} C(N,m) y^m = y^k W_{N,k}(y).
Condition (★) for (p, ε, k, ν): ε·ν_p(C(k+i, i)) < i·ν for every i ≥ 1. K_p(ν) is the set of k with (★) for ε = 1.

1. **Valuation theorem** (T2 Thm 2.1). Let O be a discrete valuation ring of characteristic 0 with residue
   characteristic p, v its normalised valuation, ε = v(p), y ∈ O with ν = v(y) ≥ 1, and assume (★). Then for every
   N ≥ k: v((1+y)^N − Σ_{m<k} C(N,m)y^m) = kν + ε ν_p(C(N,k)); equivalently W_{N,k}(y) = C(N,k)(1+z) with v(z) ≥ 1.
   The same holds for N ∈ Z_p ∖ {0, …, k−1} when O is complete, and in the "deformed hockey stick" form
   v(Σ_{r=k−1}^{N−1} C(r,k−1)x^r) = ε ν_p(C(N,k)) when v(x−1) = ν ≥ 1.
2. **Powers of p** (T2 Cor 3.2). For s ≥ 0, (p, s) ≠ (2, 0), every y ∈ Z_p with ν_p(y) ≥ 1 and every N ≥ p^s:
   ν_p((1+y)^N − Σ_{m<p^s} C(N,m)y^m) = p^s ν_p(y) + ν_p(⌊N/p^s⌋). (s = 0 is the LTE lemma.)
   Also: one-digit k = c·p^s (Cor 3.3), and the uniform criterion p^{ν/ε} > k+1 ⇒ (★), which is sharp (Cor 3.4).
3. **Exact range of validity** (T2 Thm 3.5–3.7, Prop 3.8). A digit description of K_p(ν) (carry runs); for p = 2, (★) is
   necessary and sufficient for each fixed y; for odd p, failure is strict when min_i(iν − ν_p(C(k+i,i))) < 0, ties are
   decided by a finite criterion, ties at i = 1 always fail, and (★) is not necessary (p = 3, k = 7).
4. **Alternating column sums, 2-adically** (T3 Thm 3.1, Cor 3.2). For A(n,j) := Σ_{r=j}^{n} (−1)^r C(r,j):
   ν_2(A(n,j)) = ν_2(C(n+1,j+1)) for every n ≥ j if and only if j+1 ∈ K_2(1) (j = 1, 3, 7, 9, 11, 15, 17, …); for
   j = 2^s − 1, ν_2(A(n,j)) = ν_2(⌊(n+1)/2^s⌋). (T3 Cor 3.2: the same at any x ∈ Z_2 with x ≡ 1 mod 2^ν, for j+1 ∈ K_2(ν).)
5. **Column 5 is governed by √41** (T3 Prop 3.3). With α ∈ Z_2 the root of 4t² − 3t − 2 (α = (3 − √41)/8, α even):
   ν_2(A(n,5)) − ν_2(C(n+1,6)) = ν_2(⌊n/2⌋ − λ_n), λ_n = α (n even), 1 − α (n odd). So the exceptional residue classes
   follow the binary digits of a quadratic irrational. T3 Prop 3.4: reflection symmetry and double roots in even
   columns. T3 Obs 3.5 (computation only): irrational 2-adic roots in columns 5, 13, 21, 28, 29.
6. **Other settings** (T4): ramified rings of integers and cyclotomic norms (Thm 4.1, Ex. 4.1a–b); matrices,
   (I+Y)^N − Σ_{m<k} C(N,m)Y^m = C(N,k)Y^k U with U ∈ GL_d(Z_(p)) (Thm 4.2); Abel and q-deformations (Thm 4.3, 4.4);
   Bhargava binomials, conditional on Bhargava's integrality theorem (Rem 4.5).

How this relates to the starting identity: Ishan's AB identity is the x = −1 case of the classical column-sum identity
(T1 of the branch, Thm 1.2(b)), which is equivalent to the Chaundy–Bullard identity (Chaundy–Bullard 1960;
Koornwinder–Schlosser, Indag. Math. 19 (2008), (1.1)) through the classical truncation identity
Σ_{m≤j} C(N,m) w^m (1−w)^{j−m} = Σ_{k≤j} C(N−1−j+k, k) w^k at w = 2. The consolidation re-checked this equivalence
(`consolidation_2026-10-26/t16_checks.py`, part CB). Item 4 is the arithmetic of that identity; items 1–3 are its p-adic generalization.

## Classical inputs (from the branch's T1, which is not imported)
T2–T4 cite the branch's file T1. That file is classical material, and it also transcribes Ishan's private note, so it
is left out; the inputs it supplies are these.
- **Tail lemma** (T1 Lemma 1.1). For N ≥ 0, k ≥ 1: (1+y)^N − Σ_{m<k} C(N,m)y^m = y^k Σ_{r=k−1}^{N−1} C(r,k−1)(1+y)^{N−1−r}.
  Proof: call the left side F(N); F(0) = 0, and Pascal's rule gives F(N) = (1+y)F(N−1) + C(N−1,k−1)y^k; unroll.
- **Column-sum identity** (T1 Thm 1.2). For 0 ≤ j ≤ n:
  (1−x)^{j+1} Σ_{r=j}^{n} C(r,j)x^r = x^j − x^{n+1} Σ_{m=0}^{j} C(n+1,m) x^{j−m}(1−x)^m, equivalent (through the truncation
  identity above) to the Chaundy–Bullard identity (1−x)^{j+1}Σ_{k≤n−j}C(j+k,k)x^k + x^{n−j+1}Σ_{k≤j}C(n−j+k,k)(1−x)^k = 1
  (Koornwinder–Schlosser, arXiv:0712.2125, (1.1)). T3 uses it at x = −1, through the tail form
  A(n,j) = (−1)^n W_{n+1,j+1}(−2).
- **Reflection** (T1 Prop 1.6, used in T3 Prop 3.4). T_j(N) := Σ_{m≤j} C(N,m)(−2)^m satisfies T_j(j−N) = (−1)^j T_j(N),
  from the truncation identity at w = 2 and C(M+k,k) = (−1)^k C(−M−1,k).
- **q-forms** (T1 (1.8), (1.9), used in T4 Thm 4.4): Newton interpolation of X^N at the nodes 1, q, q², …, whose divided
  differences are Gaussian binomials, and Gauss's q-binomial theorem; both classical (Hoshi–Katori–Koornwinder–Schlosser,
  arXiv:2304.10003, (1.12), for the q-Chaundy–Bullard identity).

## Proof 1: subset-of-a-subset identity plus the ultrametric inequality, with Kummer's carries
- **Mechanism.** C(N,k+i)C(k+i,i) = C(N,k)C(N−k,i) gives W_{N,k}(y) = C(N,k)·Σ_i [C(N−k,i)/C(k+i,i)] y^i; the i-th
  term has valuation ≥ iν − εν_p(C(k+i,i)), which is ≥ 1 under (★). Kummer's theorem turns (★) into a carry condition
  (digit criterion), and for p = 2 the first failing index gives the converse.
- **Files:** `source/theorems/T2_higher_order_LTE.md`, `source/theorems/T3_AB_two_adic.md`,
  `source/theorems/T4_other_settings.md`; scripts and evidence in `source/scripts/`, `source/evidence/`.
- **Lean:** none.
- **Audit:** none recorded. See the consolidation checks below.
- **Provenance and novelty.** The branch's searches L1–L26 (`source/literature/SEARCH_LOG.md`) and the consolidation's
  searches B1–B2 (`consolidation_2026-10-26/literature_searches_2026-10-26.txt`) found no statement of items 1–5 for k ≥ 2; the closest relatives are
  P. Miska, arXiv:1803.04533 Thm 2.2 (an LTE-type statement for linear combinations of exponentials near a zero) and
  M. Billal, arXiv:1509.03288. Label: NEW (no match found), with the caveat that the proof of item 1 is three lines
  and the statement may be folklore. The matrix case k = 1 of item 6 is standard (congruence subgroups) and is not
  claimed. The classical identities used (Chaundy–Bullard, Abel, Gauss's q-binomial theorem, Newton interpolation)
  are not claimed.

## Consolidation checks (2026-10-26)
`consolidation_2026-10-26/t16_checks.py` (output `consolidation_2026-10-26/t16_checks.out`; new code, exact integers):
- CB: the equivalence through the classical forms at x = −1, all 0 ≤ j ≤ n < 160 (12,880 pairs), and symbolically in x
  for n ≤ 12;
- T21: item 1, part (c) and the refined part (b), 206,784 cases (p < 30, k ≤ 40, ν ≤ 3, three units, N ≤ 160), with (★)
  decided by brute force;
- C32: item 2, 52,878 cases (p < 30, p^s ≤ 64, N ≤ 400);
- T35: the digit criterion against brute force (p ∈ {2,3,5,7}, ν ≤ 3, k < 700 or 400), 5,688 cases;
- T36: the converse for p = 2 at N = k + i*, 1,280 cases; T37: the three examples of Thm 3.7 and Prop 3.8;
- T31: item 4 for all j < 300 and n < 700, and the list of admissible j; P33: item 5 for 5 ≤ n < 20000, with α to 80 bits;
- T42: the matrix theorem, 2,160 cases (3 × 3, p ∈ {2,3,5}, k = p^s).
A mutation test (`consolidation_2026-10-26/t16_checks_MUTANT.out`) shows the checks detect a wrong formula. The proofs of T2 Thm 2.1,
Cor 3.2, 3.3, 3.4, Thm 3.5, 3.6, 3.7, Prop 3.8, T3 Thm 3.1, Prop 3.3, 3.4 and T4 Thm 4.2 were read by the consolidation
(`consolidation_2026-10-26/proof_reading_notes_2026-10-26.txt`); no gap was found.

## Open questions
- An exact formula for ν_2(A(n,j)) in the non-admissible columns (governed by 2-adic roots of the cofactors).
- A closed description of the odd-p tie case (Thm 3.7 is an algorithm).
- The density of K_2(1) (about 0.351 up to 2^20, an observation).
- Stirling and Hurwitz analogues: no analogue of the subset-of-a-subset identity (`source/notes/OPEN_PROBLEMS.md`).
