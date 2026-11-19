# Legendre-symbol Hankel determinants: eight conjectures of Z.-W. Sun (2024) and track 03's Conjecture 4.7
Status: WRITTEN (unreviewed; checked numerically by the 2026-10-26 consolidation)     Methods: not yet classified (proposed: a group-ring "inverse package" method, see Proof 1)
Generalizes: none     Generalized by: none

Origin: branch `origin/claude/math-prime-cross-domain-g1vwsn`, tip 26c718ac, folder
`cloud_agent_claude/12_legendre_determinants/` (cloud Claude programme, track 12, 18–19 Oct 2026). The branch's files
are in [source/](source/). No independent review is recorded: the two referees started on 19 Oct 2026 were stopped
before writing a report; their partial scripts and outputs are in `source/incomplete_review_evidence/`.

## Statement
Notation. p is a prime, n = (p−1)/2, χ = (·/p) with χ(0) = 0, E = χ(−1), p* = Ep, L = (2/p).
- For p ≡ 1 (mod 4): ε_p > 1 is the fundamental unit and h(p) the class number of Q(√p), and
  ε_p^{h(p)} = a_p + b_p√p with 2a_p, 2b_p ∈ Z (Sun's (1.1)).
- For p ≡ 3 (mod 4): h(−p) is the class number of Q(√−p), (h) := (−1)^{(h(−p)−1)/2}, and ⟨n!⟩ ∈ {±1} is the sign
  with n! ≡ ⟨n!⟩ (mod p).
- D^{(c)}_{[s,t]}(x,y,z,w) := det[ x + ((j+k+c)/p) + (j/p)y + (k/p)z + (jk/p)w ]_{s≤j,k≤t}.
- v₀ := wx − (y+1)(z+1), c_p := (p+1)b_p − 2, d_p := (p+1)b_p − 2(a_p+1).

**(A) Sun's conjectures** (Z.-W. Sun, *Problems and results on determinants involving Legendre symbols*,
arXiv:2405.03626v8, 28 Jul 2024; still the latest version on 26 Oct 2026, no journal reference listed). For every prime
p > 3 in the stated class:

| Conj. | determinant | p ≡ 1 (mod 4) | p ≡ 3 (mod 4) |
|---|---|---|---|
| 3.2(i) | D^{(0)}_{[0,n−1]} | L·2^{n−1}[(p b_p − 2a_p)x + (a_p − 2b_p)v₀] | 2^{n−1}(v₀ − 2x) |
| 3.5(i) | D^{(0)}_{[1,n]} | L·2^{n}[a_p(w−x) + b_p + (b_p−1)(y+z) − c_p(wx−yz)] | −2^{n}[w + x + (h)(y+z+2yz−2wx)] |
| 3.6(i) | D^{(0)}_{[1,n−1]} | L·2^{n−1}[b_p − a_p x + (a_p−2b_p)w + (b_p−1)(y+z) + d_p(yz−wx)] | 2^{n−1}[w + x + 2(wx−yz) + (h)(y+z+2yz−2wx)] |
| 4.2 | D^{(−1)}_{[1,n−1]} | L·2^{n−2}[((2−L)a_p − p b_p)x + (a_p + (L−2)b_p)(w+y+z+1) + ((p−1)b_p + (L−1)(a_p+b_p))(yz−wx)] | 2^{n−2}[(L−2)x − w − y − z − 1 + (1−L)(yz−wx)] |
| 4.3 | D^{(−1)}_{[0,n]} | L·2^{n−1}[(p b_p − 2a_p)(w+y+z+1) + (2b_p − a_p)p x + ((p−2)a_p − p b_p)(yz−wx)] | 2^{n−1}[2w(1−x) + 2(y+1)(z+1) + p((w+1)x − yz)] |
| 4.4 | D^{(−1)}_{[0,n−1]} | (not conjectured; see (C)) | 2^{n−2}[((p/2)(2−L) − 4)x + (p/2 − 2 − L)(w+y+z+1) + ((p/2)(L−1) + 2 − L)(yz−wx)] |
| 4.5 | D^{(1)}_{[0,n]} | L·2^{n}[p b_p(x + ((p−2)/2)(yz−wx)) + a_p(w(x + (p−2)/2) − (y+1)(z+1))] | 2^{n}[w((p−2)/2 − x) + (y+1)(z+1)] |
| 4.6 | D^{(1)}_{[1,n]} | with w = 0: L·2^{n−1}[(p b_p − 2a_p)x + 2(n_p + b_p − a_p)yz + (2b_p − a_p − 1)(y+z+1) + 1] "for some positive integer n_p" | 2^{n−1}[(1 − (h))(y + z + 2(yz−wx)) + (p−3)(yz − wx + w/2) − 2x + 1] |

The proof gives n_p = (p(p−1)b_p − 2(p−3)a_p − 4)/4 and shows that it is a positive integer (Sun's n_5 = 1,
n_13 = 11, …, n_53 = 256 are reproduced).

**(B) Track 03's Conjecture 4.7** (internal to the cloud Claude programme, not published). With
a_k = (χ(k) − χ(k+1))/2 and H_r = (a_{i+j+1})_{0≤i,j≤r−1}:
det H_{n−2} = s_p (p − 4(−3/p) + 3(−1/p))/12, s_p = −(2/p) for p ≡ 1 (mod 4) and s_p = (−1)^{(h(−p)+1)/2} for
p ≡ 3 (mod 4). Equivalently |det H_{n−2}| = ⌊p/12⌋ + [p ≡ 2 (mod 3)].

**(C) Further evaluations proved on the way** (all elementary, in the files named):
- T1 Thm 1.1: det H_n = (2/p) (p ≡ 1), = ⟨n!⟩ (p ≡ 3); det H_{n−1} = −E det H_n; det H_{n−2} = −((p* + 3 − 4χ(3))/12) det H_n.
- T1 Thm 1.3: det H_{n−k} = (−1)^k Δ_k det H_n for 0 ≤ k ≤ n, Δ_k a k × k determinant of K-coordinates of ζ^{−b};
  Cor 1.4 gives Δ_3 as an explicit cubic in p*.
- T2 Thm 2.2, Prop 2.3: with A = (χ(j+k+1))_{0≤j,k≤n} and Π_R := Π_{r∈R}(1 − ζ^r) = u + v g (g the quadratic Gauss sum),
  𝟙ᵀA⁻¹ψ = 1, ψᵀA⁻¹ψ = −E(p−2)/2, 𝟙ᵀA⁻¹𝟙 = Eu/v, and det A = 2^n v det H_n.
- T3 Prop 3.3: the 15 entries of A⁻¹ on the vectors 𝟙, ψ, e₀, e_{n−1}, e_n; T3 Thm 3.2: the four-parameter determinant of
  Conj. 4.6 with w, and D^{(−1)}_{[0,n−1]} for p ≡ 1 (mod 4).

**Imported theorems** (named hypotheses; nothing else is imported in the main line):
- (M) Mordell, Amer. Math. Monthly 68 (1961) 145–146: n! ≡ (−1)^{(h(−p)+1)/2} (mod p) for p ≡ 3 (mod 4), p > 3. Used only
  to write the signs of (B) for p ≡ 3, of Sun's 3.5(i), 3.6(i), 4.6 for p ≡ 3, and T1 Cor 1.5(ii), through h(−p).
- (D) Dirichlet's class number formula in product form, Π_{k=1}^{n}(1 − ζ^{k²}) = √p ε_p^{−h(p)} (p ≡ 1 mod 4; sources:
  Sun, Finite Fields Appl. 59 (2019), Thm 1.3(i), (1.12); Chapman, Acta Arith. 115 (2004) p. 237), **and**
  (G) Gauss's sign g = +√p (p ≡ 1 mod 4; Borevich–Shafarevich Ch. 5 §4 Thm 7). Used only in T2 Lemma 2.7, to identify
  (−v, u/p) with (a_p, b_p). Note: the p ≡ 1 case needs both (D) and (G), not the class number formula alone.
- Without any import, every identity holds in the class-number-free form: (h) → −⟨n!⟩ and (a_p, b_p) → (−v, u/p)
  (T3 Thm 3.1′).

So: for p ≡ 3 (mod 4), Sun's 3.2(i), 4.2, 4.3, 4.4 and 4.5 are proved with no import, and 3.5(i), 3.6(i), 4.6 as stated
under (M); for p ≡ 1 (mod 4), 3.2(i), 3.5(i), 3.6(i), 4.2, 4.3, 4.5, 4.6 as stated are proved under (D) and (G).

## Proof 1: group ring Q[Z/p] ≅ Q × Q(ζ), and the inverse package of Chapman's matrix
- **Mechanism.** A four-parameter window of the Hankel matrix (χ(j+k+c)) is a bordered determinant. By the matrix
  determinant lemma and a window formula (T3 Lemma 4.1) it is governed by a few entries of A⁻¹, A = C*_p(0). These are
  found by solving A v = t in Q[Z/p]: a column of A is the restriction of x^{−k−1}γ (γ = Σχ(t)x^t), so the system splits
  into one scalar equation (trivial character) and one equation over K = Q(√p*), written in the K-basis {ζ^w : w ∈ W},
  W = {n+1, …, 2n}. The only global quantities are det A and Π_R = u + vg. Cramer's rule ties det A to the unimodular
  Hankel matrix H_n (det A = 2^n v det H_n); the sign of det H_n comes from Euler's criterion and a Vandermonde
  factorisation modulo p; for p ≡ 3 (mod 4) the sign of v comes from a congruence modulo (1 − ζ)^{n+1}.
- **Files:** `source/theorems/T1_hankel_differences.md`, `source/theorems/T2_sun_conjecture_4_5.md`,
  `source/theorems/T3_hankel_windows.md`; the branch README and scripts s01–s11 with evidence E01–E11 in `source/`.
- **Lean:** none.
- **Audit:** none recorded. See the consolidation checks below.
- **Provenance and novelty** (consolidation searches S1–S17 of 26 Oct 2026,
  [consolidation_2026-10-26/literature_searches_2026-10-26.txt](consolidation_2026-10-26/literature_searches_2026-10-26.txt);
  the branch's own searches Q1–Q29 are in `source/literature/LIT_LOG.md`):
  - No proof of Sun's 3.2(i), 3.5(i), 3.6(i), 4.2, 4.3, 4.4, 4.5 or 4.6 was found before 26 Oct 2026. The papers that
    settle conjectures of the survey settle other items: Chen–Sun arXiv:2409.08213 (3.1, 3.7(i), 4.1, 3.9 for p ≡ 1),
    Liu–Sun–Wang 2408.07034 (3.2(ii)), Yang–Yang–Zhang 2606.22548 (3.4), Yang–Zhang 2605.19517 (4.8(i), p ≡ 1).
    Yang–Zhang 2605.19502 and Li–Yan 2608.08509 use the numbers "4.6", "4.8(i)" for conjectures of other Sun papers.
  - So the results are NEW as proofs of published open conjectures (Ishan's criterion: no earlier proof of these
    theorems exists). The method is a relative of Chapman's projection method (2004, Lemma 1 and Thm 2) and of the
    "inverse package" reduction of Chen–Sun (Thm 3.1) and Yang–Zhang (Prop. 2.1); the group-ring evaluation of the
    Hankel inverse package (T2 Thm 2.2, T3 Prop 3.3) was not found in print.
  - Known parts: T1 Thm 1.1(a) is equivalent to Chapman 2004 Thm 5 with (M); T1 Cor 1.5(ii) is Chapman's Thm 5;
    T2 Lemma 2.6 (Π_R = ⟨n!⟩g for p ≡ 3) is Chapman 2004 Lemma 4 and Sun 2019 (1.13), and its step (b) is the classical
    Stickelberger congruence for the quadratic Gauss sum, so its proof is not claimed as new; T3 Prop 3.6 is Chapman's
    Cor. 3; T3 Remark 6.3 (Sun's 3.3 for p ≡ 3) is a direct corollary of Yang–Yang–Zhang's preprint and is not claimed.
- **Dependency:** T1 uses track 03's notation, its Conjecture 4.7 and its Basis Lemma (T4 of track 03), but re-proves the
  Basis Lemma (T1 Appendix A), so the proof is self-contained. Track 03 itself is triaged separately.

## Consolidation checks (2026-10-26)
New code in [consolidation_2026-10-26/](consolidation_2026-10-26/), written without the branch's scripts (Python and numpy
only; exact integers by modular elimination with Hadamard's bound and CRT; class numbers by reduced forms and cycles,
units by continued fractions, Π_R exactly in Z[ζ]; library `t12lib.py`, tested by `t12lib_selftest.py`):
- `t12_sun_direct.py` (`sun_direct_g1.out` … `sun_direct_g8.out`; the prime groups are in `groups8.txt`): all eight
  conjectures, Thm 3.2(a), (b), n_p (closed form, positivity, Sun's values) and the controls 3.1(i), 4.1(i), 4.1(ii),
  as 16 exact integer determinants each at {0,1}⁴, for all 76 primes 5 ≤ p ≤ 397, in both the stated and the
  class-number-free forms. No failure.
- `t12_T1.py` (`t1_5_400.out`): T1 Thm 1.1, Cor 1.2 in both sign forms (track 03's Conj. 4.7), Mordell's congruence,
  |det H_{n−2}|, Cor 1.4, Cor 1.5(i), (ii), Remark 7.1, all primes 5 ≤ p ≤ 397. No failure.
- `t12_T1_allk.py` (`t1_allk_5_160.out`): T1 Thm 1.3 for every k ≤ n, with the coordinates of ζ^{−b} obtained by exact
  linear algebra in Q(ζ), not by the branch's Newton-identity formulas; all primes 5 ≤ p ≤ 157. No failure.
- `t12_inverse.py` (`inverse_5_330.out`, `inverse_331_400.out`): T2 Thm 2.2, Prop 2.3, Lemmas 2.5–2.7, T3 Prop 3.3 (all
  15 entries, exact rational A⁻¹), Prop 3.6, Lemma 4.3, integrality of n_p, Remark 6.3, and the imports (D), (G)
  themselves (u = p b_p, v = −a_p exactly; Π_R = √p ε^{−h} at 80 digits), all primes 5 ≤ p ≤ 397. No failure.
- `t12_symbolic_assembly.py` (`symbolic_assembly.out`): an independent sympy re-derivation of the only computer-algebra
  step of the proof (T3 §5), from the Γ-table of Prop 3.3 and Lemmas 4.1–4.3; all 22 (conjecture, residue class)
  identities hold exactly, and n_p comes out in closed form. The script refuses floating-point coefficients. Mutation
  test `t12_symbolic_assembly_numeric_mutation.py` (`symbolic_assembly_numeric_mutation.out`): the same functions at the
  exact data of p = 13 and p = 11 give 0 in all 22 cases, and a perturbed Γ-entry is caught exactly in the frames that
  use it.
- Proof reading: every step of T1–T3 was re-derived by hand (`proof_reading_notes_2026-10-26.txt`); no gap was found.
  This is one reader's check, not a recorded independent audit, so the status stays WRITTEN.
- Statements: three independent readings of v8 agree term by term (the branch, from the rendered pages; the stopped
  referee, from the TeX source, `source/incomplete_review_evidence/evidence_12_T2T3/r01_direct_dets.gp`; the
  consolidation, from the arXiv HTML of v8). Searches: `literature_searches_2026-10-26.txt`, S1–S17.

## Open questions
- Sun's Toeplitz-type conjectures (3.3 for p ≡ 1, 3.5(ii), 3.6(ii), 4.7, 4.8 for p ≡ 3, 4.9): the group-ring system does
  not split there (exact obstacle in `source/notes/LEADS.md`, L1).
- A closed form of Δ_k for all k (L3).
- Observation L4, not claimed: the g-coordinate v of Π_R satisfies v ≡ ((p−1)/2)! (mod p) for every prime p ≥ 5 (checked
  for p < 400 by the consolidation); probably classical (a congruence of Chowla type), source not identified.
