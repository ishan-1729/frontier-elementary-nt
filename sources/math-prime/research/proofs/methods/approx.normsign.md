# approx.normsign: integer ladder with a Wronskian sign certificate
Parent: approx     Children: [approx.normsign.all_degree](approx.normsign.all_degree.md)

## What it does (Ljunggren, proof 2)
1. Pythagoras twice turns x²+1=2y⁴ into the quartic Thue equation Q(m,n)=1.
2. Re-centring at (3,2) and two symmetries reduce the problem to a sector where Φ(U,V)=28561 or −114244.
3. An explicit integer ladder T_n/S_n approaches the root r=tan(¼·arctan(1/239)). Formal power series prove
   the integrality, and a polynomial Wronskian sign lemma proves that the fractions alternate around r.
4. A determinant gap principle gives U≤6, and a congruence modulo 13 finishes.

## With a real parameter (quartic_t_family, proof 2)
The sign lemma already holds for every real a ≥ 1
([elementary_norm_sign_lemma.md](../theorems/ljunggren/norm_sign/elementary_norm_sign_lemma.md), eq. (2)). Branch A's
track 08 added three ingredients:
- the integrality S_n, T_n ∈ Z[a] of the ladder polynomials, by Lagrange inversion and the ODE
  (1 − 16az − 64z²)Y″ − (24a + 192z)Y′ − 60Y = 0;
- the explicit bounds e_n < (1/(4a))(4/a)ⁿ and S_n ≤ (16a + 8/a)ⁿ;
- centring at a known solution, which gives an irrational parameter (a = 27√6 for t=6, a = 102√3 for t=12) and
  linear forms with values in Z ∪ √6Z (resp. Z ∪ √3Z).

With them it proves the odd-Pell-index halves of 7x⁴−6y²=1 and 13x⁴−12y²=1. File:
[02_pade_norm_sign.md](../theorems/quartic_t_family/source/08_analysis_pade_thue/theorems/02_pade_norm_sign.md).

## Limits
- The method gives the irrationality exponent 1 + log(16a)/log(a/4), which is below Liouville's 4 only when a > 32.
- The even Pell index of t=6 needs an approximation to a different root, β₃ ≈ 0.2381, and no good base point is
  known. Among the base points with y₀ ≤ 2·10⁵, only 5/21 has E > 1 and κ < 4, and it gives κ = 3.381 where κ < 3
  is needed (a measurement; [quartic_t_family](../theorems/quartic_t_family/README.md), open question 1).
- t=2 (Bumby): κ ≈ 4.06 even for the odd half, so the method does not apply.

## Owner's note
Ishan, 2026-10-14: approximation proofs of this kind are classical in spirit (the Thue–Siegel and Padé
lineage). They are not the programme's distinctive method.

## Theorems using it
- [ljunggren](../theorems/ljunggren/README.md) (proof 2);
- [quartic_t_family](../theorems/quartic_t_family/README.md) (proof 2: the odd Pell index of t=6 and t=12).
- Through its child: [cubic_angle_sector](../theorems/cubic_angle_sector/README.md) (degree 3).

## Sources
- The paper: [papers/ljunggren.tex](../../shared/prior/elementary_diophantine_paper/papers/ljunggren.tex).
- The exposition: [norm_sign/](../theorems/ljunggren/norm_sign/elementary_proof.md).
- Lean: [CompleteNormSign.lean](../../../lean/Ljunggren/CompleteNormSign.lean).
- The parameter version: branch A, `origin/claude/math-prime-cross-domain-g1vwsn`, tip `26c718ac`,
  `cloud_agent_claude/08_analysis_pade_thue/`, copied to [quartic_t_family/source/](../theorems/quartic_t_family/source/).
