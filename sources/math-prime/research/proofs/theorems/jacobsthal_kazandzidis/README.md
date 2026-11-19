# Jacobsthal–Kazandzidis congruence: C(ap, bp) ≡ C(a, b) (mod p^{3+v}), v = ν_p(ab(a−b)), p ≥ 5
Status: WRITTEN (unreviewed; checked numerically by the 2026-10-26 consolidation)     Methods: not yet classified (symmetric product + p-adic logarithm)
Generalizes: Ljunggren's congruence (v = 0)     Generalized by: none

Origin: branch `origin/claude/math-prime-cross-domain-g1vwsn`, tip 26c718ac, folder
`cloud_agent_claude/17_finitary_identities_nt/`, file `theorems/T2_jacobsthal_kazandzidis.md`, in [../wolstenholme/source/](../wolstenholme/source/) (the track-17 bundle is imported once, under the Wolstenholme card, so that its internal links stay valid). No
independent review is recorded.

## Statement
Let p be an odd prime and a > b > 0 integers, v := ν_p(a) + ν_p(b) + ν_p(a−b).
- (Thm 2.1) C(ap, bp) ≡ C(a, b) (mod p^{3+v}) for p ≥ 5, and (mod 3^{2+v}) for p = 3. KNOWN (Jacobsthal 1952;
  Kazandzidis 1968; Ljunggren 1952 for v = 0).
- (Cor 2.7) C(ap, bp)/C(a, b) ≡ 1 + p H_{p−1}·ab(a−b) (mod p^{4+v}) for every odd p, and for p ≥ 5 one may write
  −(p²/2) H^{(2)}_{p−1} for p H_{p−1}. Novelty UNCLEAR: for v ≤ 1 it follows from Zhao's congruence modulo p⁵
  (C(np, mp)/C(n,m) ≡ 1 + w_p nm(n−m)p³, p ≥ 7, as reported in Meštrović's survey arXiv:1111.3057); for v ≥ 2 it was not
  found stated, but Helou–Terjanian (2008, modulo p⁶ with Bernoulli numbers) was not read.

## Proof 1: the symmetric product f(n) = C(np−1, p−1), f(1−n) = f(n), and odd summation polynomials
- **Mechanism.** C((b+c)p, bp)/C(b+c, b) = Φ(b+c)/(Φ(b)Φ(c)) with Φ(N) = Π_{n≤N} f(n). The reflection f(1−n) = f(n)
  (pairing r ↔ p−r) makes f a polynomial in n(n−1), whose coefficients have valuation ≥ 2j, and the linear one ≥ 3
  because e₂(1, 1/2, …, 1/(p−1)) ≡ 0 (mod p) (from the Wolstenholme card). Taking p-adic logarithms, the discrete
  antiderivatives P_j(N) = Σ_{n≤N}(n(n−1))^j are odd, so P_j(b+c) − P_j(b) − P_j(c) = bc(b+c)R_j(b,c), which carries the
  three factors of v.
- **Imports:** the p-adic logarithm on 1 + pZ_p (Robert, GTM 198, Ch. V §4; Koblitz, GTM 58, Ch. IV §1); the file notes
  that a truncated logarithm suffices. Theorem 1.4 of the Wolstenholme card (m = 1, 2, 3).
- **Files:** `../wolstenholme/source/theorems/T2_jacobsthal_kazandzidis.md`; `../wolstenholme/source/scripts/t2_verify.py`, `../wolstenholme/source/evidence/t2_verify.txt`.
- **Audit:** none recorded. See the consolidation checks below.
- **Provenance and novelty.** Label: NEW PROOF OF A KNOWN RESULT, but a variant of published routes, not a new
  mechanism. Robert–Zuber (Rend. Sem. Mat. Univ. Padova 94 (1995) 235–243) use the oddness of log Γ_p and Bernoulli
  estimates; Granville ("Arithmetic properties of binomial coefficients I", CMS Conf. Proc. 20 (1997), §11) derives
  Jacobsthal's congruence from a p-adic exponential/logarithm identity for products of integers prime to p, with
  Bernoulli-number input (consolidation reading of the HTML version, whose formulas are images; searches W1, W4–W6 in
  `../wolstenholme/consolidation_2026-10-26/literature_searches_2026-10-26.txt`). This proof replaces
  Γ_p and the Bernoulli estimates by the explicit symmetric polynomial f and the single input e₂ ≡ 0 (mod p). Robbins
  (1988) gave an elementary proof with a weaker exponent.

## Consolidation checks (2026-10-26)
`../wolstenholme/consolidation_2026-10-26/t17_checks.py` (output `../wolstenholme/consolidation_2026-10-26/t17_checks.out`; exact rationals): Thm 2.1 and Cor 2.7 (both forms) in 19,470
cases (odd p < 40, 1 ≤ b < a ≤ 60); the exponent of Thm 2.1 is attained for every prime tested. The proof (Lemmas
2.2–2.5, Steps 1–3, the proof of Cor 2.7) was read by the consolidation; no gap was found.

## Open questions
- p = 2: the reflection fails (f(1−n) = −f(n)); not attempted (`../wolstenholme/source/notes/OBSTACLES.md`, O4).
- The next order (mod p^{5+v}) needs e₂ (mod p³), e₃ (mod p²), e₄ (mod p) (O5).
