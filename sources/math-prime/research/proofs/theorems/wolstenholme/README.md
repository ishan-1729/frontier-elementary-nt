# Wolstenholme's theorem: H_{p−1} ≡ 0 (mod p²), H^{(2)}_{p−1} ≡ 0 (mod p), C(2p−1, p−1) ≡ 1 (mod p³) for p ≥ 5
Status: WRITTEN (unreviewed; checked numerically by the 2026-10-26 consolidation)     Methods: not yet classified (finitary identity: interpolation / partial fractions)
Generalizes: none     Generalized by: none (the classical generalizations are listed under Provenance)

Origin: branch `origin/claude/math-prime-cross-domain-g1vwsn`, tip 26c718ac, folder
`cloud_agent_claude/17_finitary_identities_nt/` (cloud Claude programme, track 17, 18–19 Oct 2026), file
`theorems/T1_interpolation_identity_wolstenholme.md`, in [source/](source/). The track was stopped before it wrote a
README. No independent review is recorded.

## Statement
Let p be a prime and H^{(m)}_{n} = Σ_{k=1}^{n} k^{−m}, H_n = H^{(1)}_n.
- (Thm 1.4) For 1 ≤ m ≤ p−2: H^{(m)}_{p−1} ≡ 0 (mod p).
- (Thm 1.5) For p ≥ 5: (a) H_{p−1} ≡ 0 (mod p²); (b) H^{(2)}_{p−1} ≡ 0 (mod p); (c) C(2p−1, p−1) ≡ 1 (mod p³).
- (1.2) For every odd n ≥ 1, exactly: 1/C(2n−1, n−1) = 1 − 2n² Σ_{k=1}^{n−1} (−1)^k C(n−1,k)/(k² − n²).

Wolstenholme (1862) for (a)–(c); all KNOWN. The point is the proof.

## Proof 1: one interpolation identity, read through its Taylor coefficients and its values at x = ±p
- **Identity (I):** Σ_{k=0}^{N} (−1)^k C(N,k) x/(x+k) = 1/C(x+N, N) (Lagrange interpolation of 1 at 0, −1, …, −N).
- **Mechanism.** The Taylor coefficients at x = 0 give Σ_{k≥1}(−1)^{k−1}C(N,k)k^{−m} = h_m(1, 1/2, …, 1/N) (Dilcher's
  formula); at N = p−1 the weights (−1)^k C(p−1,k) = Π_{j≤k}(1 − p/j) are ≡ 1 (mod p), so h_m ≡ −H^{(m)}, and Newton's
  identity m h_m = Σ_i H^{(i)} h_{m−i} gives (m+1)H^{(m)} ≡ 0 by induction. Modulo p², the same expansion gives (a).
  The values at x = ±n give (1.2), hence (c).
- **What it avoids:** Fermat's little theorem, Wilson's theorem, primitive roots, the inverse map on (Z/p)^×, Bernoulli
  numbers.
- **Files:** `source/theorems/T1_interpolation_identity_wolstenholme.md`; `source/scripts/t1_verify.py`,
  `source/evidence/t1_verify.txt`; Lean draft `source/lean/T17_FermatFree.lean` (cases m = 1, 2 of Thm 1.4). **The Lean
  file has no recorded build: the consolidation did not compile it, and it must not be cited as LEAN.**
- **Audit:** none recorded. See the consolidation checks below.
- **Provenance and novelty.** No source deriving H^{(m)}_{p−1} ≡ 0 (mod p) from Dilcher's formula by Newton's identity,
  or deriving (1.2), was found (branch: `source/literature/SEARCH_LOG.md`; consolidation: searches W2–W3 in
  `consolidation_2026-10-26/literature_searches_2026-10-26.txt`). Closest
  relatives: MacMillan–Sondow, Amer. Math. Monthly 118 (2011) 549–551 (arXiv:1011.0076), an elementary proof of the
  power-sum congruence from Pascal's 1654 identity (positive powers, not reciprocals); K. Dilcher's partial-fraction
  method for q-harmonic sums (Electron. J. Combin. 15 (2008) R63). Label: NEW PROOF OF A KNOWN RESULT (route not found;
  every ingredient is classical). Not new and not claimed: C(cp−1, p−1) ≡ 1 (mod p³) for all p-integral c, and Morley's
  congruence as its case c = 1/2 (T1 Cor 1.6; Bencherif–Boumahdi, arXiv:1607.00700, (10) and Cor. 3).

## Consolidation checks (2026-10-26)
`consolidation_2026-10-26/t17_checks.py` (output `consolidation_2026-10-26/t17_checks.out`; new code): identity (I) and (I′) symbolically for N ≤ 14; Dilcher's
formula for N ≤ 25, m ≤ 6; the two congruences of the proof of Thm 1.4 and its conclusion for all odd p < 400 and
1 ≤ m ≤ p−2 (13,731 pairs, in F_p); (1.2) exactly for odd n < 160; (a)–(c) for 5 ≤ p < 1000; Cor 1.6 and Morley for
5 ≤ p < 120. The proofs of Lemmas 1.1–1.3 and Theorems 1.4–1.5 were read by the consolidation
(`consolidation_2026-10-26/proof_reading_notes_2026-10-26.txt`); no gap was found.

## Open questions
- Composite moduli: (1.2) gives only C(2p^r−1, p^r−1) ≡ 1 (mod p²) (`source/notes/OBSTACLES.md`, O7).
- Leudesdorf's theorem from one identity: blocked, because the nodes are no longer 1, …, p−1 (O1).
