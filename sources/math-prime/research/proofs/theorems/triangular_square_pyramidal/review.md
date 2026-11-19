# Review record: triangular_square_pyramidal

All reviewers below are Claude subagents of the same model family as the author. So this is in-session review, not
an independent audit. After round 2 the card was WRITTEN (Theorem A) and PARTIAL (Theorem B, under hypothesis (H)).
Round 3 reviewed the proof of (H) in [hypothesis_H.md](hypothesis_H.md); since then both theorems are WRITTEN
(in-session review).

## Round 1: scripts and independent verification (2026-10-15)
Two workflows ran, with ten agents (five per curve; this record covers the 8604.a1 agents, and the literature
agent). Each worked alone, wrote its own script, and reported in a structured form. The author wrote
[core.py](scripts/core.py), [01_coverage.py](scripts/01_coverage.py) and [02_exclusion.py](scripts/02_exclusion.py)
from the template's, with new data; [lifting2.py](scripts/lifting2.py) is a verbatim copy.

| Agent | Task | Verdict | What it found |
|---|---|---|---|
| fields | wrote [00_fields.py](scripts/00_fields.py) | 172 checks pass | Fields, maximal orders (also by sympy's round_two), markers, the substitution lemma for class (1,0), units, known solutions; the three fields are pairwise non-isomorphic, though (0,1) and (1,1) have the same discriminant. |
| reconstruction | wrote [03_reconstruction.py](scripts/03_reconstruction.py) | 1075 checks pass | The 2-descent re-derived; u² ≡ F_cl (mod 12) in all three classes, so g = √c is automatic; the trivial class gives x ∈ {−5, 7, 1099}; 24 points; T(n) = P(m) values 0, 1, 55, 91, 208335; Uchiyama's model 3Y² = X³ − X + 3 has 18 integral points. O_K = Z[λ], λ³ − λ + 3 = 0, is Uchiyama's cubic ring. |
| coverage | wrote [verify_coverage.py](scripts/verify_coverage.py), no project imports | 230 checks pass | Own Sturm and Rouché enclosures, own polygon by half-plane clipping. The enumeration runs to the weaker linear bound of Remark 3.2(ii) in every class (up to 672 million lattice points), so it does not use the log-concavity step; one orbit per (class, c), as in 01; the unit pairs are fundamental. |
| exclusion | wrote [verify_exclusion.py](scripts/verify_exclusion.py), no project imports | 1465 checks pass | Every number printed by 02 reproduced. Λ_q by its own method, certified by index. Prop. 17 in four bases at all 24 (q, zero) pairs: λ and ν(J) identical at 22 of them, all 11 used pairs among them; the two that differ fail in every basis. Negative controls bite. |
| brute force and audit | wrote [verify_bruteforce.py](scripts/verify_bruteforce.py) | all checks pass | Exact box \|a\|,\|b\| ≤ 500 and modular box ≤ 3000: only the certified zeros. Direct search \|m\| ≤ 10⁶ with rigorous windows: only the known solutions of G and of F_cl. Line-by-line audit of 02 and lifting2. |
| literature | wrote [literature_check.md](literature_check.md) and [lit_uchiyama.py](scripts/lit_uchiyama.py) | the script's checks pass | Read Uchiyama (1984) in full from the scan: the same fields, equations and units as R1, and a congruence-only exclusion that is incomplete in print (index-2 "integral bases", non-units in Props. 3, 7 and 8, a false congruence step); 117 of 155 tested statements are true. Found Finkelstein–London (1972), not cited before in the project. The author edited the note: its status, and a cross-reference to a sibling session's uncommitted file, replaced by a direct citation. |

**Fixes made after round 1.**
1. *Printing.* [01_coverage.py](scripts/01_coverage.py) printed three positive vertex coordinates of class (1,1) as
   0.000000 (their values are 4.0·10⁻⁷, 5.8·10⁻⁸ and 5.7·10⁻⁸). It now prints them in exponent form, with the least
   coordinate. Positivity was always checked exactly.
2. *Known solutions located.* 02 skipped a known solution that it could not locate in an orbit, silently. It now
   asserts, per class, that every known solution lies in exactly one orbit. (The brute-force agent had closed the gap
   from outside: all 11 known pairs lie in exactly one orbit.)
3. *A checkable certificate.* 02 now prints, for each modulus used, the basis of Λ_q; after each step, the Hermite
   form of the current lattice; each certified class with its lattice; and each sibling class with its witness. The
   one sibling class is (88, 44) mod (264, 176, 88), in class (1,1), c = 4, emptied by Proposition 17 at (23, (1,4)).
4. *Wording.* 02's line "each pair ±z meets the family at most once" now cites the injectivity from the cone
   certificate as well. 01's docstring states the lemma for c ≠ 0.
5. *proof.md.* The interval widths are < 10⁻⁵² (not the template's 10⁻⁵⁵); the margins of Proposition 17 are stated
   exactly (eight instances at δ = e − 1 − ι, three with slack); the basis remark is sharpened to "a positive verdict
   never depends on the basis"; the D-matrix of class (0,1) is not diagonal.
6. *Re-runs.* After the changes, 01 and 02 were re-run (only the printed format changed, and the new lines), and the
   three verifiers were re-run against the new outputs. Their parsers of 02's output were adapted to the new line
   format.

**Facts recorded by the verifiers.**
- *Candidates that persist modulo powers of 2.* In the orbit (0,1), c = 36, a second cluster of candidate classes,
  at exponents ≡ (8144558, 12701350) mod 2²⁴ with no small rational reconstruction, persists modulo 2^k for every
  k ≤ 48. So no sieve with moduli 2^k, k ≤ 48, can finish this orbit; the odd moduli 43 and 13 do. Round 1 called this
  a 2-adic zero; that is presumed, not proved (round 2).
- *Moduli.* In orbit (1,1), c = 4, the moduli 4 and 8 add nothing; [7, 23] alone certifies it.
- *Λ_q is not h₁Z × h₂Z* at 17 of the 22 (orbit, modulus) pairs used (index 4 to 22; corrected in round 2); 02 uses
  Λ_q.

**Process notes.** Several agents ran one read-only `git log`, which their briefs forbade; nothing changed. The author
converted some files from CRLF to LF line endings while the agents worked; the agents checked that the content was
unchanged.

All scripts whose inputs changed were re-run after the fixes, and all pass (outputs in [scripts/](scripts/)).

## Round 2: referees on the text (2026-10-15)
Three referees worked alone on this card, each with one lens: coverage (§§2–3 with core.py, 00 and 01), exclusion
(§4 with lifting2.py and 02, re-run) and global (§§0, 1, 5–7, the README, the literature check, this record and the
shared files changed for the card). mordell_k_minus17 had three referees with the same briefs; the six ran in two
workflows. Each re-ran the scripts it needed and wrote its own scratch checks, among them:
- a 200-digit recomputation of table 3.3, and an own enumeration of the lattice points;
- a re-implementation of §4 from the tables of proof.md alone (784 checks), with random GL₂(Z) changes of basis, the
  minimal certifying lists of moduli, and a box scan of the exponents;
- an own sympy re-derivation of §1 and §5;
- a 2-adic lifting of the orbit (0,1), c = 36 to 2¹²⁰, with a Hensel certificate.

**Result.** 21 findings: 0 critical, 1 major, 20 minor. All three verdicts were "accept with minor revisions". No
finding changes Theorem A′, A or B or any certificate, and all were fixed.

**The major finding: the 2-adic zero.**
- Round 1 stated that the markers of the orbit (0,1), c = 36 have a 2-adic zero that is not an integer zero, so that
  no sieve by powers of 2 can finish that orbit. Only the persistence of candidates up to 2⁴⁸ had been checked
  (exclusion referee). The text was first cut back to that finite fact.
- The global referee then proved the claim: a Hensel certificate after lifting to 2¹²⁰, with ν(F) ≥ 120 > 14 =
  2ν(det DF). The author read its script line by line, added asserts and a docstring, and re-ran it as
  [check_2adic_zero.py](scripts/check_2adic_zero.py). The statement of §4.5 is restored with that proof. The
  multivariable Hensel lemma is p-adic analysis, and only this remark uses it.
- verify_exclusion.py gained two checks (2^kZ × 2^{k−1}Z ⊆ Λ_{2^k}, and the separation from the known zero mod Λ₄)
  and a corrected message; it now runs 1467 checks.
- A consequence, found by the global referee: Uchiyama's Prop. 8 is this orbit, so it cannot be completed by
  congruences modulo powers of 2, although his paper says it goes "just as in the proof of Proposition 2". This is now
  in proof.md §6.3, the literature check §1.7, the README and the method card ghost.unit.multi_unit, whose Limits had
  also said that no p-adic zero occurred on 8604.a1.

**Minor findings, by group (all fixed).**
- *Notation in §§2–4.* The unit group of Lemma 3.2 is now Γ, so it no longer clashes with the field form G. The
  coefficients a₀..a₃ are named in 2.1. The bound T₁ of Remark 3.2(iii) is defined, where T(1) had clashed with the
  triangular numbers. §4 uses r for the marker index, (s, s') for coset coordinates and q₁, …, q_s in 4.4.
- *Condensed steps restored.* Remark 3.2(i) gives the reason that ℓ(u) lies in H and ends "against step 4". The
  logarithm-free route of Remark 3.2(ii) is written out: segments, then vertex to edge point. §0.3 lists the
  convex-hull fact, and says that e^y ≥ 1 + y serves only Remark 3.2(i).
- *Printed evidence.* 01 now prints the Sturm counts and root-enclosure widths, the LDLᵀ pivots and the number of
  certified Fincke–Pohst ranges (the same script serves both cards). The bullets of §4.5 cite verify_exclusion.out.txt
  by section. "Ĝ-value" is now the count of points with xᵀ(Ĝ − 4ρI)x ≤ T.
- *Rounding.* The E bounds are rounded up throughout (101.6, not a mix); the README said 5.3√c, below the certified
  5.3072√c.
- *The basis remark.* J changes by a sign modulo p^{e+min(λ,e)}, not p^{e+λ}. Since ν(J) ≥ 2λ, the conclusion
  stands. "ν(J) ≤ 1 < e" was false at the instances with e = 1, and is restated.
- *Lemma 4.4 and the steps of 02.* A paragraph now explains why classes removed modulo the step lattices suffice, and
  "the only sibling class" refers to 02's steps.
- *The phantom bullet* is restated exactly: (38,2) mod Λ₁₃, and 231 subclasses, of which 230 fail mod 43.
- *Uchiyama.* h_K = 1 is stated in his paper without a source; only the unit comes from the Delone–Faddeev table. In
  Prop. 7 the matrix M is right and the units were transported wrongly. Prop. 1's congruence is false in most tested
  cases, and its conclusion fails near the solutions.
- *Conditional statements.* The 18 points of Uchiyama's model hold under (H). "The theorem is right" now says what
  proves it: Theorem A with Uchiyama's §1, given h_K = 1 and the unit.
- *Naming the hypothesis.* CATALOGUE, GENERALIZATIONS (§3.2 item 6, which also gained the 2-adic obstruction) and
  03's output now name the class list (H), implied by rank ≤ 2.
- *Literature check.* Thomas's Monthly problem is tied to Finkelstein–London only, and LMFDB's own account of its
  integral-point computation (lmfdb_data.md) is cited.
- *The README's sibling rule* now says "for the same q".
- *This record.* The round-1 table gained the literature agent's row. The Λ_q index range is corrected: the referee
  counted 16 of 22 pairs, and the output gives 17. So is the four-basis statement, and the re-run sentence.
- *proof.md §7* lists lit_uchiyama.py and check_2adic_zero.py; §6.1 says what the in-session review was.

**Process notes.**
- The referee workflows were cut off once, when the author's context ended. They were resumed: the coverage referee
  of this card returned from cache, and the other five started again.
- The author applied fixes while referees were still reading. They noted the file times and reviewed the latest
  versions.
- Referees fetched abstracts, reviews and page scans into their scratch folders; none of it is stored in the
  repository.

After the round-2 fixes every script of the card was re-run, and all pass.

## Round 3: the proof of (H) (2026-10-15 to 2026-10-26)
After round 2 the author wrote [hypothesis_H.md](hypothesis_H.md), with the scripts 04–07 and cover_search, and changed
proof.md (§§0–1 with Lemma 1.3′, §§5–7) and the README, so that Theorem B no longer depends on the rank. Round 3
reviewed this. One workflow ran nine agents, each alone. For each card there were three verifiers, which wrote their own
scripts and imported nothing from the card, and a referee for the whole of hypothesis_H.md, the proof.md changes, the
README and the shared files; one literature agent served both cards. The referees re-ran the card's scripts in scratch
copies and wrote their own checks.

| Agent | Task | Verdict | What it found |
|---|---|---|---|
| units | wrote [verify_H_units.py](scripts/verify_H_units.py) (Sections 1 and 4) | 110 checks pass | The norm three ways (own determinant, resultant, σ₁Q); the descent identities as identities in K[a,b,c]; the five triples with T₂ ≤ C from the exact bounding box of the ellipsoid; (4.1) a second time from a fundamental domain e < \|σ₁\| ≤ 1, without the descent lemma; a search \|a\|, \|b\|, \|c\| ≤ 150 finds exactly ±ε^k, −3 ≤ k ≤ 5; eight planted errors caught. |
| cover | wrote [verify_H_cover.py](scripts/verify_H_cover.py) (Section 2) | 58 checks pass | The cover by measure and by exact recursive splitting; the bound by its own Bernstein expansion and by interval arithmetic with bisection; the maximum 32733/32768 at (−3/8, −7/32, 7/16) and the five other boxes at that corner; the reduction to D on all 6072 points of (1/m)Z³/Z³, 2 ≤ m ≤ 12; Remark 2.5 by a fundamental domain for ε (24 elements of C); eight planted errors caught. |
| classes | wrote [verify_H_classes.py](scripts/verify_H_classes.py) (Section 5) | 109 checks pass | Section 5 without floating point: the identities in X; non-association twice (quotients, and ring maps to F₃ and F₂₃₉); all 64 parity vectors, of which 8 survive; the step modulo 4 for every μ ∈ R; its own characters at 71, 37 and 59; all 12 X with \|X\| ≤ 10⁶ in their classes, with explicit S; (5.4)–(5.6) replayed on each point. |
| referee | the whole card | accept with minor revisions | 18 findings: 1 major, 17 minor. Its own scripts for the box table, for the units and the fixed classes by a fundamental domain, for Sections 5–6, and for the integral points with \|X\| ≤ 2·10⁶. It read arXiv:1202.6019v1: the row "−239 E 8/9" is in the §7 table, p. 16. |
| literature | both cards | — | [literature_check.md](literature_check.md) §6: Taylor (1975/76) first proved Z[λ] norm-Euclidean, on Lemmermeyer's attribution; Cavallar–Lemmermeyer and Lemmermeyer's 1989 table record M = 8/9; LMFDB 3.1.239.1 agrees with ε; no ideal-free or table-free proof of (H) was found. Its identification scripts are [lit_H_lemmermeyer.py](scripts/lit_H_lemmermeyer.py) and [lit_H_lmfdb.py](scripts/lit_H_lmfdb.py). |

**Result.** No finding changes Theorem H, Theorem B or any certificate. The major finding concerned the record, not the
mathematics: proof.md §6.1, the README, the CATALOGUE row and GENERALIZATIONS cited an in-session review of
hypothesis_H.md that was not yet recorded here, a verifier sat in scripts/ unlisted, and the header of this file still
called Theorem B PARTIAL. This section is that record. All findings were fixed.

**Fixes made after round 3.**
1. *The record.* This section and the header; the verifiers listed in hypothesis_H.md (Checks) and in proof.md §7.
2. *The constants of (4.3).* The printed bounds ρ_hi ≤ 14.4949 and ρ_lo ≥ 0.026894 were true but too coarse to give
   q ∈ [A, B]. The text now prints E < 8.137981, ρ_hi < 14.49483, ρ_lo > 0.02689446 and A³·14.49483 ≤ 1 ≤ B³·0.02689446,
   and 06 checks exactly these numbers. C is defined as max(h(A), h(B)) = 6.765141… ≤ 6.76515, the bound that (4.4)
   enumerates.
3. *No calculus where none is needed.* (1.6) proves that the real root is unique from f(−3/2) = 9/8 and the
   factorization of f, not from critical values; (4.3) bounds h by two endpoint identities, not by convexity. 04 and 06
   check both.
4. *K₃ named.* (2.5) and a paragraph after (4.3′) define the exact threshold K₃ = 4(E² + E + 1)³/(E²(E + 1)²), where
   the text had used 05's names μ and λ (λ also names the generator). 05 now writes e for σ₁(ε) and ρ_lo, ρ_hi as the
   text does, checks ψ(ρ_lo) = ψ(ρ_hi) = K₃ symbolically, and prints that the rounded C would also fit here.
5. *Finiteness, not a limit.* (2.5) and 05 take the element of least T₂ from the finite set
   {ε^k y : T₂(ε^k y) ≤ T₂(y)}; the text had said that T₂ tends to infinity along the orbit.
6. *T_low.* (2.5) proves T_low ≤ T₂ and that T_low is positive definite; (4.4) says that 06 uses ρ = 2⁻³⁰ and lists
   the five triples.
7. *Condensed steps written out*, as on the 1611.a1 card after its round 4:
   - (5.4): why no element of 𝒮 divides the final quotients, the definitions of A and B, and the role of (5.3);
   - (5.5): Y ≠ 0 and N(R₀) ≠ 0; (5.6): why S ∈ R;
   - (3.3): why two prime elements are associates or coprime; (3.6): m′ ∈ R, and why ξ² ∈ R gives ξ ∈ R;
   - the preamble of §3 lists (1.3b) and the characteristic polynomial.
8. *Characters.* (5.8) defines χ as a Legendre symbol after a ring map R → F_p, and restricts the constancy argument
   to nonzero values: 49 = 7² is in the class of 1, and the character at (7, 4) is 0 on it. The letter η, which the
   other card uses for its unit, is no longer used there.
9. *The step modulo 4.* (5.7) prints the eight values of μ₀S² modulo 4 for the coset εT, writes Z where T had clashed,
   and records what the step forces on X modulo 4 (07 checks it on all 12 X). (5.9) calls the step a consequence of the
   condition at 2, read modulo 4: modulo 4 it does not separate 1 from π₃ₐπ₃꜀, nor π₃ₐπ₃ᵦ from π₃ᵦπ₃꜀.
10. *Non-association by hand.* (1.5) separates π₃ₐ, π₃ᵦ, π₃꜀ and π₂₃₉ by ring maps R → F₃ and R → F₂₃₉; 04 checks it.
11. *Section 2.* The Taylor coefficients gᵢ, q_{ij} and h^m are printed in the text. The box centre is now ξ_c, since
    ξ₀ is the fixed point of (2.5). The cell count is the cover argument. (2.4) names the box of the sharp maximum and
    the five other boxes at that corner (05 now prints them), and counts 2998 divisions. (2.5) states the upper bound
    M ≤ 32733/32768 that the table gives.
12. *Wording.* "Unique factorization" is replaced by the divisibility theory of Sections 2–3. F′(θ) = 9π₂₃₉ with
    F(x) = x³ − 9x + 81, since f is the polynomial of (1.1); 04 printed the same slip. (1.1) uses (d, e, g). (4.2)
    encloses E directly. (4.6) says that ε is fundamental by (4.1). Theorem B's dependencies include proof.md §5.
13. *proof.md.* The header says that every check is exact and floating point only seeds searches. §1.2 cites (5.8) and
    (3.6) for the distinctness of the classes. §6.3 says "his class number (stated without a source) and his tabulated
    unit". §7 lists the new scripts.
14. *README.* Step 1 cites Lemma 1.3′; the elementary standard says "the divisibility theory".
15. *Provenance.* hypothesis_H.md §7 cites arXiv:1202.6019v1, §7 table, p. 16, read in this round; Taylor 1975/76 as
    the first proof that Z[λ] is norm-Euclidean, on Lemmermeyer's attribution; Lemmermeyer's 1989 table; LMFDB.
    literature_check.md has a new §6, and dated updates in §1.2 and §4.
16. *Shared files.* methods/norm.triples (what it does, the height in the complex case, the measured table length),
    methods/norm.triples.multiplier and methods/README (this card listed for the one-unit descent),
    targets/ghost_screen row 8 (Avanesov's dates, a link to this card), and triangular_tetrahedral/README (which still
    called this card PARTIAL).

**Process notes.**
- The workflow started on 2026-10-15 and was cut off when the author's session ended, after the two classes verifiers
  had finished. It was resumed on 2026-10-26 from the same run, but the cached results did not match the resubmitted
  prompts, so all nine agents ran again. The 15-October classes verifier of this card (93 checks, all passing) was
  overwritten by the 26-October one; its files are kept in
  [Archive/h_review_round3_2026-10-15/](../../../../Archive/h_review_round3_2026-10-15/).
- The units verifier ran one read-only git command, to read the line-ending setting, which its brief forbade; nothing
  changed.
- The referees downloaded arXiv:1202.6019 and the Cassels 1950 scan into their scratch folders to read them; none of
  it is stored in the repository.

After the round-3 fixes every (H) script of the card was re-run: 04–07, cover_search (the table regenerates byte for
byte) and the three verifiers. All pass.
