# Review record: triangular_tetrahedral

All reviewers below are Claude subagents of the same model family as the author. So this is in-session review, not
an independent audit. After rounds 1–2 the card was WRITTEN (Theorem A) and PARTIAL (Theorem B). After rounds 3–4,
which proved (H), it is WRITTEN (in-session review only).

## Round 1: verification of the computations (2026-10-15)
Two workflows ran, with six agents. Each worked alone, wrote its own script, and reported in a structured form.

| Agent | Task | Verdict | What it found |
|---|---|---|---|
| fields | wrote [00_fields.py](scripts/00_fields.py) | 186 checks pass | Fields, maximal orders, markers and units confirmed. Avanesov's fields are ours, with explicit isomorphisms. His six equations are six of ours. His units generate our group modulo ±1. His p. 413 remark on discriminants is false for one field. |
| coverage audit | general Lemma 3.2 with proof; audit of 01; wrote [verify_coverage.py](scripts/verify_coverage.py) | no bug affecting coverage | Two rigour gaps and one cosmetic issue, all fixed (below). The independent enumeration with the weaker linear bound finds the same 16 orbits. By-product: O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩. Provenance: Koecher's reduction theory. |
| reconstruction | wrote [03_reconstruction.py](scripts/03_reconstruction.py) | 880 checks pass | 28 points on the monic model and 22 on 1611.a1; C(n,2) = C(m,3) values {0, 1, 10, 120, 1540, 7140}. New fact: u² ≡ F (mod 6). It names N1 and N2 as further inputs of the reduction. |
| literature | sweep for other proofs | — | Proofs by Stroeker–de Weger (1999) and Gallegos-Ruiz et al. (2020). No elementary proof found. Uchiyama (1984) gives a congruence-only two-index argument for 8604.a1. |
| exclusion re-implementation | wrote [verify_exclusion.py](scripts/verify_exclusion.py), with no imports from core, lifting2 or 02 | CONFIRMED | All 16 orbits match. Negative controls bite. λ and ν(J) agree in three bases at every zero used. |
| exclusion audit | line-by-line audit; wrote [verify_bruteforce.py](scripts/verify_bruteforce.py) | CONFIRMED | Exact box \|a\|,\|b\| ≤ 40 and modular box ≤ 1000: only the certified zeros. Direct search \|m\| ≤ 10⁵: only the known solutions. 810 planted two-zero systems: Proposition 17 fails on all, as it must. |

**Fixes made after round 1.**
1. *Signature check.* [core.py](scripts/core.py) `Embeddings` now asserts by Sturm's theorem that f has exactly two
   real roots. The coverage audit noted that the S, P formula for |σ₃|² silently needs this.
2. *Enumeration ranges.* [01_coverage.py](scripts/01_coverage.py) now certifies exactly that each float-guided
   Fincke–Pohst range contains every solution: both endpoints fail the exact test, and they lie on either side of the
   centre.
3. *Dead code.* Removed from 01, together with a comment that did not match the code.
4. *Unit-group by-product.* 01 now prints it.
5. *Proposition 17 instances.* [02_exclusion.py](scripts/02_exclusion.py) now prints every instance used, including
   those behind sibling classes: (13, (0,1)) in orbit (1,1) c = 1, and (25, (−1,1)), (25, (2,1)) in orbit (1,1)
   c = 9. Both verifiers had flagged that the printed certificate was incomplete; the step itself was right.
6. *Assertions.* 02 now asserts the facts that the correspondence "zeros ↔ solution pairs" needs: unit norms 1, the
   shape of the bases (marker criterion), and −1 ∉ ⟨ε₁, ε₂⟩ (−1 is not a residue modulo 4 or 8). Its docstring now
   points to proof.md.
7. *Basis of Λ_q.* proof.md §4.3 said that λ and ν(J) do not depend on the basis of Λ_q. That is true when ν(J) < e,
   which covers every instance used. In general it is false: the re-implementation found different ν(J) at failing
   instances. The text is corrected.
8. *Inputs of Theorem B.* proof.md §0.3, §1.2 and README now name N1 and N2 beside (H).
9. *Threshold instances.* Two instances meet Proposition 17's bound with equality: (16, (−1,1)) and (25, (2,1)). The
   strict inequalities inside its proof still hold. This is recorded in proof.md §4.5.

**Process notes.**
- The fields agent ran one read-only `git diff`, which its brief forbade. Nothing changed.
- The exclusion re-implementer's `py_compile` left `scripts/__pycache__/verify_exclusion.cpython-313.pyc`. It is
  gitignored.

All scripts were re-run after the fixes, and all pass (outputs in [scripts/](scripts/)).

## Round 2: referees on proof.md (2026-10-15)
Three referees worked alone, each with a distinct lens: coverage (§§2–3 with core.py and 01), exclusion (§4 with
lifting2.py and 02, re-run), and global (§§0–1, 5–7, README, the literature check and the method card, with
Avanesov's page images). They ran their own scratch checks. These included:
- 45,000 random tests of the interval arithmetic;
- an mpmath recomputation of every vertex and of b_min and T;
- a brute-force enumeration of the boxes;
- greedy descents from random orbit elements;
- dependent unit pairs, as a negative control.

**Result.** 28 findings: 0 critical, 2 major, 26 minor. No finding changes Theorem A, Theorem B or any certificate.
The referees confirmed that every number in tables 3.3 and 4.5 matches the outputs, apart from one rounding slip.

**Major findings, both fixed.**
1. *Remark 3.2(i) (independence of ε₁, ε₂ from the cone certificate) was asserted without proof.* Several statements
   depended on it: injectivity, −1 ∉ ⟨ε₁,ε₂⟩, and the fundamental-unit by-product. A proof is now written in
   Remark 3.2(i), as the referee suggested.
   - If the log size vectors lay on a line, a point of P⁺ with a zero coordinate would exist, by e^y − 1 ≥ y.
   - The referee checked numerically that dependent pairs such as (ε₁, ε₁²) do fail the certificate.
   - Theorem A never needed independence: Z_w = Z₀ accounts for every solution, and −1 ∉ ⟨ε₁,ε₂⟩ is also certified
     modulo 4 or 8.
2. *§6.3 and the literature check said the units were "not proved fundamental".* This contradicted Remark 3.2(iii).
   Both now say that the units are not assumed fundamental, and that fundamentality is a by-product. That by-product
   confirms Avanesov's p. 411 claim that his units are fundamental.

**Minor findings, by group (all fixed unless stated).**
- *Notation in §§3–4.* There were collisions of c₀, c₃ with c; S, P with the set and the polygon; b with the
  exponent; b₀ with the class representative; and "T₂-minimal" with E. Fixes:
  - the coefficients are now a₀..a₃;
  - the pair is t₃, n₃ and the set Σ_z;
  - the vertex minimum is b_min;
  - Lemma 4.1 uses (h, k, j);
  - ε^{(a,b)} is defined.
- *Unstated standard facts.* These are now stated:
  - Minkowski's convex-hull theorem, with vertices as feasible basic solutions;
  - the Vandermonde reason that the embedding is an isomorphism;
  - the irreducibility behind the marker lemma.

  Step 4 now uses concavity of log π, so the chord argument and its tacit 2-dimensionality are gone.
- *Enumeration bullet.* The exact Gram matrix G*, the symmetrized Ĝ and the positive pivots are now named.
- *Rounding.* 12.8·√c is corrected to 12.7·√c.
- *Koecher 1960* is now marked as not accessed, with the attribution following Porter–Bali–Leibak.
- *"Two instances sit exactly at the bound" was false.* 12 of the 14 do. The margins are now stated exactly. The
  referee checked that (16, (−1,1)) can be dropped, while (25, (2,1)) is needed.
- *"On the basis" paragraph.* It is restated: λ is invariant when λ < e, ν(J) is invariant when ν(J) < e, and
  λ = 0 at every zero used.
- *Phantoms (§6.2).* Orbit (1,1) c = 1 never uses q = 4 or 8. Its phantom died as a sibling at q = 13, because −1 is
  not a residue modulo 27·13. Now stated.
- *Citation for −1 ∉ ⟨ε₁,ε₂⟩.* §4 now cites 02's own check (q = 4, 4, 8), with 19/11/47 as the cross-check.
- *Counting.* "Solutions ↔ Z_w" is now "solution pairs ↔ Z_w"; "bijection onto the class" is now "onto the exponent
  class".
- *Lemma 4.4* now excludes q = 2.
- *Vacuous modulus.* q = 3 after q = 9 in orbit (1,1) c = 1 adds nothing. It is footnoted in table 4.5, and the
  certificate is left as it is.
- *Overstatements in README and the method card.*
  - Sibling removals are named as such.
  - The sibling rule requires Proposition 17 at the same q.
  - Each zero may use its own modulus.
  - The card now restricts G to norm +1, or enumerates norm ±c.
  - The card says §3 proves the weights-1, signature-(2,1) case.
- *(H).* It was stated two ways. README now states proof.md's (H) verbatim, with the explicit classes, and says
  "(H), implied by rank ≤ 2" throughout. N5 (Mordell–Weil) is named. Avanesov's reduction is said to cover only the
  points of his model.
- *§1.* The homomorphism is proved by the one-line identity; the 109-pair check is extra. "Re-run" is corrected to
  "independently re-derived". The rank claim in §0.1 is hedged as an mwrank value.
- *§0.3.* N1 is a determinant identity and N2 the standard integrality lemmas. Theorem A's list now includes Sturm's
  theorem and the real-variable facts of Lemma 3.2.
- *Avanesov's coordinates.* Now X_E = 6X, Y_E = 36Y, consistent with 00_fields.
- *Avanesov's case analysis* (checked again against pp. 414–419). Corrected in literature_check.md §1:
  - cases 3 and 6 fail mod 2 and case 5 fails mod 4;
  - cases 4 and 8 fail by 3f₂ = ±1, and case 8 is only asserted;
  - case 1 uses coprime factorization, not Lemma 1;
  - Mordell 1923 is cited only for finiteness;
  - his Lemma 1 as printed misses the primitive solutions with v odd, such as (1,1,2). This is harmless, because he
    uses it with v even.
- *"Every lemma is proved"* overstated the case. It now reads: all twelve Thue equations are proved, and every
  computation is re-runnable.
- *lifting2.py docstring.* It now describes the precision cap correctly.
- *Not done.* Avanesov's claim (5), that ε̄₁, ε̄₂ generate Z[η]^× modulo ±1, was confirmed by a referee's scratch
  check. No script in this folder records it, so the literature check calls it unchecked.

All scripts were re-run after round 2, and all pass.

## Round 3: hypothesis (H), exploration and build (2026-10-15)
A side session proved (H) for integral points: [hypothesis_H.md](hypothesis_H.md). The same caveat applies: all agents
were Claude subagents of the author's model family.

**Exploration by the lead** (scratch scripts, not kept; their results are re-proved by the scripts below).
- *The suggested route, a Euclidean algorithm, fails.*
  - A float covering of the torus with small translates left residual points near two places, ±ξ*.
  - ξ* is fixed by ε modulo Z[λ], and a complete search shows min |N(ξ* − κ)| = 121/109.
  - This is now Remark 2.5, proved in 05.
- *What works.*
  - A Dedekind–Hasse division with multipliers {1, 2}. Float prototypes needed about 830 boxes in Minkowski
    coordinates and about 1900 in (a,b,c) coordinates.
  - The norm condition leaves 8 classes, and a congruence modulo 8 removes the coset εT, both found by quick checks.
  - The box principle needs a table only up to a constant; the first estimate was 171.

**Build workflow** (five agents in parallel, each writing its own files):

| Agent | Wrote | Result | Notes for the text |
|---|---|---|---|
| ring and units | 04_ring.py, 06_units.py | 87 + 33 checks pass | Q = \|σ₂\|² is only positive *semi*definite on real triples, so the plan was corrected; T₂ is definite (determinant 716). Q = σ₁(α*) exactly. Unit descent with threshold 114; the units with T₂ ≤ 114 are ±1. |
| division | 05_dh_cover.py, dh_boxes.txt, cover_search_K.py | all pass | 398 boxes (an optimal split tree among five rules), largest Taylor bound 8191/8192, sharp maximum 4093/4096. Remark 2.5 proved rigorously; of the 654 ε-fixed classes only ±ξ* have minimum ≥ 1. Planted errors (box removed, κ or s changed, corner moved) all make 05 fail. |
| classes | 07_classes.py | all pass | The coset εT already dies modulo 4, with a hand proof. Coprimality cannot be read off norms: at 7 of the 14 points the S-free parts have norms with a common factor, so the script uses explicit Bézout identities. (3.2′) is needed for several primes at once. |
| box principle | 08_ufd_box.py | all pass | Another basis and a Bernstein bound give C₀ = 535/24 < 23 in place of the plan's 171, so the table needs only p ≤ 22 and 179. The descent was run for all primes 23 ≤ p < 3000 (781 lattices). |
| literature | (report only) | — | Cavallar–Lemmermeyer (1998) and Lemmermeyer's survey list −716 as "N 121/109", by computation. No multiplier certificate for −716 and no earlier elementary (H) was found. See literature_check.md §5. |

**Lead's fixes before verification.**
- §1.6 now says Q is semidefinite, and states the T₂ forms.
- Remark 2.5 now uses the sharp constant C₀ ≤ 113.6238: the rounded 114 of (4.3′) would not fit inside T_e.
- (3.2′) was added for products, and §5.4 now cites it.
- The step in §5.7 is now modulo 4, with the hand proof; the counts modulo 8 are recorded as well.
- Line endings: a Python edit had turned six LF files into CRLF. They were restored to LF, and the diffs are minimal
  again.

## Round 4: verification and referees of hypothesis_H.md (2026-10-15)
Seven fresh agents worked in parallel. Each worked alone and was told nothing of the others. The four verifiers
imported nothing from the author's scripts. The session was cut off while four of the seven were still running. The
lead then:
- re-ran verify_H_box.py to completion (its script was complete, its run had been killed);
- re-ran the claims referee as a fresh agent;
- used the stored outputs of verify_H_units and verify_H_dh, which had both finished.

| Agent | Wrote | Verdict | Method and numbers |
|---|---|---|---|
| units | verify_H_units.py | CONFIRMED, 101 checks | Norm form three ways (resultant, Sylvester, Newton). Unit theorem by a different route: a fundamental domain for ε. All 291 lattice points with \|σ₁\| ≤ 606 and Q ≤ 1 were enumerated, 290 of them with \|σ₁\| ≥ 1; the units there are ±1 (and ±ε⁻¹ on the excluded end). Brute force on \|a\|,\|b\|,\|c\| ≤ 60 finds ±1 and ±ε; ε⁻¹ lies outside that box. Control: ε² in place of ε exposes ε. |
| division | verify_H_dh.py | CONFIRMED, 61 checks | Every box re-bounded by the factored form max\|σ₁\|·max Q, which is not the Taylor bound of 05. 189 boxes needed splitting (2722 sub-boxes). The largest value is 4093/4096 up to 2⁻²⁰. Cover checked by volume, disjointness and a sweep. 10⁵ random points, 2000 divisions. Remark 2.5 by a different finiteness argument (fundamental domain for ε on the coset): the minimum is 121/109, and of the 654 fixed classes only ±ξ* fail. Planted errors are caught. |
| classes | verify_H_classes.py | CONFIRMED, 99 checks | Own arithmetic: polynomial division, Leibniz determinant, Cramer division, q-adic square roots. The 8 admissible μ₀ found by enumerating all 128 vectors. Counts mod 4: 0,0,0,0 and 8,16,16,8; mod 8: 0,0,0,0 and 32,48,48,32. The identities of (5.8) hold. Distinctness by characters at other primes. A search over \|X\| ≤ 2·10⁶ gives exactly the 14 X, each with X − θ = μ₀S² explicitly in the predicted class. Planted errors are caught. |
| box principle | verify_H_box.py | CONFIRMED, 103 checks | The constant by exact critical-point analysis (vertices, edges, faces, interior): max \|F\| on the box is 49531/4374 + (439/2187)√1756 ≈ 19.74 < 535/24. Lemma 6.3 for all 9592 primes below 10⁵. Lemma 6.2 for p ≤ 60 and p = 179 by all residues. The table and its HNFs, the descent for 23 ≤ p < 3000, and 179 = unit·π₁₇₉²·Q₁₇₉. |
| logic referee | (scratch only) | 11 findings, all minor | Theorem H is proved. Theorem B no longer needs the rank, N5, N1 or N2. Its dependencies outside this file are Theorem A, §1 (including the conic parametrization) and §5 of proof.md. |
| rigor referee | (scratch only) | 8 findings, all minor | Re-ran 04–08: identical outputs. Re-derived the key certificates. No numerical or logical gap. |
| claims referee | (scratch only) | 1 major, 10 minor | Re-ran 04–08: identical outputs. Spot checks of the numbers. Section 2 is free of the §6 kernels. The major finding is a provenance error outside the proof (fix 11). |

**Fixes made after round 4.**
1. *(4.3).* "Both ρ's increase with E" was wrong: ρ_lo decreases. The bound used, ρ_lo > 5453·10⁻⁹, was right and is
   checked by 06 with the correct identity. Corrected. (logic, rigor)
2. *Remark 2.5.*
   - The E of K₃ is now defined: E = Q(ε), and K₃ is symmetric under E ↦ 1/E, as 05 uses.
   - K₃ is bounded above by its value at the end of the enclosure.
   - The proof that T_low ≤ T₂ is written out. (rigor, logic)
3. *Remark 2.5, last paragraph.* "The multiplier is needed" now reads: a plain Euclidean division is impossible;
   Section 2 uses a multiplier, and Section 6 avoids division. (logic)
4. *§3 preamble.* The proofs use (1.2), (1.3a), (1.4) and (3.1). The text now gives the numbering map ((3.6) here is
   mordell_1963 (3.7)) and the changes. (3.5) now has its base case, and (3.6) writes out the argument and the
   characteristic polynomial. (logic)
5. *Citation.* [04: 1.3] is corrected to [04: 1.2]. (logic, rigor)
6. *§5.4.* The reason why no prime of 𝒮 divides the final quotients is corrected. The exponents are renamed e_a, e_b,
   e_c, e₁₇₉ and their parities ν, so that they no longer clash with the coordinates (a, b, c) or the gcd g. In (5.8),
   r is renamed τ. The text now says why S ∈ R: ε⁻¹ ∈ R. (logic)
7. *§7.* R ⊆ Ξ_μ = −Ξ_μ, the line (μξη)² = (μξ²)(μη²), and "Ξ_μ is a full lattice" are added. The descent script's
   N1/N2/N5 labels are now said to be superseded for 1611.a1. (logic)
8. *proof.md §1.3: the conic parametrization is now proved (Lemma 1.3′).* Before, it was only checked by script. The
   lemma shows that every primitive point of a nondegenerate integral ternary conic with a primitive point P₀ is
   ±X(m,n)/g with (m,n) coprime. The tangent case is excluded by the nonzero determinant. The gcd bound g | 36 is
   referred to the identities the descent script checks. proof.md §1.1 now cites (3.6) for O_K = Z[λ]. (logic,
   finding 9: the one gap that concerned Theorem B rather than Theorem H)
9. *Wording.*
   - "All arithmetic is exact" now reads: every check is exact, and floats only seed searches.
   - "Ten coefficients" now reads Σ|n_m| = 113.
   - The sharp-bound remark now says "each box on its own closed box".
   - §6 now says "no real root of f".
   - §8 gives 239 with unit weights and 171 with optimized weights for the power basis.
   - proof.md §6.3 and GENERALIZATIONS now write M ≥ 121/109.
   (rigor, logic)
10. *README.* Step 1 cites Lemma 1.3′. (logic)

11. *Major (claims referee): Hurwitz and integer multipliers.*
    - The method card and §8 said that Hurwitz (1919) guarantees a finite set of integer multipliers for every field
      of class number 1. That is false for Dedekind–Hasse division. In Z[(1+√−19)/2], no integer s works at ω/2: the
      minimum of |N| on its coset is 5/4, and an even s lands in R.
    - Hurwitz's lemma allows jx = y. Both texts now say so, and say that a certificate must decide each set of
      multipliers, including the points with sξ ∈ R.
    - norm.triples.md no longer generalizes {1, 2} to every ring of class number 1.
    - §8 now says what Cámpoli and Wilson use at the 2-torsion points (θ or θ̄).
    - Theorem H is not affected: (2.1)(b) treats all seven half-lattice classes of Z[λ].
12. *Stale statements.* The following still said that Theorem B was conditional or that (H) came only from the rank.
    They now point to hypothesis_H.md:
    - this file's header;
    - targets/literature_screen_2026-10-15.md and targets/rank2_signatures_2026-10.md;
    - methods/ghost.unit.multi_unit.md;
    - GENERALIZATIONS.md §4 item 7.
13. *Provenance.*
    - Stroeker–de Weger's reading marks are made consistent across literature_check.md §2, the not-accessed list and
      §5.4, and the mwrank detail is recorded in §5.4.
    - "Billevich" now reads "computed by Billevich's method".
    - Dedekind's criterion is cited as Werke II 373–375 [not accessed].
    - The citation is now Satz I (box lemma) and Satz II (finiteness).
    - −335 is added to the two-stage list.
    - "Unique factorization by Sections 2–3" now reads "the divisibility theory that class number 1 provided".
    - The contrast with Cavallar–Lemmermeyer is fairer: they made exact checks at the exceptional points, and the part
      of their value that is not a proof is the upper bound.
14. *Elementarity disclosure.* The header and the README now say that Section 6 is not part of the proof of Theorem H,
    and that its P(p, g) are prime ideals in modern language. Only Sections 1–5 meet the no-ideal-theory standard.
15. *Real numbers.*
    - (1.6) now proves the uniqueness of the real root by 04's factorization, without critical values.
    - (4.3) cites 06's endpoint identities, with convexity as an aside.
    - Remark 2.5 uses finiteness instead of a limit, and the header says that Remark 2.5 is not used by Theorem H.
16. *Numbers.*
    - "Largest |κᵢ| is 22."
    - The five other boxes at the sharp point carry the referee's values: 2829/4096 three times, −3267/4096 and
      −3603/4096.
    - "3302 random and planted pairs, 3301 with β ∤ α."
    - "171 with optimized weights."
    - "781 prime lattices, 428 of them built by the descent."
17. *Notation.* T is now defined in (5.6), and εT is described as the products up to the square ε².
18. *proof.md header.* It now reads: every check is exact, and floats only seed searches.

All scripts 04–08 and the four verifiers were re-run after the session cut, with output to a scratch folder. All exit
0, and the outputs agree with the stored ones except for timing lines. Status stays WRITTEN
(in-session review only): every reviewer was a Claude subagent of the author's model family.
