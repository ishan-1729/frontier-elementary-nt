# Octahedral numbers, bounded sums: every N ≥ 0 is a sum of eight octahedral numbers (given the three-square theorem)
Status: WRITTEN-AUDITED (two referees of the authoring cloud session, both PASS; not re-derived in round/t009)     Methods: not yet classified
Generalizes: none     Generalized by: none

Origin: branch `origin/claude/math-prime-cross-domain-g1vwsn`, tip 26c718ac, folders
`cloud_agent_claude/05_octahedral_bounded_sums/` (track 05: T1, T4) and `cloud_agent_claude/10_eight_octahedral/`
(track 10: T10, without its `seven/` subfolder), with their review files from `cloud_agent_claude/reviews/`. They are
copied unchanged into [source/](source/) and keep their own layout. Consolidated 2026-10-26.

This card records a result outside the Waring lane objective: the objective is seven octahedral numbers
(`research/waring/problem.md`), and nothing here changes a ROUTE.md row. O(n) = (2n³+n)/3.

## Statement
**Named hypothesis GL3(1 mod 4).** Every integer m ≥ 1 with m ≡ 1 (mod 4) is a sum of three integer squares.
This is the only case of the three-square theorem that the proofs use. It is the theorem proved on pp. 91–94 of
Dirichlet's 1850 memoir (Werke II), which is stored in the repository at
`research/waring/literature/sources/retrieved_by_local_web/dirichlet_1850_werke2_89_96.pdf`: every positive integer
that is neither 4k nor 8k+7 is a sum of three squares without common divisor. An m ≡ 1 (mod 4) is neither, so no
4-descent is needed. Dirichlet's upstream inputs (primes in progressions, Jacobi reciprocity, the reduction bound for
ternary forms) are the ones named in ledger W007 and in `research/waring/MECHANISMS.md`, M-W03.

**Theorem A (eight).** Assume GL3(1 mod 4). Every integer N ≥ 0 is O(n₁)+…+O(n₈) with integers nᵢ ≥ 0.
- (a) *Large range* (T10, Theorem 1). For N ≥ X0 = ⌈41·8³²/10⌉ = 324,835,466,308,483,784,133,530,196,378
  (about 3.25·10²⁹), N = Σᵢ₌₁³ [O(c+dᵢ)+O(c−dᵢ)] + O(u) + O(v) with c a power of 2, |dᵢ| ≤ c and u, v ≥ 0.
- (b) *Small range* (no hypothesis). Every N with 0 ≤ N ≤ U8 = 20,626,541,664,582,647,929,372,758,803,668,289,246,952
  (about 2.06·10⁴⁰) is a sum of at most eight octahedral numbers.
- (b) needs no new computation. Ledger W011 proves seven terms for every N ≤ 1,971,099,043,927,191,059,080,960,614.
  Its accepted general implication (paragraph 2 of its claim), applied once more with s = 7, L = 65,285,684, gives
  J = 31,393,463,045,092 and U8 = L + O(J+1) − 1. See [finite_part_check.py](finite_part_check.py) and its output
  [finite_part_check.out.txt](finite_part_check.out.txt).
- Since U8 > X0, (a) and (b) together give Theorem A.

**Theorem B (nine; T1).** Assume GL3(1 mod 4). Every N ≥ N1 = 4,395,799,615 is
Σᵢ₌₁³ [O(c+dᵢ)+O(c−dᵢ)] + O(A+τ) + O(A−τ) + O(x). Here c = 2ᵃ, A = 16α with α odd, 0 ≤ τ ≤ A, 0 ≤ x ≤ 63 and
|dᵢ| < c. With W011 (seven terms for every N ≤ 1.97·10²⁷ > N1), every N ≥ 0 is a sum of nine octahedral numbers.
Theorem B is weaker than Theorem A. It is kept because its carrier lemma (T1, Lemma 2 and Corollary 2.1) is used by
case C1 of Theorem A, and because its large range is checked by the Lean kernel.

**Lemma C (T4; limits of the method).** 6(O(u)+O(v)) = 2(u+v)(2(u²−uv+v²)+1), so v₂(O(u)+O(v)) = v₂(u+v). As a
consequence, two free values cannot compress 2-adic residues: covering every class mod 2ᵏ needs indices up to
2^(k−2). Three free values compress (Proposition 4.3).

## Proof 1 (Theorem A)
- Files:
  - [source/10_eight_octahedral/theorems/T10_eight_octahedral.md](source/10_eight_octahedral/theorems/T10_eight_octahedral.md),
    §§1–7: 2-adic exponential sums along the curves O(u)+O(v) ≡ R (mod 2ᵐ) and F_r(σ,τ) ≡ R″ (mod 2ᴹ), by
    elementary 2-adic stationary phase (Lemmas 2.2, 2.5); counting in boxes with Fejér weights (§§3–4); the carrier
    without small term when v₂(N) ≥ 5 (§5, from T1 Lemma 2); scales, boxes and windows (§6); four cases on v₂(N)
    (§7).
  - [source/05_octahedral_bounded_sums/theorems/T1_nine_octahedral.md](source/05_octahedral_bounded_sums/theorems/T1_nine_octahedral.md):
    Lemmas 1.1, 1.2, 2 and Corollary 2.1, which T10 imports.
  - The scripts of `source/10_eight_octahedral/scripts/`. `constants_check.py` and `identities_check.py` are part of
    the proof: every numerical inequality is checked in exact rationals at c = 2³², and every polynomial identity in
    sympy. The others are consistency checks and the executable construction. Their outputs are in
    `source/10_eight_octahedral/evidence/`.
  - The finite part is (b) above. It replaces T10 Corollary 2's use of track 05's T2(e), whose ranges W011 outclasses.
    Track 05's T2 is not imported.
- Lean (partial). [source/10_eight_octahedral/lean/EightCore.lean](source/10_eight_octahedral/lean/EightCore.lean)
  - Declarations: `shift_half`, `shift_lin`, `pair_param`, `sep_A`, `I_le_8`, `sep_B`, `J_le_4`, `carrier_C2`.
  - These are the arithmetic cores: Lemma 1.3, the identity of Lemma 1.4, Step 1 of Lemmas 2.2 and 2.5 with the
    bounds I_n ≤ 8 and J_n ≤ 4, and Lemma 5.2.
  - Not formalised: the exponential sums, the Fejér counting, the sizes and the case assembly.
  - Built in the branch's environment: Lean v4.19.0 with that branch's pinned Mathlib; see
    [source/tools/LEAN.md](source/tools/LEAN.md) and evidence E9. Referee B recompiled it. It was not rebuilt in
    round/t009, and it is not part of the verified project in `lean/`.
- Audit:
  - [source/reviews/10__A_exponential_sums.md](source/reviews/10__A_exponential_sums.md) (scope §2) and
    [source/reviews/10__B_counting_cases.md](source/reviews/10__B_counting_cases.md) (everything else) both PASS.
  - The referees are separate agents of the same cloud Claude session as the author: not the repository's two-review
    standard, and not a referee report of `docs/referee/`.
  - Points raised by referee B:
    - the below-X0 remark is wrong for case C2 (it works from a = 25);
    - the novelty note must credit Basak–Dong–Saettone–Zaharescu §9.
  - The round/t009 branch review of 2026-10-26 (rev1) rated Theorem A CONSISTENT, not re-derived
    ([branch_review_2026-10-26.md](../../../shared/consolidation_2026-10-26/branch_review_2026-10-26.md), §3.1, item 3).
- Provenance and novelty:
  - Searches: the branch's (`source/10_eight_octahedral/literature/SEARCH_LOG.md`), its referees', and the
    consolidation's of 2026-10-26, listed below.
  - No earlier statement that every integer is a sum of at most K octahedral numbers for an explicit K was found:
    - Brady (J. London Math. Soc. 93 (2016) 244–272) proves seven for N > e^(10⁷);
    - Dong–Nguyen–Zaharescu (arXiv 2412.13974) state only Brady's bound for {3,4};
    - Basak–Dong–Saettone–Zaharescu (arXiv 2608.08019; IMRN 2025, rnaf180) treat icosahedral and dodecahedral
      numbers only.
  - The architecture (three symmetric pairs, two free values, a curve counted in boxes) closely adapts
    Basak–Dong–Saettone–Zaharescu §9, which uses a prime modulus with Weil–Bombieri. Here an elementary 2-adic
    stationary-phase bound replaces Weil. The pair backbone (three symmetric pairs reduced to three squares) is
    Maillet 1895/96.
  - By the novelty criterion of 2026-10-15 (Ishan): no earlier proof, or rewording of one, of this theorem, an
    equivalent or a generalization was found. Reuse of parts proved elsewhere does not count against it.

## Proof 2 (Theorem B)
- Files: [source/05_octahedral_bounded_sums/theorems/T1_nine_octahedral.md](source/05_octahedral_bounded_sums/theorems/T1_nine_octahedral.md).
  - The mechanism is a 2-adic carrier pair at centre 16α plus one small term x ≤ 63, with a monotone-gap covering
    across dyadic scales (Lemma 3).
  - The scripts used by T1 (`interval_chain.py`, `construct_nine.py`, `check_gl3_small.py`) are in
    `source/05_octahedral_bounded_sums/scripts/`, with outputs in `.../evidence/`.
  - The small range comes from W011, not from track 05's T2. T1's own text cites T2, which is not imported.
- Lean: `nine_octahedral_large (hGL3 : GL3) (N : ℤ) (hN : 4395799615 ≤ N)`, in
  [source/05_octahedral_bounded_sums/lean/Octahedral.lean](source/05_octahedral_bounded_sums/lean/Octahedral.lean).
  - It states the nine-term conclusion, with only |dᵢ| ≤ c exported. GL3 is the explicit hypothesis
    ∀ m ≥ 0, m % 4 = 1 → sum of three squares.
  - Evidence E5. The reviewer compiled it and added an ℕ restatement (`source/reviews/evidence/05_T1T4/`).
  - Not rebuilt in round/t009.
- Audit: [source/reviews/05__T1_nine_octahedral.md](source/reviews/05__T1_nine_octahedral.md): PASS, MATCH (one
  in-session referee). rev1, 2026-10-26: the ranges were CONFIRMED and the theorem is CONSISTENT.

## Limits (Lemma C)
[source/05_octahedral_bounded_sums/theorems/T4_two_term_valuation.md](source/05_octahedral_bounded_sums/theorems/T4_two_term_valuation.md),
reviewed PASS ([source/reviews/05__T4_valuation.md](source/reviews/05__T4_valuation.md)).
- Lemma C obstructs only the deterministic monotone-gap method. Eight-term representations of the same shape exist
  for every odd N tested, and T10 handles them by averaging over the curve.
- Seven terms (the lane objective) would need a single free value. No averaging is available then, and this method
  does not reach it.
- The seven-term work is in `research/waring/experiments/CONS26_single_form_17_18_19/`.

## What this does not prove
- Pollock's conjecture (seven terms), or anything for the Waring lane's obligations.
- Anything in the large range without GL3(1 mod 4).
- A Lean proof of Theorem A: only its arithmetic cores are formalised.

## Stale lines in the frozen source
- `source/10_eight_octahedral/README.md` says "Independent review: NOT YET DONE". The two referee reports above were
  filed later (rows "10" of the branch's `reviews/REVIEW_LEDGER.md`, which is not imported).
- T1 and T10 cite track 05's T2 for their finite ranges. This card uses W011 instead.

## Searches of 2026-10-26 (novelty)
WebSearch:
- "every positive integer is a sum of at most octahedral numbers explicit bound";
- "\"octahedral numbers\" sum of eight OR nine every integer theorem";
- Kim, "On regular polytope numbers" (Proc. AMS 131 (2003) 65–75): the abstract only; the full text was refused
  (HTTP 403);
- the OEIS for the twelve numbers needing seven terms; A282349 gives no verification range.

Pages fetched:
- arXiv abstracts 2608.08019, 2608.03309 and 1904.02787;
- the HTML introduction of 2412.13974.

No match.

## Open questions
- Lower X0. The construction succeeds numerically from about 10^11.5 (T10, notes/N1).
- Formalise the exponential sums, the Fejér counting and the case assembly.
- Have the eight-term proof reviewed by the repository's two-review process if it is ever needed by a lane.
- Does the same 2-adic method give eight values of every cubic that permutes Z₂ and has a three-square packet,
  as a Watson-type theorem?
