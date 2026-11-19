# ζ(7): consolidated holdings

Consolidated on 2026-10-26 from three branches, on Ishan's instruction of that day to bring the other branches into
`round/t009` and keep what is new, leaving out what is redundant or outclassed. The branches are not deleted. Anything
not imported here is cited as `<branch>:<path>` and can be read there.

## 1. Objective and status

- Target: prove that ζ(7) = Σ_{n≥1} n^−7 is irrational.
- Standing: Ishan has directed this work since 19 Oct 2026. The cloud Claude programme ran it as "lane Z" on 19 Oct;
  two Codex runs continued it from 20 to 25 Oct. It is NOT an objective of AGENTS.md. It has no ROUTE.md, no draft
  bundle and no ledger entry. Nothing here has passed `tools/accept.py` or been claimed to the referee.
- Status: OPEN. No branch proves ζ(7) ∉ Q, and none claims to.
- What is kept (RESULTS.md, seven entries):
  - R-Z01: a conditional improvement of the irrationality exponent of ζ(5);
  - R-Z02 and R-Z04: two new constructions of positive pure ζ(7) integrals on the seven-cube, with exact
    arithmetic;
  - R-Z05 and R-Z06: two proved sufficient criteria, one from a lattice and one modular. The first is shown to be
    equivalent to the target;
  - R-Z03 and R-Z07: reviewed exclusions and a counterexample.
  Every approach that was tested and stopped is in MECHANISMS.md, with its exact stopping point (21 entries).
- Review level: everything was reviewed at most by threads of the same run that produced it.
  - Branch B: two Codex reviews (individual and cohort) for most items.
  - Branch A: one internal review in the same cloud Claude session, for two items only.
  - Branch C: the captain's own scrutiny only.
  RESULTS.md records each item's actual reviews. It also records where our reviewer of 2026-10-26 ("rev4") or this
  consolidation recomputed something.

## 2. The benchmark preprint (citation and hash only)

- A. Fauzan, "ζ(5) is irrational", Zenodo, version 1, 17 Sep 2026, doi:10.5281/zenodo.22826419.
  - File: ZETA5_IS_IRRATIONAL.pdf, 476,337 bytes.
  - SHA-256 a42c05ccc46b9a7c36428913453c4617202c59109713ba20f464ab9ec640d0da. Recorded by branch C in
    `root/fauzan_functional/sources.json`. Branch A's `notes/05_fauzan_construction.md` records the same first and
    last eight hex digits.
  - MD5 30dafe4adbb13a79f202e94a046c8c8d, as listed on the Zenodo record (read 2026-10-26).
  - Licence: all rights reserved. The text is not stored here and must not be.
- Status: an unrefereed preprint that claims a proof.
  - Two third-party Lean 4 formalizations are reported, and neither was rebuilt here:
    - D. Romik, `github.com/danromik/zeta5-irrationality`, with the prime number theorem as an added axiom;
    - M. Firsching, `github.com/mo271/Zeta5`, with PNT imported from a dependency, standard axioms only, and a
      changed normalization.
  - Fauzan's own exponent bound μ(ζ(5)) ≤ 260 (his Corollary 1.2, Appendix C) is not formalized (branch A,
    `notes/05`).
  - Wikipedia and MathWorld describe the proof as claimed or announced.
- Correction applied throughout: branch A's track README says the preprint "proves" ζ(5) irrational, and branch C
  counts it among "valid proof controls". Read both as "claims to prove (unrefereed)".

## 3. Literature status (searched 2026-10-26; queries in the consolidation notes)

- Unconditional results:
  - Ball–Rivoal: infinitely many odd zeta values are irrational. T. Rivoal, C. R. Acad. Sci. 331 (2000).
  - Zudilin: one of ζ(5), ζ(7), ζ(9), ζ(11) is irrational. Russ. Math. Surveys 56 (2001).
  - Rivoal–Zudilin: two irrationals among ζ(5), …, ζ(69).
  - Fischler–Sprang–Zudilin, Compositio 2019; Lai–Yu, Compositio 156 (2020); L. Lai, arXiv:2501.05321: lower bounds
    of order √(s/log s) for the number of irrationals among ζ(3), …, ζ(s).
  - Lai–Zhou, arXiv:2103.00904: at least two of ζ(5), ζ(7), …, ζ(35) are irrational.
  - Fischler, "Linear independence of odd zeta values using Siegel's lemma" (version of March 2026): non-explicit
    small forms.
  - Brown–Zudilin, arXiv:2210.03391: cellular approximations 0 < |ζ(5) − p/q| < q^−0.86.
- New in 2026:
  - F. Brown, arXiv:2604.20741 (22 Apr 2026): a multi-parameter irrationality criterion for Mellin integrals.
  - Fauzan's preprint, Section 2 above.
  - L. Lai, arXiv:2610.06838 (5 Oct 2026; version 2 on 15 Oct): 1, ζ(3), ζ(4) are linearly independent, using
    Fauzan's determinant method and Dauguet–Zudilin's simultaneous approximations (arXiv:1401.5322).
  - F. Calegari's blog Persiflage, "zeta(5) is irrational" (galoisrepresentations.org, 24 Sep 2026). It reads
    Fauzan's construction through Prévost, Uvarov and Krattenthaler, and states that the method gives ζ(k) for
    k ≤ 5 "but not (at least directly) for larger k".
- Prior measurement of the Fauzan route at ζ(7): C. Toledo, "The Hankel-determinant method past ζ(5), measured"
  (carlostoledo.co, sources pinned 2026-10-05). The page describes itself as machine-derived and not peer-reviewed.
  It covers 323 families and finds the zero crossing near weight 5.3.
- Claimed proofs of ζ(7) ∉ Q, none accepted:
  - P. Anand, Zenodo 10.5281/zenodo.22920911, 23 Sep 2026. It follows Fauzan's route. Toledo's audit recomputes
    A₂₀₀ + U ≈ +4.86, where the paper prints −11.99, and finds that the proof fails.
  - A. Chattopadhyay, preprints.org, Jan 2026; and S. Suman, arXiv:2609.22316 (math.GM), "one of ζ(5), ζ(7), ζ(9)".
    Both are unrefereed and were not assessed here.
- No further claim about ζ(7) was found for September–October 2026.

## 4. Where everything lives

| What | Where |
|---|---|
| Kept results, with exact statements, true status, origin and checks | `RESULTS.md` (R-Z01 to R-Z07) |
| Every approach tested, where it stops, and its exclusions with their quantifiers | `MECHANISMS.md` (M-Z01 to M-Z21) |
| Cloud Claude lane Z, the curated mirror of `cloud_agent_claude/14_recurrences_wz_apery/` on `origin/claude/math-prime-cross-domain-g1vwsn` at 26c718ac | `experiments/A_claude_lane_z/` (index and corrections in its `CARD.md`) |
| Codex cloud run, Z side: curated mirror of `cloud_agent/claude_continuation/` on `cloud-agent/elementary-connections` at 48352281 | `experiments/B_codex_cloud_z/` (`CARD.md`) |
| Codex captain campaign, zeta side: curated mirror of `cloud_agent/goal_20261023/` on `codex/zeta7-pollock-research` at d263e40a | `experiments/C_codex_captain_z/` (`CARD.md`) |
| Consolidation checks (rev4 pointer and the new R-Z02 check) | `experiments/consolidation_2026-10-26/` |
| Everything else (Pollock material, briefs, logs, snapshots, third-party texts, owner originals, superseded or finite-only packages) | the three branches, which are kept |

Inside each experiment folder the source sub-tree is preserved, so the relative links between files of one branch
still resolve. Links from C to B (`../../claude_continuation/...`) and from B to its copy of A (`source_snapshot/`) do
not; each CARD.md lists where they now point.

## 5. The blocking quantitative obligations

Every route here produces nonzero real numbers Λ_n: a linear form A_n + B_nζ(7), or the value Q_n(ζ(7)) of an
integer polynomial of degree h_n. An irrationality proof needs, along an unbounded sequence:
(i) an explicit normalization N_n making N_nΛ_n an integer form;
(ii) Λ_n ≠ 0;
(iii) N_n·|Λ_n|·b^{h_n} → 0 for every fixed b.
Each family has a natural scale s_n. The table writes log|Λ_n| ≈ −r·s_n (real side) and log N_n ≈ a·s_n (arithmetic
side); success needs a < r together with nonvanishing.

| Id | Family (origin; entry) | Scale s_n | Real side r | Arithmetic side a | What is missing | Equivalence |
|---|---|---|---|---|---|---|
| O1 | Fauzan-type Hankel determinants at ζ(7) (A; M-Z01–M-Z03) | K² | ≈ 1.3855 for the equilibrium measure (numerical; independent of the weight s). Certified r ≥ 1.385008 at Fauzan's shape for ζ(5) only (R-Z01); for ζ(7) the same field applies, but the adapted §6.2 argument is not written out (A, `notes/09` §5) | true content exponent 1.931, 2.104, 2.149, 2.176 (K = 40…160, exact); extrapolated 2.23–2.26 | an extra e^{−0.87K²}. Measured totals a − r: +0.767, +0.846, +0.855, +0.862 (K = 40…160, Fauzan's shape); +0.555 is the minimum over 280 shapes; +0.744 for the best other design. Is there a family with a − r < 0? | none; "all families fail" is unproved, since it needs exact content upper bounds |
| O2 | Canonical seven-cube alternant Gram determinants (B; R-Z02, R-Z03) | n·h_n, h_n ≈ (841/129600)n⁶ | r ≤ 20.642: the real decay can be no faster (proved) | a ≤ 23.70 with PNT and a ≤ 26.20 elementarily (proved upper bounds) | a normalization of cost a < 20.642 for the full family, or a subfamily or lift that beats its own real bound. The thin corner δ = 1/1000 and fixed shifts τ ≥ 16 are proved to fail | none |
| O3 | Morris even family J_c = A_c/B_c + L_cζ(7) (B; decay from C; R-Z04) | n = 6c | r = α, 8.4080693218 ≤ α ≤ 8.4080693441 (proved) | a = limsup log S_c / n, with S_c = B_c/gcd(A_c, L_c); the generic bound is S_c ≤ lcm(1,…,13n)^7, i.e. a ≤ 91 | a < α, i.e. almost all of the generic large-prime denominators must cancel. Universal large-prime cancellation already fails at (c, p) = (1, 67), (1, 71), (2, 139) | none |
| O4 | Ball–Rivoal V(7,1) cubic-filter pure lattice (C; R-Z05) | n | r = σ = 9.1748678470 (proved) | log covolume/n → μ = 13.0634881489; limsup log D_n/n ≤ c∞ = 2.4451818849 | D_nλ_{2,n} ≤ e^{(σ−ε)n} along an unbounded sequence, i.e. log(λ2/λ1)/n < 2σ − 2c∞ − μ = 0.3958837752; the threshold is 0.8410656601 for filter degrees 2n ≤ d ≤ 3n+2 | EQUIVALENT to ζ(7) ∉ Q, given C's inputs (switch theorem). This is a reformulation, not a reduction |
| O5 | Modular prime-marked congruence lattices (C; R-Z06) | prime p | r = η* = log 2 − ¼ log(2 + 2√2) = 0.2995169887 (as a supremum over contours; proved), up to a factor ‖s_p‖∞ | a = limsup log‖s_p‖∞/p for the chosen marked vector s_p. The clearing factor 209·p⁷ is only polynomial, and the mark P ≡ 1728s₀ (mod p) gives nonvanishing under rationality | marked vectors s_p ∈ Λ_p∖Λ'_p with a < η* on an unbounded set of primes | sufficient only. Under rationality every marked vector would have ‖s_p‖ ≥ e^{(η*−o(1))p}/poly(p) |
| O6 | Brown's multi-parameter criterion (imported, arXiv:2604.20741, as read by A; M-Z06, M-Z07) | parameters n_i | Sup^rec of the image of the domain | d^7_{max n_i} | a positive Mellin family with Zariski-dense image, values in Q + Qζ(7) for all parameters and these denominators, with r ≥ 10 parameters, or r = 7 with Sup^rec < 0.153 | open: purity and max-type denominators have never been obtained together |

Notes on the table:
- O1: the K² rates are finite-K measurements with exact arithmetic. The asymptotic values are extrapolations, and
  "mapped exactly" in branch A's RESUME.md means exactly this.
- O4: the shape threshold comes from λ1λ2 ≍ covolume:
  - D_nλ2 ≈ exp(n(c∞ + μ/2 + shape/2)), which is below e^{σn} exactly when shape < 2σ − 2c∞ − μ;
  - c∞ is an unreviewed limsup bound, so the threshold is sufficient;
  - C's `root/fermat_descent` shows that rationality would force the shape rate to equal the threshold exactly.
- O2–O5 rest on author-level proofs: two Codex reviews for B's items, captain self-scrutiny for C's. rev4 recomputed
  the constants of O4 and O5 numerically.

## 6. Using this folder

- Read MECHANISMS.md before starting anything. An approach repeated unchanged will stop where its entry says.
- A route should first close completely on a known control: ζ(3) by Apéry, and Fauzan's claimed ζ(5) proof. C's
  matched controls (M-Z21) show that many ζ(7) "failures" are uninformative, because the same mechanism also fails at
  ζ(3) or ζ(5).
- Nothing here is a citable dependency yet. If Ishan makes ζ(7) a lane, a statement needed by a route must first be
  filed as a draft with two fresh reviews (AGENTS.md, rules 2 to 4).
