# ζ(7) lane: kept results

Seven entries, consolidated 2026-10-26 from three branches:
- A = `origin/claude/math-prime-cross-domain-g1vwsn:cloud_agent_claude/14_recurrences_wz_apery/`;
- B = `cloud-agent/elementary-connections:cloud_agent/claude_continuation/`;
- C = `codex/zeta7-pollock-research:cloud_agent/goal_20261023/`.

Imported copies are under `experiments/A_claude_lane_z/`, `experiments/B_codex_cloud_z/` and
`experiments/C_codex_captain_z/`, with the same relative paths. The imported files keep their original status words
("PROVED" and the like). The statuses below replace them.

Status vocabulary:
- **author proof**: a written proof by the producing thread, not checked by anyone else unless stated.
- **reviews**: independent reviews that exist, with who wrote them. "Codex review" means a fresh thread of the same
  Codex run. None of these meets AGENTS.md rule 3 (two fresh local threads, the individual reviewer of another model
  family).
- **certificate**: exact or interval computation (finite).
- **measurement**: a finite exact or numerical computation, with no asymptotic claim.
- **conditional on X**: X is an unproved, unrefereed or imported hypothesis.
- **reformulation**: the missing hypothesis is as strong as ζ(7) ∉ Q (AGENTS.md rule 6).

Nothing below is referee-ACCEPTED or a lane draft. "rev4" is our branch reviewer of 2026-10-26
(`research/shared/consolidation_2026-10-26/review_checks/rev4/`).

---

## R-Z01. Irrationality-exponent bounds for ζ(5), conditional on Fauzan's Theorem 2.1

**Statement.**
Notation is Fauzan's: K = 40n, N = 3n, h = 37n, α = 3/40, λ = 37/40, field V, F_K = S_KΔ_K, Q_{K,M} = m_{K,M}F_K.
Let (F) denote Fauzan's Theorem 2.1 for M = 100000: Q_{40n,M} ∈ Z[X], deg Q = 37n, Q(ζ(5)) > 0.

- (a) *Certified comparison measure* (A, `theorems/T6_zeta5_real_exponent.md`, Thm 6.1). The 160-layer nested arcsine
  measure ρ' in `data/rho_zeta5_qp160.txt` satisfies
  U' := λM₀' − I(ρ') + C* ≤ −1.385008.
  For comparison, Fauzan's measure gives −1.366996 and the numerical equilibrium optimum is −1.385510.
- (b) *Decay* (Thm 6.2, Cor 6.3). For every K ∈ 40Z:
  0 < F_K(ζ(5)) and log F_K(ζ(5)) ≤ U'K² + 24K log K + 200K.
  Hence 0 < Q_{40n,100000}(ζ(5)) < exp(−107.89n²) for large n. Fauzan's (2.8) gives exp(−79.07n²).
- (c) *Roots* (A, `theorems/T5_zeta5_exponent.md`, Thm 5.1, for any admissible nested-arcsine comparison measure σ):
  log‖B_K‖ ≤ K·max Φ + 5 log K + 69. This is a Christoffel-function bound at the poles for the roots of Q_K.
- (d) *Exponent, tier 1* (Cor 6.4). Under (F), for every μ > 103.3555 there is b₀ with |ζ(5) − a/b| > b^{−μ} for
  all b ≥ b₀. So μ(ζ(5)) ≤ 103.36. Fauzan's Corollary 1.2 gives 260.
- (e) *Exponent, tier 2* (Cor 6.5, with A `z3_arithmetic/theorems/T1_inner_extension.md` and `T2_parity.md`). Assume
  in addition Fauzan's Lemmas 3.1, 3.2, Propositions 4.3, 5.2, (5.7) and §5.2 on the extended ranges used by Z3:
  - inner range K/M < p with K/p > 61/37;
  - parity blocks for K/p ∈ [5/2, 6].
  Then the decay is exp(−151.92n²) and μ(ζ(5)) ≤ 73.41.

**Status.**
- (a): interval certificate, python-flint/arb (author; reproduced byte for byte in the author's clean rerun).
- (b)–(d): author proofs conditional on (F). The proof of (b) re-uses Fauzan's §6 argument with properties (i)–(v)
  of the measure, which the author re-read against the paper.
- (e): additionally conditional on Z3's T1 and T2. These extend Fauzan's local lemmas beyond the ranges he proved
  and beyond the Lean formalizations.
- Reviews: one internal review by agent Z2 in the same cloud Claude session (A `z2_fauzan/reviews/R2_Z1_T5_T6.md`):
  - PASS for T5, for Thm 6.1 (constants re-evaluated independently) and for Cor 6.4 conditional on (F);
  - OPEN for Cor 6.5, pending Z3's T2.
  Z2's `R1_Z3_T1_inner_extension.md` gives PASS for Z3's T1.
- Not reviewed outside that session, and not rechecked by rev4.
- (F) itself is Fauzan's claimed theorem. Per A's reading of 19 Oct, Romik's Lean formalization covers Theorems 1.1
  and 2.1 for every M ≥ 40 (prime number theorem as an axiom). Fauzan's own Appendix C (μ ≤ 260) is not formalized.

**Novelty.**
- No improvement of Fauzan's μ ≤ 260 was found: A's searches S8, S9, S11, S12 (19 Oct, `notes/00_literature_log.md`),
  and this consolidation's searches on 26 Oct (MathWorld still lists 260).
- New, as an improvement of the constants of an unrefereed preprint. It inherits that preprint's status.

**Superseded and left out.** Z2 T4 (μ ≤ 170.004), T5 Cor 5.2 (μ ≤ 141.89), Z3's μ ≤ 91.19.

**Limits.**
- ρ' is within 5.1·10⁻⁴ of the numerical optimum.
- The measured root growth ρ_K → ≈ 7.47 is within 1% of the certified 7.535. So, at this decay, the Christoffel
  argument cannot give much below ≈ 102 (A, T6, "Numerical cross-check").
- Z2's predicted exponents ≈ 61–76 are not bounds.
- This is a ζ(5) result; its relevance to ζ(7) is only through M-Z01.

---

## R-Z02. Vertex-permutation purity of cube integrals, and a pure Gram family of rank Θ(n⁶)

**Statement** (B, `round2/growing_purity/purity_proof.md`).
Setting:
- Fix d ≥ 2, X = x₁⋯x_d and I_d(F) = ∫_{[0,1]^d} F(x)/(1 − X) dx.
- For a monomial x^α let Z_j(x^α), for 2 ≤ j ≤ d, be the sum over poles of the (t + a)^{−j} coefficients in the
  partial fractions of ∏_i (t + α_i)^{−1}. Then I_d(x^α) = Σ_j Z_j(x^α)ζ(j) + (an explicit rational number). This is a
  formal identity: no independence of zeta values is used.
- Z_j is invariant under adding a constant to all α_i. It is therefore a functional on the group algebra of
  Z^d/Z(1,…,1).
- Write x_i = y_i/y_{i+1}, cyclically.

Claims:
- (a) If a Laurent polynomial F represents a class invariant under all permutations of the vertex variables
  y₁, …, y_d, then Z_j(F) = 0 for 2 ≤ j < d. Hence every polynomial lift F̃ satisfies I_d(F̃) ∈ Q + Qζ(d).
- (b) Take d = 7 and strictly increasing a ∈ Z⁷ ∩ [−n, n]⁷ with Σa_i = 0. Let P_a be the canonical lift of the
  alternant det[y_i^{a_j}]; its coordinate degrees are at most 3n − 3 (n ≥ 3). Then
  G_n = [I₇(P_aP_b)] = R_n + ζ(7)B_n, where:
  - B_n = −7! times the permutation matrix of a ↦ a* = (−a₇, …, −a₁);
  - G_n is positive definite;
  - the size is h_n = (841/129600)n⁶ + O(n⁵).
- (c) L_{6n}^7 R_n is integral. The crude criterion L_{6n}^7(det G_n)^{1/h_n} → 0 is impossible (R-Z03(a)).

**Status.**
- Author proof (B thread 01a15a1b-e1c6-7303-9e44-8c1e0a416064).
- Two Codex reviews, individual and cohort, both PASS/NARROWER. They recomputed the stated coefficient vectors and
  checked dimensions 3 and 5; the individual reviewer also checked 2, 4, 6 and 8.
- Consolidation check (non-Codex, independent code; `experiments/consolidation_2026-10-26/`). Exact verification for:
  - d = 3: every pair with n = 4 and n = 6 (36 + 171 pairs);
  - d = 5: n = 3 and n = 4 (6 + 78 pairs).
  In every case Z_j = 0 for 2 ≤ j < d and Z_d = sgn(reversal)·d!·[b = a*]. A negative control shows the checker
  detects impure inputs.
- rev4 did not check R-Z02. (It checked the five-dimensional space of M-Z15.)

**Novelty.**
- No source with this mechanism was found (queries in the consolidation notes). The closest known mechanisms are
  Vasilyev's nested integrals (rational multiples of ζ(n) at N = 0), the Rhin–Viola and Fischler symmetry groups of
  Beukers/Sorokin integrals, and Cresson–Fischler–Rivoal's symmetry phenomena (J. reine angew. Math. 617, 2008).
  These give parity or hypergeometric symmetry, not vertex-permutation averaging.
- Plausibly new; not confirmed.

**Relations.**
- The rank grows with n. This escapes the fixed-rank limit of the coordinate-symmetric construction (M-Z15).
- The rank attains the order of A's lattice-class bound (A `zx_wall` Lemma L3).
- B's README states that the monomial-only rank lemma of A `zx_cellular` (claim 11) does not apply here, because
  the cancellation is between monomials.

---

## R-Z03. The canonical alternant family: normalization lemma, arithmetic and real rates, two exclusions

**Statement** (B; all files under `round2/`).
- (a) *Gram domination* (`automaton_cover/alternant_bounds/GRAM_NORMALIZATION.md`). Let P₁, …, P_h be linearly
  independent real polynomials, G = [∫P_iP_j/(1 − X)] and U = [∫P_iP_j] over [0,1]⁷. If D > 0 and DU is integral,
  then D^h det G > D^h det U ≥ 1.
  - Consequence: no denominator that also clears the unweighted Gram matrix can make (det G)^{1/h} small;
    L_{6n}^7 is such a denominator for R-Z02.
  - With purity, any pencil normalization needs integer content larger than det(DU).
- (b) *Prime strata* (`PRIME_STRATA.md`). For primes p with p² > 6n:
  - p > 3n gives v_p(R_ab) ≥ −5, and p > 4n gives v_p(R_ab) ≥ −3. Hence
    D_n* = L_{6n}^7 / (∏_{3n<p≤4n} p² · ∏_{4n<p≤6n} p⁴) clears R_n.
  - For p > 2n, p⁷R ≡ C_p(a)·[b = a*] (mod p) (top-pole matrix).
  - For 2n < p ≤ 3n the content is divisible by p^{h−ρ_p}.
- (c) *Rates* (`CONTENT_THRESHOLDS.md`, `DENSITY_BOUND.md`, `ELEMENTARY_RATE.md`). Let C_n = (D_n*)^{h_n}/M_n, with
  M_n the proved content divisor. Then Q_n(T) = C_n det(R_n + TB_n) ∈ Z[T], and
  - limsup log C_n/(n h_n) ≤ 435978278/18396875 = 23.69850 using PNT;
  - limsup log C_n/(n h_n) ≤ (6A/5)·23.69850 = 26.19988 elementarily, with A = log2/2 + log3/3 + log5/5 − log30/30.
    The Chebyshev bound ψ(x) ≤ (6A/5)x + O(log²x) is proved in the file.
- (d) *Real lower bound* (`real_determinant/CANONICAL_OBSTRUCTION.md`).
  liminf log det G_n/(n h_n) ≥ −(175314/5887) log 2 > −20.642.
  The proof bounds det G_n from below via shifted Legendre polynomials, the triangular minor of leading monomials and
  Cauchy–Binet, and uses the exact mean degree 87657/11774.
  - So the full canonical family cannot satisfy a real bound of rate 23.7 or 26.2.
  - It could succeed only with an arithmetic normalization of cost below 20.642, which is not excluded.
- (e) *Exclusions* (`real_determinant/ADAPTATIONS.md`, `automaton_cover/alternant_bounds/CORNER_PRIME_FLOOR.md`,
  `UNIFORM_SHIFT_PRIME_FLOOR.md`). In two subfamilies, every integer-polynomial normalization of det G grows
  exponentially at ζ(7), because the proved arithmetic floors exceed the proved real lower bounds:
  - the thin corner family with δ = 1/1000;
  - the uniform shifts X^mP_a with m/n → τ, for every fixed τ ≥ 16.
  Moderate δ, moderate τ, other lifts and other index sets are open.
- (f) *Finite instances* (`root_alternant/`):
  - for an explicit P, ∫P²/(1 − X) = r − 5040ζ(7) with r an explicit rational number; the integral is about
    1.2753·10⁻¹¹, but its cleared form is about 1.53·10³⁷;
  - a fourth-power Pfaffian instance has sole zeta coefficient 681080400, with cleared value about 5.8·10⁸⁹.
  Neither is a small integer form.

**Status.**
Author proofs (B threads 01a15a1a… and 01a15a46…). Reviews, all by threads of the same Codex run:
- (a): the cohort review only (`reviews/cohort/gram_bounds`, PASS).
- (b): two reviews (`reviews/new_claims`, PASS; `reviews/determinant_bounds`, PASS).
- (c), (d), (e): one review (`reviews/determinant_bounds`, PASS); DENSITY_BOUND is marked conditional on PNT.
- (f): the first integral has two reviews (individual and cohort, PASS; all 116 numerator terms recomputed); the
  second has one (`new_claims`, PASS; the Pfaffian normalization was derived by the reviewer).
- `real_determinant/PROPOSED_EXTENSION.md` is an unreviewed lead, not a result.
- rev4 rechecked none of this.

**Novelty.**
- (a) is elementary. The general principle, that a fixed positive measure's Hankel determinants are limited by
  arithmetic, is known in spirit (Zudilin's transfinite-diameter criterion; cited by A `z4_design` N5). The exact
  domination statement was not found.
- (b)–(f) are new computations for a new family.

---

## R-Z04. Morris/Rodrigues positive pure integrals with explicit nonzero ζ(7) coefficient

**Statement** (B `round2/growing_purity/morris_rodrigues.md`, `morris_norm_bound.md`, `continuation_notes.md`;
`round2/arithmetic_lattice/RODRIGUES_CONTENT.md`, `CONTENT_PRIME_CUTOFF.md`, `DENOMINATOR_STATUS.md`;
`round2/root_alternant/PRIMITIVE_NORMALIZATION.md`; decay from C `zeta_diagonal/report_decay.md`).
Setting:
- F = ∏_{l=1}^{3} ∏_{i mod 7} (1 − x_i⋯x_{i+l−1}). Its class modulo X = 1 is the normalized Vandermonde R(y).
- D_i = x_i∂_i − x_{i−1}∂_{i−1} and T_n = ∏_i binom(D_i, n).

Claims:
- (a) Dyson's constant term, and Morris's constant term including negative b, are proved by elementary recurrences.
  - The coefficient of A_b in R^{2n+1} is (−1)^n L_{n,c}/7!, where L_{n,c} = ∏_{l=1}^{6}
    (lk − 1)!((l + 1)k)!/[(lk − c − 1)!(lk + c)!k!], k = n + 1 and b = (−c−3, −c−2, −c−1, −c, −c+1, −c+2, 6c+3).
- (b) T_n preserves integer polynomials. For exact lifts, I(P_aT_nQ) = λ_aI(P_aQ) with λ_a = ∏_i binom(−a_i, n);
  this is integration by parts with boundary vanishing, for Q = F^{2n+1}.
- (c) For c ≥ 3, n = 6c + 3, with the explicit common-shift lift P_a > 0 and Q = F^{2n+1}:
  - I(P_aT_nQ) = λ_aI(P_aQ) > 0, and this number is a rational plus the positive integer λ_aL_{n,c} times ζ(7);
  - for c ≥ 4, 0 < I(P_aQ) ≤ 5040ζ(7)ρ^{c−3}, with ρ = ∏_{j=1}^{6} j^j·12^{12}/(j + 12)^{j+12} ∈ (0, 1) and
    −(1/6) log ρ = 7.92105….
- (d) *Even family*. Take n = 6c and M = (x₂x₃²x₄³x₅⁴x₆⁵x₇⁶)^c. Then J_c = ∫MF^{2n}/(1 − X) = A_c/B_c + L_cζ(7),
  with L_c > 0 given by Morris.
  - The primitive integer form is Z_c = S_cJ_c, where S_c = B_c/gcd(A_c, L_c). It is unchanged by every nonzero
    rational rescaling or content normalization.
  - The corrected lemma needs (A, L) ≠ (0, 0): the unrestricted wording was REJECTED in review and corrected in
    `PRIMITIVE_NORMALIZATION.md`.
  - For p > 7n, p⁷A_c/B_c ≡ −Σ_{j≥p} h_j (mod p), where the h_j are diagonal coefficients.
- (e) Every prime divisor of the content of T_nF^m, for m ∈ {2n, 2n + 1}, is at most m. Universal and Frobenius
  divisors, and a nondivisor at a neighbouring prime, are also proved.
- (f) From C:
  - log J_c/(6c) → −α, with 8.4080693218 ≤ α ≤ 8.4080693441, by an exact rational certificate
    (`zeta_diagonal/report_decay.md`);
  - h_j = h_{13n−j}, and the large-prime residues are nonzero at (c, p) = (1, 67), (1, 71), (2, 139), so universal
    large-prime cancellation fails (`zeta_diagonal/report.md`).

**Status.**
- (a)–(e): author proofs. Reviews, all by Codex threads of one run:
  - `reviews/positive_integral`: PASS for (a)–(c) and the content lemmas;
  - `reviews/latest_cohort/morris_rodrigues.*`: PASS for (a)–(c);
  - `reviews/determinant_bounds/addendum_*`: PASS for (d) and the content cutoff;
  - `reviews/positive_integral/addendum_*`: REJECT, fixable, for the unrestricted normalization wording;
    `correction_*`: PASS for the corrected lemma.
- (f): author and captain only.
- rev4 rechecked none of this.

**Missing.** Obligation O3 (README): limsup log S_c / n < α. The generic bound is e^{91n}. Not proved.

**Novelty.** Morris's identity is classical, and equivalent to Selberg's integral. Its use to produce a positive
pure ζ(7) cube integral with an explicit nonzero target was not found. New construction; the arithmetic it would need
is open.

---

## R-Z05. The Ball–Rivoal V(7,1) cubic-filter pure lattice: rates, a sufficient criterion, and an equivalence

**Setting** (C, `zeta_astra/report_quotient.md`, `report_transition.md`).
- W_n(t) = (2t + n)·n!⁶·∏_{k=1}^{n}(t − k)(t + n + k) / ∏_{j=0}^{n}(t + j)⁸, and z = t(t + n).
- For a cubic P: Σ_{t≥1} W_n(t)P(z) = V + C3·ζ(3) + C5·ζ(5) + C7·ζ(7), in formal partial-fraction coordinates.
- M_n is the 4×4 coordinate matrix: rows V, C3, C5, C7; columns 1, z, z², z³.
- Λ_n = ker(C3, C5) ∩ Z⁴ is saturated of rank 2. The norm is ‖p‖² = Σ_j ((3n²)^j p_j)².
- D_n is the least common denominator of (V, C7) on Λ_n; d_n is the covolume and λ_{1,n} ≤ λ_{2,n} the successive
  minima.

**Statement.**
- (a) M_{n+1} = M_nT_n, with an explicit integer-polynomial certificate, and
  det M_n = (−1)^n·(8/(3n))·(3n)!/(n!)³ (`report_transition.md`). The limiting characteristic polynomial
  a⁴ − 828a³ − 132246a² + 260604a − 27 is prior art: Zudilin, J. Théor. Nombres Bordeaux 15 (2003) 593–626, §7, per
  C; the abstract confirms "low-order recursions for linear forms in odd zeta values".
- (b) *Analytic bound* (`report_norm_and_tail.md` §3, elementary). For n ≥ 8 and every integer cubic pure P:
  |V(P) + C7(P)ζ(7)| ≤ K·n⁴⁰·e^{−σn}·‖P‖_n, with K = 3¹⁷·2¹⁶ + 41.
  Here σ = −φ(τ) = 9.17486784697…, where τ is the critical point of
  φ(u) = 9u log u + (u + 2) log(u + 2) − (u − 1) log(u − 1) − 9(u + 1) log(u + 1). Equivalently σ = −log of the
  smallest characteristic root.
- (c) *Denominators* (`report_denominator_cutoff.md`, `report_determinant_cone.md`). (−1)^{n+1}Δ_n > 0 for all n, by
  an exact exterior-cone certificate, and limsup log D_n/n ≤ c∞ = 5 − π/(2√3) − (3/2) log 3 = 2.44518188488….
  - The bound imports PNT through the stored proof of Selberg (Ann. of Math. 50, 1949).
  - `root/coupled_image/report.md` further claims that the limit exists and equals c∞.
- (d) *Covolume* (`root/exterior_frobenius/report.md`). lim log d_n/n = μ = log ρ + γ₂ = 13.0634881489…, with
  ρ ∈ (133854, 133855) a root of X⁴ − 131976X³ − 251413578X² + 3563352X + 729 and γ₂ = 2 + π/(2√3) − (3/2) log 3.
  The proof uses the finite-field residue lemma of `root/owner_integrals/report.md` §12.
- (e) *Sufficient criterion* (`report_norm_and_tail.md` §4). If D_nλ_{2,n} ≤ e^{(σ−ε)n} on an unbounded sequence of
  n, then ζ(7) ∉ Q.
  - With (c) and (d) this follows if log(λ2/λ1)/n < 2σ − 2c∞ − μ = 0.395883775… along such a sequence.
  - `root/fermat_descent/report.md`: under rationality the shape rate would equal this threshold exactly.
  - With filter degrees 2n ≤ d ≤ 3n + 2 the threshold rises to 0.841065660126… (`root/full_degree_prime_transfer`).
- (f) *Equivalence* (`root/horizontal_symmetry/shortest_switch_report.md`). For any fixed choice of shortest vectors:
  ζ(7) ∉ Q if and only if the projective image line of the shortest vector changes for infinitely many n.
  - A change at n forces λ2/λ1 = O(n¹⁶). So, given (b)–(d), the criterion (e) is EQUIVALENT to ζ(7) ∉ Q.
  - It is a reformulation, not a reduction.
- (g) *Scope* (`zeta_astra/report.md` §4). With unrestricted signed filters every coprime integer pair is already
  realized at n = 1, degree 3 (explicit telescoping identities). Unrestricted small values are therefore ordinary
  rational approximation.
- (h) *Finite data* (`report_transition.md`, `root/horizontal_symmetry/switch_checks.json`):
  - at all 16 tested n ≤ 1024, both Gauss-reduced basis vectors give nonzero primitive forms of absolute value below 1
    (rational interval certificates);
  - the shortest line switches at every one of the 127 steps n ≤ 128.

**Status.**
- Author proofs: Astra thread 01a16a72… for zeta_astra, captain thread 01a16916… for the root packages. Captain
  self-scrutiny only; no independent review of any proof.
- rev4 (independent code) recomputed:
  - the closed forms and values of σ (also as −log of the smallest root of the quartic), c∞ (two ways), γ₂ and μ;
  - σ − c∞ − μ/2 = 0.19794 and 2σ − 2c∞ − μ = 0.39588;
  - the whole lattice from the definitions for n = 1, 2, 4, …, 1024 and 400, 700: λ1, λ2, D_n and the primitive
    forms, all below 1 in absolute value, matching C. It also confirmed det M_n against the formula in (a).
- rev4 did not check the proofs of (b)–(f).

**Novelty.**
- The family and its recurrence are prior art: Ball–Rivoal, and Zudilin 2003.
- Non-explicit small forms with a separate nonvanishing problem are prior art: Fischler's Siegel-lemma method.
- The lattice rates (c), (d) and the criterion and equivalence (e), (f) were not found elsewhere.
- (f) shows the route is exactly as hard as the target.

**Corrections to imported text.**
- C's `LANE_HANDOFF.md` calls σ "established" and treats c∞ as settled. In fact (b) is an author-level proof and (c)
  is an unreviewed limsup bound; its exactness, claimed in coupled_image, is also unreviewed.
- A's `theorems/T8_apery_limit_map.md` §8.2 says "the pure ζ(7) forms do not even tend to zero". That holds only for
  the bounded-window (Casoratian) eliminations of its Lemma 8.1. It is false for lattice-reduced forms, by (h), which
  rev4 confirmed.

---

## R-Z06. A prime-marked modular criterion

**Statement** (C, `root/congruence_filter/report.md` §§3–6, `root/modular_period_kernel/report.md`; definitions in
`root/cusp_divisor/report.md`).
Setting:
- E = ζ(7)A₀ − B₀ is the level-12 Fricke/Ramanujan period. A₀ = (hQ)⁶ is an integral unit series; B₀ = A₀f, where
  f_p = (1728/209)(1 + p⁻⁷) for primes p ∉ {2, 3, 11, 19}.
- For odd p put λ = (p − 1)/2, d = ⌊(p − 1)/4⌋ and S(T) = Σ_{j≤d} s_jT^j. Let α = [T^p]A₀(1 + T)^λS and
  β = [T^p]B₀(1 + T)^λS.

Claims:
- (a) For each fixed c < 4: |αζ(7) − β| ≤ C_c(d + 1)‖s‖∞·e^{−η(c)p}, where η(c) = (1/2) log c − (1/4) log R_c and
  R_c = (c + √(c² + 4c))/2 (a lemniscate contour). η(c) increases to
  η* = log 2 − (1/4) log(2 + 2√2) = 0.299516988665073….
- (b) Let Λ_p = {s : Σc_js_j ≡ 0 (mod L_{p−1}⁷)} and Λ'_p = {s ∈ Λ_p : p | s₀}. For s ∈ Λ_p∖Λ'_p, P = 209p⁷β and
  Q = 209p⁷α are integers, and P ≡ 1728s₀ ≢ 0 (mod p).
- (c) *Sufficient criterion.* If, for an unbounded set of admissible primes, there is s_p ∈ Λ_p∖Λ'_p with
  limsup log‖s_p‖∞/p < η*, then ζ(7) ∉ Q.
  - Conversely, if ζ(7) = a/b with p ∤ b, then for each fixed c < 4 every marked vector has
    ‖s‖₂ ≥ e^{η(c)p}/(C_c·b·p⁷·√(d + 1)), by a dual-vector argument.
  - `modular_period_kernel` keeps the criterion after projecting onto the two-dimensional span of the period
    columns.

**Status.**
- Author proofs (captain thread), captain self-scrutiny only; no independent review.
- rev4 recomputed η(c) and its limit η* numerically.
- The certified finite forms carry no information:
  - for example 0 < |Qζ(7) − P| < 10⁻³³⁷ with a 356-digit Q at p = 509;
  - the convergents of the continued fraction of any irrational number give |qx − p| < 1/q, which is below
    10⁻³⁵⁵ for such sizes;
  - C's controls show equally small finite forms at weights 3, 5 and 7.

**Missing.** Obligation O5. No equivalence with irrationality is proved, unlike R-Z05(f).

**Novelty.** Modular-form irrationality proofs for ζ(2) and ζ(3) are known (Beukers, 1980s). That reference is from
memory and was not re-checked here. This weight-8 construction with a prime mark was not found elsewhere.

---

## R-Z07. The class law for Fauzan-family denominators, and its strictness

**Statement.**
- (a) (A, `z5_intercept/theorems/T1_excess_identity.md`, `T3_class_law_bound.md`.) Write the Fauzan-family pencil in
  nodal form Δ ∝ det(D + S).
  - Exactly: excess_p := v_p(den χ) − θ-price_p = v^G_p(det D) − v^G_p(det(D + S)).
  - Upper bound: for p > s + 1 with (H1) (p-integral moments),
    excess_p ≤ Σ_{small a} κ_p(a) − 2·#{reflection pairs} + Σ_{big b} max(0, v^G_p(D_b)).
    When p² > 2 max(S ∪ Z) this is the class law.
  - Equality holds if and only if the pair-reduced small-node block is nonsingular mod p (Prop. 2).
- (b) (B, `lane_z/intercept/T1_counterexample.md`.) Take (s, K, N, σ, p) = (7, 36, 4, 6, 19), S = {5, …, 36} and
  Ω = ∏_{j≤4}(t + j²)⁵.
  - The class-law bound is 47 and the actual excess is 46.
  - The 20×20 pair-reduced residual matrix mod 19 has rank exactly 19, with an explicit kernel vector, and
    det C(0) ≡ 57 (mod 19²).
  - So universal equality fails even in the standard family. Here N/K = 1/9; the ratio 3/40 is not settled.
  - The certificate also corrects T3 Remark (1): κ_b ≤ σ fails, since κ₂₁ = 7 > 6.

**Status.**
- (a): author proof, unreviewed. Equality was checked at every computed prime p > K/2 in Fauzan-family designs
  (K ≤ 160); there are 9 failures among 2080 random (design, prime) cases.
- (b): exact finite certificate (Python standard library). Two Codex reviews (individual PASS, cohort PASS). rev4
  CONFIRMED it by a different exact method.

**Consequence.**
- A's `z5_intercept` "item 8" (the family cannot reach ζ(7) even if all primes p ≤ K/3 were free) used the class law
  as a lower bound. That lower bound is unproved, and (b) shows it fails pointwise; the no-go remains a conjecture.
- Aggregate savings (a positive density of saving primes, of total size ~K²) are open.

**Novelty.** No match in A's recorded searches. Toledo measures true contents but states no such law.
