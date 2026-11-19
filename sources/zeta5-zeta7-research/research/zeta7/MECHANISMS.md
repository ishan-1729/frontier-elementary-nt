# ζ(7) lane: register of mechanisms already tried

Purpose: no prover should rediscover a mechanism that has already been tried, and no reviewer should miss one.
Every entry names the exact point where the mechanism stops and keeps the quantifiers of each exclusion. Entries are
context. A proof dependency must use a statement in RESULTS.md, or a lane draft that passes the checker.

Consolidated 2026-10-26 from three branches, deduplicated across them:
- A = cloud Claude lane Z, `origin/claude/math-prime-cross-domain-g1vwsn:cloud_agent_claude/14_recurrences_wz_apery/`;
- B = Codex cloud run, `cloud-agent/elementary-connections:cloud_agent/claude_continuation/`;
- C = Codex captain campaign, `codex/zeta7-pollock-research:cloud_agent/goal_20261023/`.
Paths are relative to those roots. "(imported)" means the file is also under `experiments/<A|B|C>_…/` with the same
relative path; "(branch)" means it is cited in the branch only.

Statuses:
- DEAD (why): the mechanism fails in the stated scope.
- OPEN (what is missing).
- REFORMULATION: the remaining hypothesis is equivalent to ζ(7) ∉ Q.
- METHOD: a working practice, not a mechanism.

Review levels are those in RESULTS.md. Nothing in this file is referee-checked.

Fields: Id and name; Idea; Origin; How far it got; Where it stops; Status; Prior art and overlaps; Last checked.

## Family I: one-dimensional positive kernels (Fauzan-type Hankel determinants)

## M-Z01 Fauzan's ζ(5) Hankel scheme carried to ζ(7)
- Idea: replace Fauzan's weight by the positive ζ(7) Lambert weight w₇(y) = 2y⁷f⁽⁶⁾(y)/6!, f = 1/(e^{2πy} − 1).
  - Moments: (−1)^e B_{2e+2}(2e+3)⋯(2e+7)/720.
  - Pole values: j⁶(X − H_j^{(7)}) + 1/(2j) − 1/6.
  - Same Hankel matrices G_K(X) = [μ_X(D_N^r t^{i+j}/D_K)], with r, N/K and h/K free.
- Origin:
  - A (imported): `notes/09_zeta7_obstacle.md`; `z2_fauzan/` (T1 functional, T2 balance identities, T3 ζ(7)
    decision, notes N1–N7); `z3_arithmetic/` (exact polynomials and true contents, `notes/N02`–`N04`);
    `z4_design/notes/N01`; `zx_analytic/` A2, A3, A11.
  - C (imported): `root/cross_weight_controls/` (rates at K = 40: −0.994536, −0.165706, +0.766950 at ζ(3), ζ(5),
    ζ(7), equal to A's numbers).
  - C (branch): `root/fauzan_functional/`, `root/hankel_local_sharpness/`.
- How far it got:
  - The functional for every odd s, with moments, pole values, positivity and the pullback (A Z2 T1 and Z3 N01:
    two independent derivations; KNOWN in substance).
  - Exact balance identities for all (s, e, α, λ) (A Z2 T2, author proof): the K² log K terms cancel; the real
    constant does not depend on s; ∂R/∂s = λx.
  - Exact polynomials F_K with complete contents: ζ(5) for K ≤ 240, ζ(7) for K ≤ 160 (A Z3).
  - Fauzan's local estimates are sharp to leading order (measured): the outer bound equals the true valuation for
    p > K/2 (checked on 209 pairs, K = 40…240); the inner bounds are off by O(1) per prime.
  - Raising the weight never improves the real value (A ZX-A A2, author proof): Δ_h(ν_{s'}) ≤ C(s',s)^h·Δ_h(ν_s)
    for odd 3 ≤ s' < s.
  - The naive price of one unit of weight is at least (1/2 − o(1))h² (A ZX-A A11, author proof conditional on PNT).
    That this price is actually paid is only measured.
- Where it stops:
  - Every measured family has a positive ζ(7) total.
    - Fauzan's shape: +0.767, +0.846, +0.855, +0.862 per K² (K = 40, 80, 120, 160; exact contents).
    - Over 280 shapes (r = 2…10, N/K ≤ 0.3, h/K ∈ [0.25, 0.95], exact) the minimum is +0.555, and total/λ ≥ 0.83.
    - B₇ ≥ +0.868 on the scanned region λ ≥ 1 − α (model).
  - The critical weight is s* ≈ 5.2–5.4 (finite K).
  - Measured mechanism:
    - each prime p < K/3 costs about 2h extra powers, from the p^{−(s−1)} distribution loss on all rows;
    - each pole costs s + 2 powers for p > K/3.
  - There is no obstruction theorem. It would need an upper bound for the content, i.e. the exact p-adic valuation
    of Δ_K for every K, and no method for that exists (A `notes/09` §4).
- Status: DEAD as a direct transfer, as a finite-K measurement. "Positive for all parameters" is a conjecture.
- Prior art and overlaps:
  - C. Toledo's report (cert-machine, pinned 2026-10-05; machine-derived): 323 families, crossing near weight 5.3.
  - F. Calegari's Persiflage post of 24 Sep 2026: the method gives ζ(k) for k ≤ 5 "but not (at least directly)
    for larger k".
  - P. Anand's claimed ζ(7) proof by this route (Zenodo 10.5281/zenodo.22920911) fails: Toledo's audit finds
    A₂₀₀ + U ≈ +4.86, and C `fauzan_functional` notes that its normalization table contradicts its own definition.
  - Do not cite A's "the wall is mapped exactly" (RESUME.md) as a theorem.
- Last checked: consolidation 2026-10-26.

## M-Z02 The design space of one-dimensional kernels and residue-sum functionals
- Idea: change the kernel (y^a f^{(k)}, levels q, half-integer or shifted poles, zeros, windows, blocks, multiple
  orthogonality, residue sums) so that the ζ(7) balance becomes negative.
- Origin:
  - A (imported): `z4_design/` T1–T4 and N1–N6; `z3_arithmetic/notes/N05` (Lemmas M, R, F); `zx_algebraic/` B4;
    `zx_analytic/` A22, A24, A26 and the dead-bridge ledger `IDEAS.md`; `z5_intercept/theorems/T2_mixed_lattice.md`.
  - C (branch): `root/cyclotomic_poles/`, `root/weighted_zero_blocks/`.
- How far it got: author proofs of design restrictions.
  - Kernels y^a f^{(k)} (A Z4 T1). Pole values b^{s−1}ζ(s, b) + rational lie in Q + Qζ(s) for b ∈ (1/2)Z. Other
    rational b bring non-principal L(s, χ) and are excluded as designs. A level q ties the confinement constant 2π to
    the harmonic index of the nodes ("index lock").
  - One-sided residue-sum pencils without polynomial parts have deg_X Δ ≤ s (A Z4 T2). Their X-rank is at most
    (s + 1)/2 (A Z3 Lemma R).
  - Every block, multiple or vector Gram determinant with values in Q + Qζ(s) reduces to a Gram determinant of
    polynomial subspaces for one weight (A Z4 T4, found independently as Z3 Lemma M). For a fixed node set, the
    design acts only through the coupling E (Lemma F).
  - A pole of order r brings in ζ(s + r − 1) (A ZX-A A24). Purity on a Gram space therefore needs simple poles.
  - Parity-pure residue-sum functionals force the zero-block ratio N/K ≤ 1/3 (F2) or ≤ 3/10 (F3), by a degree count
    (A ZX-B B4.2). This excludes the Ball–Rivoal regime (ratio above e).
  - A fixed shift of the Hurwitz poles changes the raw determinant by only e^{O(K log K)}, so there is no K² real
    gain (C `cyclotomic_poles`).
  - Basis changes inside the same space leave primitive determinants unchanged (C `weighted_zero_blocks`).
  - Mixed lattices Z ∪ (Z + 1/2) remove the Bernoulli coupling (intercept +0.15) but double the slope (s* ≈ 2.8)
    (A Z5 T2).
- Where it stops:
  - The best measured ζ(7) total in Hurwitz-node designs is +0.744 per K² (two-block numerator, K = 100, exact).
  - Every measured family crosses zero at s* ≤ 5.36.
  - A fixed-measure Hankel determinant with the natural ζ(7) measure (−log t)⁶/(720(1 − t)) has true content
    9.41h² against a real part of −1.51h² (h = 60) (A Z4 N5; the criterion is Zudilin's).
- Status: DEAD for Hurwitz-node (Lambert, Eisenstein, residue-sum) designs at the measured sizes. This is not a
  theorem about all designs.
- Last checked: consolidation 2026-10-26.

## M-Z03 Removing the weight-independent intercept: the class law
- Idea: the arithmetic exponent splits as (harmonic node price, proportional to s) + (excess E₀ ≈ 0.90–0.93K²).
  Remove E₀ while keeping Fauzan's decay.
- Origin:
  - A (imported): `z5_intercept/` (T1, T3, N01); `z3_arithmetic/notes/N05` §2; `zx_wall/` (target C3b).
  - B (imported): `lane_z/intercept/`.
  - C (imported): `root/inner_reflection_blocks/`. C (branch): `root/hankel_local_sharpness/`.
- How far it got:
  - R-Z07: the excess identity, the class-law upper bound with its equality criterion, and the strict 46 < 47
    instance.
  - E₀ comes from crowding of nodes and numerator zeros in residue classes. Von Staudt–Clausen accounts for only
    ≈ 0.015K² at s = 7.
  - For Fauzan's own shape the decay equals the weight-7 harmonic price to within 2.5·10⁻⁵K² (numerical
    equilibrium).
- Where it stops:
  - A family with Fauzan's decay, a slope of ≈ 0.54 per unit weight per h² and zero intercept would have
    w* ≈ 7.0–7.5 (target C3b, marginal). None is known.
  - Singular reduced blocks mod p at a positive density of outer primes would remove ≈ 0.40K². Every computed
    Fauzan-family design is nonsingular at every p > K/2, except B's single prime.
  - Isolated reflection-block coranks do not add up (C `inner_reflection_blocks`, K = 120, p = 29). Five 9×9 blocks
    have rank 8 each, but the full 111×111 ranks are 110, 110, 111 at s = 3, 5, 7.
  - When B ≡ 0 (mod p^s) and the largest Smith exponent of A is below s, the full Gauss valuation equals v_p(det A).
- Status: OPEN only as aggregate savings (a positive density of primes, of size ~K²). A Z5's no-go (item 8) is a
  conjecture: it used the class law as a lower bound, which is unproved and fails pointwise.
- Last checked: consolidation 2026-10-26.

## M-Z04 Sharpening Fauzan's ζ(5) constants
- Idea: better comparison measures (near-equilibrium layers), Christoffel functions at the poles for the roots,
  local denominator savings (inner extension, parity blocks).
- Origin: A (imported): `theorems/T5`, `T6`, `notes/07`, `z3_arithmetic/theorems/T1`, `T2`, `z2_fauzan/notes/03`.
- How far it got: R-Z01 (μ(ζ(5)) ≤ 103.36, or ≤ 73.41, conditional).
- Where it stops:
  - The decay is within 5.1·10⁻⁴ of the equilibrium optimum.
  - The Christoffel route is within ≈ 1% of its limit at this decay.
  - The remaining proved-versus-true inner gap is ≈ 0.007K² (measured, K ≤ 240).
- Status: OPEN only for constants; not a ζ(7) route. The arithmetic grows linearly in s (M-Z01).
- Last checked: consolidation 2026-10-26.

## M-Z05 Two-species (Cauchy/Loewner) bimoment determinants; Taylor-contact purity
- Idea: the two-species Cauchy determinant geometry, as read from the OpenAI Catalan manuscript, with purity by
  Taylor contact.
- Origin: A (imported): `zx_wall/` M2, `notes/03_catalan_and_roth.md` (C4, T1).
- How far it got: measured margins per n² are
  - −0.70 (w = 2), −0.875 (w = 3), −0.10, +0.63, +1.36 (w = 4, 5, 6; n = 30), +2.19 (w = 7; n = 60);
  - the intercept is ≈ 0, but the decay is ≈ 2.9 per n², against Fauzan's 4.06 per h².
- Where it stops:
  - The critical weight is ≈ 4.1. Catalan-style weights never help (136 shapes at w = 5).
  - Taylor-contact purity needs convolution-type impurity symbols, which are not verified for any kernel whose
    decay grows with the weight.
- Status: DEAD at the measured sizes (n ≤ 60). Taylor contact is an OPEN idea.
- Prior art and overlaps: the OpenAI manuscripts in `github.com/openai/math` ("Catalan's constant is irrational";
  "The irrationality exponent of π is 2", 24 Sep 2026) are AI-generated and unrefereed. They are read only by A
  `zx_wall/notes/01_literature.md`, and not verified here. A separate claimed Catalan proof (Z.-W. Sun,
  arXiv:2609.04176) has a reported gap (D. Wachs, arXiv:2609.22339).
- Last checked: consolidation 2026-10-26.

## M-Z06 Dimension: lifted Beukers measures and Brown's multi-parameter criterion
- Idea: decay grows with the number of integration variables (Brown's criterion), while purity might survive in
  lifted Beukers measures g_k(x₁⋯x_d)dx.
- Origin: A (imported): `zx_wall/theorems/L1_lifted_measures.md` (L1–L3), M1, `notes/02_dimension.md`;
  `zx_cellular/notes/N02_purity_is_the_wall.md` (claims 10–13).
- How far it got: author proofs.
  - L1: the moments are divided differences of the Hurwitz ζ function.
  - L2: pure ζ(7) Gram families exist exactly in d = 1, 2, 3 (for d = 3, self-reciprocal modules). For d ≥ 4 the
    reflection keeps ζ(5) and ζ(3); this is proved for the reflection mechanism only.
  - L3: the rank of the ζ(w)-part equals the number of lattice classes modulo the diagonal: (n + 1)^d − n^d,
    i.e. 2n + 1 for d = 2.
  - Brown's threshold is T(w, r) = e^{−w(2r+1)/(r(r+1))}. T(7, 7) = 0.153, and the criterion is automatic for
    r ≥ 10.
- Where it stops:
  - Margins are positive already at w = 3 in the measured pure lifted families; the arithmetic outgrows the decay by
    a power of n.
  - Purity stops at dimension 3 in this class, while Brown's automatic regime needs at least 10 parameters.
  - A positive, pure, Zariski-dense 7-parameter family exists (A ZX-D claim 13). Its denominators are
    d⁷_{n₀+7Σn_i}, so it fails Brown's (1.7): it would need Sup < e^{−64}.
- Status: OPEN for the combination of purity and max-type denominators (obligation O6). The tested classes are DEAD.
- Prior art: F. Brown, "Mellin transforms, transfinite diameter and rational approximations of integrals",
  arXiv:2604.20741 (22 Apr 2026).
- Last checked: consolidation 2026-10-26.

## M-Z07 Cellular integrals on M_{0,n} and motivic purity
- Idea: choose cells or forms whose motivic coaction leaves only ζ(7). Couple them with positivity.
- Origin: A (imported): `zx_cellular/` (claims 1–19, `notes/N01`–`N04`).
- How far it got:
  - The π_odd^{10} cellular family is the very-well-poised family F₉(b). Its only points with A₃ = A₅ = 0 are the
    tails 2(ζ(7) − H_n^{(7)}) (b₀ ≤ 10: 106,445 classes; b₀ ≤ 16 in the standard domain). This is a measurement.
  - Brown's three "tantalising" ζ(5) examples are two tails and one isolated point (b₀ ≤ 12).
  - Casoratian purification reaches worthiness γ ≈ 0.60 at ζ(7), with margin +6.4 per n. At ζ(5) it reaches
    0.80–0.84 (Brown–Zudilin: 0.778 and 0.86). Purifying higher-weight forms is worse.
  - Twist elimination (Fischler–Sprang–Zudilin type) needs at least as many zeros as poles (a counting lemma).
  - No positive single-valued Mellin family exists on M_{0,n}.
  - A true-weight engine (Dupont's spectral sequence with Keel's Chow groups, validated for n ≤ 8) finds no pure
    cellular motive at n = 8.
  - Claim 18: for 8π8∨, no basic cellular integral contains ζ(3). It follows from gr₆ = 0 (computed mod p) and
    Brown's Cor. 8.2. This is an author-level computational proof, and the agent's self-review was never written.
- Where it stops:
  - Purity is the wall at every computed size.
  - M_{0,10} (the ζ(7) case) is beyond the implementation.
  - The unfinished lead D11: search M̄_{0,8} for pairs (A, B) with gr^W = Q(0) ⊕ Q(−5).
- Status: OPEN only through D11. The scans are DEAD.
- Prior art: Brown, arXiv:1412.6508; Brown–Zudilin, arXiv:2210.03391.
- Last checked: consolidation 2026-10-26.

## M-Z08 Arithmetic holonomy, p-adic and Galois/motivic structure
- Idea: use structural information (p-adic interpolation, differential Galois groups, coaction, monodromy) as a
  proof mechanism.
- Origin:
  - A (imported): `zx_algebraic/` (B2, B3, B9, D1, holonomy notes); `zx_wall/notes/04_dessins_conformal.md`.
  - C (imported): `root/horizontal_symmetry/report.md` (Theorem A).
  - C (branch): `root/sinc_galois/`, `root/owner_integrals/` §§2–9, `root/modular_padic_extension/`,
    `root/spinor_quantum/`, `root/modular_eichler/`, `root/complex_beta/`, `zeta_astra/report_sinc_coaction.md`.
- How far it got:
  - Fauzan's distribution constant equals −p^sL_p(s, ω^{1−s}) (A ZX-B B2). The proof relies on Washington's
    Thm 5.11 as recalled, with no stored source. It is checked numerically to O(p^{s+6}).
  - The kernel decays exactly like e^{−2π|Im z|} (KNOWN).
  - Calegari–Dimitrov–Tang holonomy bounds would need conformal radius above e⁷ (literature reading).
  - Belyi pullbacks cannot raise the conformal size (at most 16 at a cusp) (A ZX-C D1, known in substance).
  - The V(7,1) recurrence has no rational horizontal endomorphism other than scalars, and no nondegenerate rational
    bilinear pairing with a scalar shift multiplier (C Theorem A).
  - The 128-dimensional spinor connection is irreducible (Frenkel–Gross B₇ via D₈ folding, a named import).
  - Differential Galois groups computed: G_a⋊μ₄ (complex beta) and G_a⁷⋊PSL₂ (level-12 Eichler).
  - A common p-adic period exists, whose defect equals the target coefficient.
- Where it stops:
  - Functional or motivic independence does not transfer to numerical values.
  - A numerical Galois coefficient extractor on Q + Qζ(7) is equivalent to irrationality
    (`zeta_astra/report_sinc_coaction.md`).
- Status: DEAD as proof mechanisms. These are structural facts only.
- Last checked: consolidation 2026-10-26.

## Family II: recurrences and Apéry-type linear forms

## M-Z09 Apéry-like binomial sums with three-term recurrences
- Origin: A (imported): `theorems/T8_apery_limit_map.md` §8.3, `scripts/s23_binomial_sum_map.py`,
  `scripts/s23b_no_three_term.py`, `scripts/s24_power_sum_limits.py` and evidence.
- How far it got:
  - Theorem 8.5, an exact rank computation mod 2⁶¹ − 1 (author, unreviewed). For a ≥ 1 and 5 ≤ a + b ≤ 8,
    Σ_k C(n,k)^aC(n+k,k)^b satisfies no recurrence of order ≤ 2 with polynomial coefficients of degree ≤ 28 valid for
    all n ≥ n₀, for any n₀ ≤ 40. Theorem 8.5′ is the same for the C(2k, k)^c extension (80 members).
  - Three-term recurrences occur exactly for (a, b) = (1,1), (3,0), (2,1), (4,0), (2,2).
  - The Apéry limits of ΣC(n,k)⁵ and ΣC(n,k)⁶ lie in Qζ(2) + Qζ(2)².
  - For a = 7, 8 there is no small relation with zeta products of weight ≤ a − 1 (PSLQ, 275–320 digits).
  - ΣC(n,k)⁸ has order 4 and minimal root 0.11775; it fails against d_n⁷ by 4.86.
- Where it stops: only the bounded box (order ≤ 4, degree ≤ 28) was searched.
- Status: DEAD in the stated box.
- Last checked: consolidation 2026-10-26.

## M-Z10 The derivative lift of Apéry's double-pole function
- Origin: C (imported): `root/apery_multiplicity_reset/`. A `theorems/T2_apery_inverse.md` is left out (classical).
- How far it got:
  - The classical ζ(3) proof is rebuilt as a control, PNT-free, using ψ(x) ≤ 1.10555x + O(log²x) and
    A_{n+1} ≥ 31A_n. rev4 confirmed it. Classical, so not kept as a result.
  - For every odd s ≥ 5 and every n ≥ 1, the boundary forcing term of the derivative lift is strictly negative,
    with an explicit formula (10). So the recurrence becomes inhomogeneous, and applying it kills the target
    coefficient. rev4 checked the weight-5 formula exactly for n ≤ 30.
  - Repairing zeros to order 2r while keeping double poles and O(t⁻²) decay forces K ≥ rn (a degree count).
  - The repaired forms fail at ζ(5) and ζ(7) (n ≤ 64; at n = 32, log₁₀ of the error is +67.986 and +207.128).
- Where it stops: as stated. This does not exclude other Apéry-like constructions.
- Status: DEAD for this lift, at every odd s ≥ 5.
- Last checked: consolidation 2026-10-26.

## M-Z11 V(7,1) with bounded-window (Casoratian) purification
- Origin:
  - A (imported): `theorems/T8_apery_limit_map.md` §§8.1–8.2, `scripts/s21_vwp_pure_forms.py`,
    `scripts/s22_vwp_recurrence.py`, evidence `s22_*`.
  - A (imported): `zx_cellular/` claims 7–8.
- How far it got:
  - Lemma 8.1, the price of purity (author proof): forms obtained by eliminating q values from an order-r recurrence
    over a bounded window decay at most like |λ₁λ_{r−q+1}⋯λ_r|.
  - V(7,1) has characteristic roots 0.00010361, 1.94682, −138.739, 964.792, with product −27. Its pure Casoratian
    forms grow like 27/1.94682 = 13.87ⁿ (measured e^{2.35n} at n = 127).
  - V(7,2), V(7,3) and V(9,1) are worse.
  - With all denominator savings, the true margin is +6.4 to +7.0 per n (n ≤ 22).
- Where it stops: a fixed-window elimination pays the two largest roots.
- Status: DEAD for bounded-window elimination.
- Corrections:
  - A's sentence "the pure ζ(7) forms do not even tend to zero" holds only for this elimination. Lattice-reduced
    forms are below 1 at every tested n ≤ 1024 (M-Z12, R-Z05(h), confirmed by rev4).
  - A's novelty note ("no ζ(7) analogue of Zudilin's recursion") is withdrawn. The recurrence and its characteristic
    polynomial are in Zudilin, J. Théor. Nombres Bordeaux 15 (2003) 593–626, §7 (per C's
    `zeta_astra/report_transition.md`). Lemma 8.1 is probably folklore (A).
- Last checked: consolidation 2026-10-26.

## M-Z12 The V(7,1) cubic-filter pure lattice and its refinements
- Idea: replace fixed eliminations by the integer lattice of pure filters, and control its second successive
  minimum.
- Origin:
  - C (imported): `zeta_astra/` reports `report`, `quotient`, `transition`, `norm_and_tail`, `prime_scaling` and
    `prime_scaling_stop`, `endpoint_translation`, `microblock_transfer`, `denominator_baseline`,
    `denominator_cutoff`, `determinant_cone`, `coaction_height`, `inhomogeneous_boundary`, `remainder_divisibility`,
    `second_minimum_stop`; and `root/` packages `exterior_frobenius`, `coupled_image`, `horizontal_symmetry`,
    `fermat_descent`, `full_degree_prime_transfer`, plus `owner_integrals` §12.
  - C (branch): `root/enlarged_cusp_arithmetic`, `two_level_image`, `quotient_duality`, `residue_tail_transfer`,
    `integral_weight_filtration`, `full_pure_image`, `old_prime_filtration`, `metric_comparison`,
    `weighted_pearson`, `zero_range`, `marked_purity`, `frobenius_excursion`, `integer_spectral`,
    `harmonic_operator`, `alternating_weyl`, `ghost_descent`, `ghost_return`, `owner_quotient_moments`,
    `hankel_moment_obstruction`.
- How far it got: R-Z05.
- Where it stops: the criterion is equivalent to the target (R-Z05(f)). Each refinement of the arithmetic image
  moves the threshold to exactly the value that rationality would force:
  - 0.395883775 for cubic filters (`fermat_descent`);
  - 0.841065660 for full-degree filters (`enlarged_cusp_arithmetic`);
  - a two-level window at a = 7/6 gains 19/60, with threshold σ + η − 521/60 (`two_level_image`).
  Other refinements:
  - quotient duality removes the lift charge and extends the decay to degree 3n + 2 (`quotient_duality`);
  - the upper-degree band from ≈ 2.7319n costs only polynomially (`residue_tail_transfer`).
- Exclusions, each in the scope its report states:
  - Unrestricted filters realize every rational pair already at n = 1 (`zeta_astra/report.md`).
  - The quartic-norm method misses the first-minimum budget by ≈ 3.77 (`report_norm_and_tail` §1).
  - The C7 coordinate-spacing argument gives only the kernel-vector bound (`report_coaction_height`).
  - The inhomogeneous C7 inverse column has height at least e^{σn}/(Kn⁴⁰) and exhausts the plane
    (`report_inhomogeneous_boundary`).
  - R_{n,h} ∈ (n!)^{2h−2}Z[y] settles divisibility but not λ₂ (`report_remainder_divisibility`).
  - There is no rational horizontal endomorphism other than scalars (`horizontal_symmetry`, Theorem A).
  - Every quadratic-unit ghost orbit has at most one pure return, for all signed exponents (`ghost_descent`).
  - A genus-5 obstruction (no point mod 7) excludes all nontrivial pairs of rational norm-one units in the two n = 1
    return directions (`ghost_return`).
  - Prime reduction sends (a, b) to (p⁷a, b), and the boundary term of the fixed-target descent is nonzero
    (`fermat_descent`).
  - Changing representatives costs at most √ν, with ν ≤ d − 2 for n ≥ 16 and d ≤ ⌊n/6⌋ + 4
    (`metric_comparison`).
  - A coefficient scale θ < e^{(σ−6)/3} pays an analytic penalty at d = 3n + 2 (`weighted_pearson`).
  - A negative 7×7 Hankel direction excludes positive-moment models for this error sequence
    (`hankel_moment_obstruction`).
  - For p > max(n, 7), the integer-polynomial pure image already equals the cubic one
    (`integral_weight_filtration`).
  - Row scaling cannot extend the map theorem, in the literal case n = 7, m = 0, p = 7 (`old_prime_filtration`).
- Status: REFORMULATION. The refinements are DEAD as stated. The only remaining opening is an independent
  arithmetic reason why the shortest image line cannot stabilize (C `shortest_switch_report.md` §4).
- Prior art: the family is Ball–Rivoal's; the recurrence is Zudilin's (2003); Fischler's Siegel-lemma approach.
  C stored a copy of Fischler's paper, which was not imported.
- Last checked: consolidation 2026-10-26.

## M-Z13 Modular (Ramanujan/Eichler) periods with prime-marked lattices
- Origin:
  - C (imported): `root/cusp_divisor/`, `root/congruence_filter/`, `root/modular_period_kernel/`.
  - C (branch): `root/modular_eichler/`, `modular_lagrange/`, `modular_marked_quotient/`, `modular_metric_lift/`,
    `modular_local_saturation/`, `modular_frobenius_kernel/`, `modular_frobenius_descent/`,
    `modular_padic_extension/`.
  - A (imported): `zx_analytic/` A0, A22; `zx_wall/notes/03` (the Eichler–Nesterenko "dream").
- How far it got:
  - R-Z06.
  - Five level-12 Ramanujan moments cancel π⁷ and keep 209ζ(7); the coefficient-ratio limit ζ(7) is proved
    (`modular_eichler`).
  - Integral Lagrange filters reach denominator 209·lcm(n, lcm(1,…,n−k−1))⁷, but the maximal filter is excluded
    along primes (`modular_lagrange`).
  - For every prime p ≥ 8s − 3, the integer zero-image first coordinates form exactly p^sZ
    (`modular_marked_quotient`).
  - Projected vectors lift with error O_s(p^s) (`modular_metric_lift`).
  - After Frobenius descent the image index has prime support {2, 3} at weights 3 and 5, and {2, 3, 5} at weight 7
    (`modular_frobenius_descent`).
- Where it stops:
  - The uniform marked minimum (obligation O5) is unproved.
  - Small finite forms occur at weights 3, 5 and 7 alike, so a finite success carries no information.
  - The certificates are weaker than continued-fraction approximations.
- Status: OPEN. R-Z06(c) is a sufficient criterion, not proved.
- Last checked: consolidation 2026-10-26.

## M-Z14 Multiple-harmonic-star, tensor, Stieltjes, quadrature and reflection approximations
- Origin:
  - C (imported): `zeta_astra/report_stieltjes_tensor.md`, `root/hahn_cauchy_descent/`.
  - C (branch): `root/tensor_extrapolation/`, `stieltjes_descent/`, `moment_quadrature/`, `fixed_prime_moments/`,
    `prime_fold_lifting/`, `reflection_pade/`, `reflection_fermat/`, `reflection_acceleration/`;
    `gram_arithmetic/` and `multipoint_bernoulli/` (superseded by `hahn_cauchy_descent` where r = o(N)).
- Exclusions (author or captain proofs, unreviewed):
  - Positive tensor regularization: S_N = (1/2)H_N^⋆(2,2,2,1) = a_N/q_N with q_N | 2·lcm(1,…,N)⁷,
    v_p(q_N) = 7 for every odd prime N/2 < p ≤ N, and 1/(N+1) < ζ(7) − S_N < 5/(4N). The primitive forms grow
    exponentially.
  - Ordinary Richardson extrapolation of order r = 1, …, 8: the fully reduced denominator rate is 7r (112 for the
    five-point log-aware filter); every fixed r ≥ 2 has rate at least 7.
  - Polynomial quadrature, for every r = o(N): the denominator and the primitive form both grow at rate at least
    7 log 2 per N; for 1 ≤ r ≤ N/32 the rate exceeds 1351/2400. This uses the all-r Cauchy–Hahn determinant formula
    proved in `hahn_cauchy_descent`.
  - Fermat prime arithmetic at a fixed base N gives lower rate 8/(N + 1) per order (`stieltjes_descent`).
  - Reflection series: v₃(q_{3^h−1}) = 7h for every h ≥ 1 (`reflection_fermat`).
  - The half-shift WZ series of Chen–Liu–Lu (arXiv:2610.16844, as cited by C; not checked here) has growing
    q·error for every fixed odd s ≥ 5 (`reflection_acceleration`).
- Where it stops: the fixed-base growing-order regime, larger order/base ratios and signed growing-order filters are
  open (finite failures only).
- Status: DEAD in the stated regimes.
- Last checked: consolidation 2026-10-26.

## Family III: seven-dimensional cube integrals ∫ P(x)/(1 − x₁⋯x₇) dx

## M-Z15 Coordinate-permutation purification (fixed rank 5)
- Origin: B (imported): `lane_z/purity/` and its two reviews.
- How far it got:
  - An explicit five-dimensional space E = 1 + X, P₁, P₂, Q₁, Q₂, built from H_ij(15, 7, −6) and H_ij(15, 3, 1). Its
    Gram matrix is diag(4ζ(7) − 385/128, (4440ζ(7) − 137895/32)K, (4440ζ(7) − 137895/32)K); all coefficients of
    ζ(2), …, ζ(6) vanish.
  - A rotation lemma, and the full signatures of the 64-dimensional Boolean forms: (36, 28, 0) for ζ(3) and
    (35, 29, 0) for ζ(5).
  - Two Codex reviews PASS; rev4 recomputed the Gram matrix.
- Where it stops:
  - In the whole Q[X]-module generated by the space, every Gram matrix has ζ(7)-coefficient rank at most 5, because
    Z_j(X^mR) = Z_j(R).
  - Damping by f(X) with f(1) = 0 kills the target.
  - The denominators are L⁷_{N+2}.
- Status: DEAD (fixed rank). Superseded as a source of pure spaces by R-Z02.
- Last checked: consolidation 2026-10-26.

## M-Z16 Vertex-permutation purity: the canonical alternant Gram family
- Origin: B (imported): `round2/growing_purity/` (purity_proof), `round2/root_alternant/`,
  `round2/automaton_cover/alternant_bounds/`, `round2/real_determinant/`, and their reviews.
- How far it got: R-Z02 and R-Z03.
- Where it stops:
  - Obligation O2: the full canonical family fails unless its true arithmetic cost is below 20.642·n·h_n.
  - The thin corner δ = 1/1000 and fixed shifts τ ≥ 16 are proved to fail.
  - `real_determinant/PROPOSED_EXTENSION.md` (unreviewed) suggests narrow-band prime floors for families near fixed
    templates.
- Status: OPEN for subfamilies, other lifts and better normalizations.
- Last checked: consolidation 2026-10-26.

## M-Z17 Single positive pure forms (Morris/Rodrigues) with multipliers, lifts and densities
- Origin:
  - B (imported): `round2/growing_purity/morris_*`, `continuation_notes.md`, `selberg_denominator_state.md`,
    `round2/arithmetic_lattice/`.
  - C (imported): `zeta_diagonal/` reports `report`, `report_decay`, `report_endpoint`, `general_line_report`,
    `report_variable_shift`; `zeta_lifts/` reports `report`, `report_density_cone_geometry`.
  - C (branch): the other `zeta_diagonal/` and `zeta_lifts/` reports; `root/plane_exhaustion_report.md`,
    `density_low_compactness_report.md`, `density_positivity_cost_report.md`, `prime_blocks_report.md`,
    `arithmetic_preflight_report.md`; `exponent_line_report.md` (superseded by `general_line_report.md`).
- How far it got: R-Z04. Exclusions (author or captain), in the scope stated:
  - Exponent lines: for every seven distinct integer directions with |a_i| ≤ m (shift and directions may vary with
    n) and every n ≥ 2⁴⁰, the logarithm of the positive primitive form is at least n/2 (`general_line_report`,
    `report_variable_shift`).
  - Fixed base endpoint: the constrained minimum is Θ(degree⁻⁹⁶), so small forms need B/(degree + 1)⁹⁶ → 0
    (`report_endpoint`).
  - The repeated-slope M5 split has log Z > d/10 for d ≥ 2³⁹ (`report_uniform_split_obstruction`).
  - All-vector metric floor: Z > e^{d/4}/G_eff for d ≥ 2¹¹⁰⁰ (`report_balanced_metric_obstruction`).
  - On the fixed planes d = 6 and 8, the minimum primitive positive integral over all integer representatives is
    exactly 1 (`plane_exhaustion`, two independent complete sieves over 264,666,820 candidates).
  - For fixed binomial (d, h) there is a positive primitive floor for every rational coefficient vector and every
    nonnegative strict-pure density (`density_low_compactness`, `report_density_dual_floor`). So a vanishing
    sequence needs h → ∞.
  - The triple-pole sparse family f = 32 − t^r + 32t^{2r} with h = 4r has log Z ≥ h for h ≥ 2⁴⁰
    (`report_triple_purity`).
  - Reflected sparse squares: every coefficient-one member with H ≥ 2⁴⁰ is excluded (`report_reflected_sparse`).
  - Every fixed limiting ratio h/d ≥ 5/3 is excluded in the positive mixture (`report_mixture_binomial`).
  - The unrestricted rational cone has primitive forms tending to 0 if and only if ζ(7) ∉ Q
    (`report_density_cone_geometry`; a scope diagnosis).
  - Universal large-prime cancellation fails at (c, p) = (1, 67), (1, 71), (2, 139) (`zeta_diagonal/report.md`).
  - Finite searches: every bounded, congruence-guided, isotropic, quadric, chart and CRT candidate fails after full
    primitive normalization (the other `zeta_lifts/` and `zeta_diagonal/` reports).
- Where it stops: obligation O3, limsup log S_c / n < α ≈ 8.40807 for the actual even family.
- Status: OPEN for the actual family. The listed subclasses are DEAD.
- Last checked: consolidation 2026-10-26.

## Family IV: other connections

## M-Z18 Gamma classes, quantum cohomology and spinor periods
- Origin: C (branch): `root/gamma_grassmannian/`, `root/spinor_quantum/`.
- How far it got:
  - For every N ≥ 13, a primitive degree-seven functional on Gr(5, N) pairs its Gamma class to
    N(N − 12)(5N − 38)ζ(7). This uses named Gamma and quantum-Pieri inputs.
  - Exact recurrences come from Schubert seeds.
  - Im L(Γ Ch(S*)) = 8cπ⁷/(45N) for every homogeneous formally pure class.
  - Pure functionals exist on OG(8, 16) and LG(7, 14); OG4 reproduces Apéry's ζ(3).
- Where it stops:
  - All retained ζ(7) primitive errors exceed 1 after full reduction: 30 on the Grassmannians, 24 + 88 in the spinor
    case.
  - Constant spectral splitters are excluded (Lie(A, B) = sl₁₃).
  - Frenkel–Gross irreducibility excludes rational horizontal projectors.
- Status: DEAD as tested.
- Prior art: Apéry constants of Grassmannians and Fano varieties (Golyshev; Galkin, arXiv:1604.04652), per A's
  searches S16–S17.
- Last checked: consolidation 2026-10-26.

## M-Z19 Shuffle/Pfaffian positive families, moment-pencil lattices and Stieltjes cones
- Origin: C (branch): `root/pfaffian_confluence/`, `moment_pencil_lattice/`, `primitive_pencil_cones/`,
  `bingo_mechanisms/`, `bingo_arithmetic_followup/`. These were inspired by Ishan's notes; the notes themselves are
  not imported.
- How far it got:
  - Positive pure ζ(7) Pfaffians with pair interaction (t_i − t_j)²((t_i − t_j)² + 4c), whose leading coefficient
    is an explicit product of quadratic norms.
  - For every rational c ≥ 0 at fixed (K, h), the primitive value grows with the height of the parameter (two-chart
    Bézout). The deformation parameter alone cannot give irrationality.
  - Complete rational moment pencils have an integral ζ-coefficient matrix and a nonzero determinant.
  - Equal-volume lattice changes are supported only at primes ≤ the dimension, at content primes and at squared
    discriminant primes.
  - Nested Stieltjes intervals J_K ∋ ζ(7) have width at most e^{−(2 + 2 log π)K + O(log K)}.
  - The natural half-order Euler bounds have reduced denominators of exact p-valuation p + 5, and fail along primes.
- Where it stops:
  - The unrestricted best primitive value tends to zero if and only if ζ(7) ∉ Q (`primitive_pencil_cones`).
  - Fixed-size families are excluded; growing-size arithmetic is open.
- Status: DEAD as tested. The unrestricted cone is a REFORMULATION.
- Last checked: consolidation 2026-10-26.

## M-Z20 Ishan's integral ideas and related exact kernels
- Origin: C (branch): `zeta_astra/report_sinc_integrals.md`, `report_sinc_deformations.md`,
  `report_sinc_coaction.md`; `root/sinc_galois/`, `owner_integrals/`, `owner_polylog/`, `owner_quotient_moments/`,
  `complex_beta/`, `mellin_resolvent/`, `rational_resolvent/`, `cellular_parity/`, `rational_star_probe/`,
  `cyclotomic_reflection/`. Ishan's originals are not imported.
- How far it got:
  - Ishan's sinc family satisfies B_n + iA_n = i^{n−1}2^{1−n}n!ζ(n) for n ≥ 2; for example B₇ = −315ζ(7)/4.
  - A complex-beta plane integral gives I_n = 2ζ(n) for odd n and 0 for even n.
  - Positive Hankel minors give ζ(3)ζ(7) − ζ(5)² > 9/128.
  - ∫F₆/x = 2ζ(7).
  - All-weight convolution identities for polylogarithms.
  - The exact reflection identity F₇(1/6) − 190F₇(1/2) − 720F₅''(1/2) − 180F₃''''(1/2) = 274320ζ(7)
    (`reflection_pade`).
- Exclusions:
  - Boundary-vanishing radial filters remove the whole target (`complex_beta`).
  - No nonzero sign-definite polynomial multiplier of the F₇ kernel removes ζ(2), in any degree (`owner_polylog`).
  - Positive polynomial purifiers of F₆/x need B = o(degree⁴) (`mellin_resolvent`).
  - Positive pure rational weights need B(a/(a + 1))² = o(degree⁴) (`rational_resolvent`).
  - A fixed nonconstant polynomial observable of the bipartite star has error at least c_f/k² (`cellular_parity`).
  - Positive fixed-observable moment closure forces φ to be constant (`owner_quotient_moments`).
  - A numerical Galois coefficient extractor is equivalent to irrationality (`report_sinc_coaction`).
- Where it stops: each construction either reduces to the original question or fails after primitive
  normalization.
- Status: DEAD as tested.
- Last checked: consolidation 2026-10-26.

## M-Z21 Matched odd-weight controls (method)
- Origin:
  - C (imported): `root/cross_weight_controls/`, `root/inner_reflection_blocks/comparison.md`.
  - C (branch): `root/modular_frobenius_descent/comparison.md`.
- Finding:
  - Each mechanism was applied at ζ(3), ζ(5) and ζ(7), with complete primitive arithmetic.
  - Fauzan's determinant separates 3 and 5 from 7. Its rates at K = 40 are −0.994536, −0.165706, +0.766950, and
    about 88% of the loss from 5 to 7 is arithmetic.
  - Ordinary level-12, trace, band and most fixed quantum forms fail at the known irrational benchmarks too.
  - Small finite modular forms occur at every weight.
  - OG4 reproduces the Apéry control.
- Use: before a ζ(7) mechanism is credited, it should close on a known control with its full arithmetic. The
  controls are ζ(3) (Apéry) and Fauzan's claimed ζ(5), which is unrefereed.
- Status: METHOD.
- Last checked: consolidation 2026-10-26.

## Left out of this register

Outside this lane's scope (recorded in triage.md):
- A's recurrence-to-sum programme for Ishan's MSE question 2965106: T1, T2, T7 and the log 2 / ζ(2) pipeline;
- all Pollock material in B and C;
- process files and logs;
- snapshots (B's `source_snapshot/`);
- third-party texts (`*/sources/`);
- Ishan's originals.
