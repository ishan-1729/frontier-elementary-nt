# 02 — Results, precision of claims, limitations and open gaps

**IMPORTANT:** These are *historically recorded* outcomes in the copied dossiers. This handoff has not rerun every Lean build, script or original-paper comparison. A label describes **the source's recorded status**, not a new global certification. Read [AGENTS.md](../AGENTS.md) before using any result as a lemma.

## A. Direct Diophantine and elementary arithmetic results

| Topic | Statement / exact scope | Recorded status | Local evidence |
| --- | --- | --- | --- |
| Ljunggren | For x,y positive integers, x²+1=2y⁴ iff (x,y)=(1,1) or (239,13). Separate Ghost-descent and integer norm-sign routes. | **LEAN / two proofs** reported; separate written audit | [theorem card](../sources/math-prime/research/proofs/theorems/ljunggren/README.md), [norm-sign manuscript](../sources/math-prime/research/proofs/theorems/ljunggren/norm_sign/elementary_proof.md), [audit receipt](../sources/math-prime/lean/FINAL_AUDIT.md), [Ghost audit](../sources/math-treasure/research/ghost-classification-audit/proof-audit.md) |
| Mordell k=7 | x²+7=y³ has integer solutions (x,y)=(±1,2),(±181,32). | **LEAN** reported | [card](../sources/math-prime/research/proofs/theorems/mordell_k7/README.md) |
| Cubic x²+x+1 | x²+x+1=y³ in integers has (0,1),(−1,1),(18,7),(−19,7). | **WRITTEN**, complete source argument | [proof](../sources/math-prime/research/proofs/theorems/x2_x_1_y3/proof.md) |
| Fixed Mordell families | For each fixed nonzero k, equation x²+k=y³ has finitely many integral solutions. Source splits positive k via binomial approximants, negative k via unit Ghosts. | **LEAN** finiteness entries reported, not full solution algorithm for each k | [positive k card](../sources/math-prime/research/proofs/theorems/mordell_positive_k/README.md), [negative k card](../sources/math-prime/research/proofs/theorems/mordell_negative_k/README.md) |
| Mordell 1963 | Integer y(y+1)=x(x+1)(x+2): x∈{0,−1,−2}, y∈{0,−1}; also (x,y)=(1,2),(1,−3),(5,14),(5,−15). | **WRITTEN-AUDITED** norm-triple proof; second Ghost proof written | [complete proof](../sources/math-prime/research/proofs/theorems/mordell_1963/proof.md) |
| Elementary quartic family | (t+1)x⁴−t y²=1 for positive integers. At t=2: (1,1),(3,11); t=3,7,11,26: only (1,1). At t=6,12 odd-Pell-index classifications only; even index and t=19 remain unclosed in this programme. | t=2 **LEAN**; t=3,7,11,26 **WRITTEN-AUDITED**; t=6,12 **PARTIAL**, t=19 **OPEN** | [exact card & boundaries](../sources/math-prime/research/proofs/theorems/quartic_t_family/README.md) |
| Three small Mordell–Thue closures | All integer solutions of X³=Y²−4 are (0,±2), of X³=Y²−3 are (1,±2), of X³=Y²−2 are (−1,±1). | **WRITTEN-AUDITED** local independent reviews, exact unit/3-adic proof | [full dossier](../sources/math-treasure/research/mordell-thue-completion/README.md) |
| Three figurate rank-two cases | Complete classification of integer points/finite Thue subsystems for triangular-tetrahedral, triangular-square-pyramidal and y²=x³+17 families, with a total 12+12+8 signature-(2,1) Thue equations in the method; see source catalogue for exact integer lists. | **WRITTEN** with in-session review; NOT implied LEAN or journal-certified | [main catalogue](../sources/math-prime/research/proofs/CATALOGUE.md), [unit method](../sources/math-prime/research/proofs/methods/ghost.unit.multi_unit.md) |
| Lucas perfect-power triples | For Lucas U_n(A,±1), fixed r≥2, source reports that U_n,U_{n+rk},U_{n+2rk} cannot all be r-th powers, and gives an obstruction showing every fixed modular sieve admits false positives. | **WRITTEN-AUDITED**, exact original hypotheses in card | [card](../sources/math-prime/research/proofs/theorems/lucas_power_triples/README.md) |
| Weighted Vieta and fixed offsets | Complete integer-orbit and finite-seed classifications for various integral-reflection conics, including positive offsets without a uniform small-offset assumption. | **WRITTEN-AUDITED** / exact scripts, depending on particular claim | [Vieta results](../sources/math-treasure/research/vieta-descent/README.md), [positive offsets](../sources/math-treasure/research/positive-offset-classification/proof.md) |
| Sum-and-square-sum | For d≥8, integer y_i with Σ y_i=s and Σ y_i²=N exist iff N≥balanced minimum M_d(s) and N≡s mod2. | **WRITTEN-AUDITED** recorded | [card](../sources/math-prime/research/proofs/theorems/sum_and_square_sum/README.md) |
| n² divides 2ⁿ+1 | For positive n, n²|2ⁿ+1 iff n∈{1,3}. | **LEAN** reported | [radical descent card](../sources/math-prime/research/proofs/theorems/radical_descent/README.md) |
| Binomial divisibility and harmonic congruences | Truncated-binomial valuation criteria, Wolstenholme congruences, Jacobsthal–Kazandzidis, some Legendre-Hankel determinant claims. | **WRITTEN, branch-unreviewed** in original catalogue; finite numerical checks are not independent full proof audit | [catalogue](../sources/math-prime/research/proofs/CATALOGUE.md) |
| Integral matrix cubes | Each 2×2 integer matrix is the sum of three integer matrix cubes; w₃(M₂(Z))=3. Local two-cube obstruction classified in source; integer two-cube question not resolved. | **LEAN** reported for three-cube theorem; two-cube global **OPEN** | [card](../sources/math-prime/research/proofs/theorems/matrix_cubes/README.md) |

### Example-specific caveats

- The preliminary [MathAgent example](../sources/MathAgent/knowledge/examples/vierzahlensatz/kieren_mse_x2_plus_1_eq_y3.md) for X²+1=Y³ contains an **intermediate reduction and method outline**, not a complete argument. Do not cite it as an independent proof.
- A prior soft evaluation purported to finish Ljunggren by asserting a hard Pell square-term exclusion; the original primer explicitly rejects that shortcut. The later norm-sign/Lean-audited proof is **separate**. [Primer](../sources/Ljunggren/Vierzahlensatz_Primer/README.md).
- For all source result cards, mathematical **truth**, claim **novelty/elementarity**, local **review status**, and worldwide **priority** are four distinct questions.

## B. Methods established as reusable statements, but with conditions

| Technique | Strongest useful proven/claimed scope | Remaining requirement |
| --- | --- | --- |
| Ghost rational orbit | Complete signed-index coverage and marker exclusion in selected rank-one rational Ghosts | Rank-two ghost lattices still lack a general all-index denominator marker theorem |
| Integer norm triples | Construct polynomial norm, adjugate, exact box division, Bézout and units in specified cubic orders | A new ring needs its *own* exact box/multiplier/unit certificate |
| Multi-unit Ghost coverage | Finitely many seeds for norm-c classes under explicit independent units; two-marker finite lifts on quartic signature (2,1) | A termination certificate for the sieve is instance-specific |
| One-index finite lifting | Exact valuation law for specified integral polynomial/Somos period maps | A suitable prime/period and nonzero primitive first increment are **not automatic** |
| Two-index zero isolation | Coprime reductions, saturation or frame theorem; Jacobian local zero isolation | Number of independent markers must be adequate; local p-adic roots can be false positives |
| Wronskian norm-sign ladder | Alternating integer/rational approximation to selected algebraic angle roots with explicit constants | Denominator normalization and approximation exponent can be too weak at small parameters or other sectors |
| Vieta fundamental arcs | Finite seed classification for prescribed integer conics and prescribed positivity regions | Not an unconditional closed form over all parameters; mixed signs require new branch analysis |

**Primary sources:** [method catalogue](../sources/math-prime/research/proofs/methods/README.md), [generalization map](../sources/math-prime/research/proofs/GENERALIZATIONS.md), [two-index full proof](../sources/math-prime/research/proofs/methods/lifting.binomial.multi/proof.md).

## C. Derived BINGO / combinatorial / cross-domain work

1. **Indexed subset products:** The coefficient of x^k y^m in ∏(u_i+v_i y x^{a_i}) equals the weighted number of m-index subset representations of k; the original source distinguishes repeated indices from a set of distinct values. Classical finite representation criterion and powers-of-two uniqueness proven. [Dossier](../sources/math-treasure/research/sum-products/README.md).
2. **q-Euler–Gauss divisibility:** Source claims a complete independently reviewed proof of Narayan–Kiran's Conjecture 4.10 for repeated primes and vanishing polynomial cases, exact statement in its [checked note](../sources/math-treasure/research/connections/breadth-screening/q-euler-gauss/CHECKED.md). **WRITTEN-AUDITED**, no external priority claim.
3. **All-k cyclovariety torus intersection:** Source gives a finite-lattice character-sum/winding proof answering a specific published open question, **WRITTEN-AUDITED**, priority not established; it is cross-domain/topological rather than strict all-integer elementary NT. [Proof](../sources/math-treasure/research/connections/breadth-screening/cyclovariety/R281-candidate.md).
4. **Gauss three triangular numbers:** The derived archive contains an independently reviewed proof via Ankeny's geometric method, importing Dirichlet. This is **not** the owner-requested *new* elementary descent. [Proof](../sources/math-treasure/research/three-triangular-proof/proof.md).
5. **Polygonal number theorem reduction:** The source reconstructs a proof that nonnegative N is a sum of at most m ordinary m-gonal numbers, with at most five terms exceeding one, using the known three-triangular theorem. The needed entirely new elementary descent remains open. [Dossier](../sources/math-treasure/research/figurate-polygonal/README.md).
6. **Shifted-power Toeplitz total positivity:** A locally reviewed theorem outside strict NT gives its complete parameter interval, with endpoint limit unknown; useful methodological inspiration only. [Status](../sources/math-treasure/research/connections/depth-total-positivity-r346/CHECKED.md).

## D. Known deep open problems and unresolved programme targets

**These are not proof results.** They are permissible goals unless they overlap the explicit excluded conjecture.

- **A genuinely new, elementary three-triangular descent:** all n≥0=T_a+T_b+T_c, with a descent mechanism not merely repackaging a standard three-squares proof. The theorem is known, but this specific proof-design target is open in the programme.
- **Seven nonnegative octahedral numbers:** whether every N≥0 is a sum of seven O(n)=(2n³+n)/3, n≥0; source programme has an eight-term reduction and extensive partial range certificates, but no complete seven-term theorem.
- **Rank-two rational Ghost marker isolation:** an all-index theorem when there are two independent rational-point indices and a single denominator-one marker. Examples of difficulty occur in the t=6 and t=12 quartic family.
- **A uniform two-index zero count independent of per-orbit finite hypotheses:** the existing Jacobian/coprime/saturation criteria are usable but must be verified for each orbit; an unconditional criterion would be substantial.
- **General integer matrix two-cube classification:** three cubes are internally reported proved, two cubes over all of M₂(Z) unresolved.
- **Catalan elementary reproof programme:** Catalan's theorem is a historically **proved theorem**, but the sought self-contained elementary substitute for cyclotomic/class-group machinery is **not established** in the copied work. Explicit deficiencies, local-source overlaps and norm barriers remain.
- **ζ(7) irrationality:** still **OPEN** in the archived programme. The archive's ζ(5) irrationality-exponent claims are **conditional** on explicitly cited external inputs, not an unconditional new proof of irrationality.

### Catalan-specific barrier register — never mistake for an elementary proof

The [C35–C38 documents](../sources/math-prime-pro-catalan/INDEX.md) are **PARTIAL**. C35 studies finite harmonic-period ranks and a missing deficiency bound d≤H/4; exact scans do not prove it. C36 distinguishes the actual source-unit overlap J=V∩U₀ from a Fourier-rank deficiency, gives h=dim J and an exact conditional class count p^{n−h}, but **does not** prove h=0 or the uniform bound h≤H/4. C37 exhibits why standard Padé and generator-smallness shortcuts do not force the missing arithmetic conclusion, even with strong exact algebra. C38 uses the second cyclotomic allocation to derive dual necessary inequalities, but their combined contradiction remains missing. The **lesson** is to identify an arithmetic property of the actual hypothetical integer solution, not merely constrain abstract models that also admit false positives.

### ζ-programme barrier register

The [ζ(7) result register](../sources/zeta5-zeta7-research/research/zeta7/RESULTS.md) **supersedes** stronger historical labels. It has seven entries covering ζ(5) conditional irrationality-exponent bounds, alternating integrals, purified matrix families, modular tests and denominator classes. The [mechanism register](../sources/zeta5-zeta7-research/research/zeta7/MECHANISMS.md) records exact failure points for single-species Hankel, cube-integral symmetry, Apéry-type recurrences, V(7,1) lattices and multi-moment constructions. No completed ζ(7) proof exists in this handoff.

## E. Universal verification hazards seen repeatedly

1. An unproved intermediate lemma may **equal the original goal** in disguise (Pell square terms, rank/source-overlap, denominator-one marker).
2. Unit principality does **not** imply an actual global p-th root. C36 contains an explicit counterfixture.
3. A local zero modulo every p^k need not be an integer zero; changing p may help, but no finite collection alone proves global solvability.
4. “No counterexample in checked range” does not prove a global exclusion or effective finiteness.
5. A multi-unit counting theorem often requires **independent markers ≥ orbit rank**; the inequality cannot be dropped.
6. A source's "complete" tag may refer to a finite lemma, a conditional implication, or a proof **under a named hypothesis**.
7. Different versions of an original research paper, or a proof after a dated local review, may alter world novelty/priority. Recheck dated external sources.

**Hand-off principle:** Start each new problem by listing its exact necessary components and the source of each. Do not claim the conclusion until the final missing global step has its own independent proof.