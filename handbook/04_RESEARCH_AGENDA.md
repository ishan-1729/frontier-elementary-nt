# 04 — Research agenda for the frontier mathematics model

**The following is a launchpad, NOT a closed list.** You are specifically instructed to **independently discover additional eligible problems** across elementary number theory, Diophantine equations, recurrences, figurate numbers and elementary combinatorial arithmetic, then rank them against the methods in this repository. **Never pursue the Schäffer/Schaffer Conjecture or any of its subcases, even when a copied reference-only notebook mentions them.**

**Snapshot:** source-repository record through late October–November 2026. Contemporary external problem-status checks were made on 2026-11-20; recheck the latest papers before claiming an open target or an unprecedented elementary proof.

## Priority 1 — A second marker for rank-two rational Ghosts

**Concrete starting point:** the source's quartic equation
\[
(t+1)x^4-ty^2=1,\qquad t=6,12,
\]
with positive integers x,y. In the archive the even-Pell-index halves remain **OPEN**, while the odd-index halves are covered by a norm-sign argument. The Ghost enlargement is
\[
(t+1)B^4-tQ^2=C^4,\quad \gcd(B,C)=1,\quad C=1\ \text{for original points}.
\]
The associated rational Ghost curve carries two independent indices; one marker C=1 does not meet the existing Jacobian dimension condition.

**Your tasks:**
1. Reconstruct the exact two-index Ghost orbit from an elementary integer parent system. Do not assume Mordell–Weil rank as a proof dependency merely because a database reports rank two.
2. Invent a new **genuinely forced** marker: transport across a covering, a fixed numerator/denominator gcd, a conjugate norm coordinate, or index divisibility after doubling. Prove it holds for the original slice and is not algebraically dependent on C=1.
3. On each residue orbit, seek two-marker finite Jacobian isolation, an exact valuation-drift law or a global height gap; use the full finite-polynomial two-index lemma rather than generic Hensel assertions.
4. Solve a meaningful intermediate result first: a new infinite class of indices excluded, a primitive-index reduction, or a uniform statement for a subfamily of t. Then attempt the global classification **with a thorough literature review** (classical results may already classify specific t).

**Why promising:** the underlying methods already solve many rank-one Ghosts and signature-(2,1) two-unit Thue systems. The *rank/marker mismatch* is an explicit, well-defined bottleneck. **High leverage; serious difficulty.**

**Sources:** [quartic family](../sources/math-prime/research/proofs/theorems/quartic_t_family/README.md), [Ghost generalizations](../sources/math-prime/research/proofs/GENERALIZATIONS.md), [two-index lifting proof](../sources/math-prime/research/proofs/methods/lifting.binomial.multi/proof.md).

## Priority 2 — Prove a reusable two-index zero-exclusion theorem

Do not merely hand-check another Thue equation. Improve the *engine* that makes all such checks possible.

Current setup: two markers f(m,n),g(m,n) from multiplication by independent algebraic units, indexed by (m,n)∈Z². The archive proves zero counts when stable reductions are coprime, after a saturation step, or with a frame condition; Jacobian certificates isolate *known* zeros. There exist real instances where a p-adic marker zero never comes from an integer exponent, and the sieve lacks guaranteed termination.

**Targets:** a criterion guaranteeing a finite constellation of primes whose marker conditions isolate every integer root; an auxiliary archimedean inequality closing any residual profinite components; or a structural method creating a third independent marker by passing to a covering.

**Proof test:** state exact field-signature, unit-independence, marker-degree and nondegeneracy hypotheses. Prove the conclusion for *all* integer indices. Show counterexamples if the new theorem fails without one hypothesis.

**Why promising:** one general theorem would apply across quartic Thue, simultaneous Pell, norm triples and potentially rank-two elliptic Ghosts. **Highest methodological payoff.**

**Source:** [full two-index result](../sources/math-prime/research/proofs/methods/lifting.binomial.multi/proof.md).

## Priority 3 — Elementary finite-norm certificates for new small Mordell/Thue classes

**Objective:** Create an algorithm that, for selected cubic quotient rings given by monic f(θ)=0, automatically generates a *human-checkable exact* Euclidean or finite-multiplier norm cover, complete Bézout/prime allocation, and a finite unit generator proof.

1. Search among small discriminant cubic equations and Mordell \(Y^2=X^3+k\) cases **not already classified elementarily** in the archive.
2. Construct explicit 3×3 multiplication matrices and norm polynomials.
3. Find a rational box cover with strict norm margin; prove every row symbolically, including half-lattice exceptions.
4. Derive all unit representatives by trace-square descent, then prove an all-exponent coefficient exclusion by exact binomial valuations.
5. Distill a **uniform transfer theorem** and finish at least one previously unproved-or-non-elementarily-proved case.

**Guard:** A class-number-one assertion is not permitted until it is itself reconstructed from exact finite division; do not assume UFD or use a package's algebraic integer solver as the proof.

**Why promising:** three sharply distinct verified proof prototypes are archived. A new arithmetic example could be a tractable independent theorem, and a generator of norm certificates is a major method extension.

**Sources:** [norm triples](../sources/math-prime/research/proofs/methods/norm.triples.md), [non-Euclidean multiplier](../sources/math-prime/research/proofs/methods/norm.triples.multiplier.md), [small Mordell completion](../sources/math-treasure/research/mordell-thue-completion/README.md).

## Priority 4 — Erdős Problem 364: three consecutive powerful numbers

A positive integer N is *powerful* if every prime p dividing N also satisfies p²|N. The question asks whether **any** n−1,n,n+1 can all be powerful. The problem page remains open in its 2026-11 snapshot: [Erdős 364](https://www.erdosproblems.com/364). Proved special configurations and extensive finite verification do not settle it.

**New route:** write every powerful number uniquely as a²b³ with b squarefree and a,b positive. Consecutive differences give
\[
a_+^2b_+^3-a_0^2b_0^3=1,\qquad
a_0^2b_0^3-a_-^2b_-^3=1.
\]
Treat **both** equations jointly. Introduce Ghost squarefree kernels b_−,b_0,b_+ and a forced prime-support allocation; try norm-pair/triple relations, a signed graph contradiction, or a second nontrivial marker on a simultaneous Pell family.

**First outputs:** prove a novel infinite subclass excludes triples (e.g. a condition on gcd, small squarefree kernel, or a prime divisor in the middle), or derive a new unconditional inequality whose hypotheses are strictly weaker than published 2025 cases. **Do not propose brute-force enumeration as a route to the universal negative statement.**

**Adjacent target:** [Erdős 365](https://www.erdosproblems.com/365), whether every pair of consecutive powerful integers has a square among them, is also structurally suited to two-term Pell/Ghost methods; verify current status.

**Why promising:** adjacent factors, square/cube allocations and multi-constraint Ghosts are the archive's specialty, though the full question is difficult.

## Priority 5 — Goormaghtigh: one repunit in two distinct bases

Let R_n(x)=1+x+\cdots+x^{n-1}=(x^n−1)/(x−1). Seek distinct x,y≥2 and lengths m,n≥3 with R_m(x)=R_n(y). The two famous examples are
\[
R_5(2)=R_3(5)=31,\qquad R_{13}(2)=R_3(90)=8191.
\]
Whether these are the only ones remains open in the available 2026 literature; [original-style status and bounds](https://arxiv.org/abs/2410.03677).

**Attack:**
- Multiply through: (y−1)(x^m−1)=(x−1)(y^n−1), which is an **AB=CD** signature.
- Allocate gcd(x−1,y−1), cyclotomic factors and prime powers before invoking Vierzahlensatz. Test differences of near powers and growth separation.
- For one exponent fixed or in a restricted congruence class, seek a uniform elementary Vieta or norm descent.
- Derive a new classification of a genuinely broader parameter family, not an already known m=3 or n=3 case.
- Check literature very carefully: partial solutions for fixed exponents are extensive and some involve hard logarithmic-form estimates.

**Why promising:** perhaps the clearest fresh use for the user's sum→product→four-factor→descent framework. Full theorem is a long-shot; **new fixed-parameter elementary families** may be realistic.

## Priority 6 — Nagell–Ljunggren repunit perfect powers

The equation
\[
R_n(x)=y^q,\quad x>1,\ n>2,\ y>1,\ q\ge2
\]
has the classical examples \(R_5(3)=11^2\), \(R_4(7)=20^2\), \(R_3(18)=7^3\). The unrestricted classification remains unresolved in current source summaries; [OEIS status](https://oeis.org/A208242).

**Method match:** cyclotomic/gcd factor splitting, finite-index Ghosts, prime-support allocation, exact unit exponent lifting, and a Wronskian Padé gap when cyclotomic factors yield Thue equations.

**First theorem targets:** a new exponent/residue subclass, an elementary replacement for a known argument using class numbers, or a lifting criterion forcing a primitive prime to land in a forbidden factor.

**Guard:** Do not incorrectly identify the full theorem from the three examples or assume a primitive-divisor theorem without reconstructing it elementarily.

## Priority 7 — Markov uniqueness via Vieta orbit labels

Markov triples are positive integers satisfying
\[
a^2+b^2+c^2=3abc.
\]
The Vieta involution c↦3ab−c maps a solution to a solution. The Markov uniqueness conjecture asks whether the **largest** Markov number determines its triple up to permutation. The general conjecture remains open according to a November 2026 trillion-value computational report; [recent report](https://hexagonmath.org/pdf/2611.00410v1).

**New route:** Treat the tree of Vieta involutions as noncommuting Ghost moves, find an integral orbit invariant that distinguishes paths, or a prime-support residue tag whose collisions lead to a strict descent. Prove uniqueness for a new **structural infinite family**, not just another larger computation.

**Why promising:** the archive's weighted Vieta, signed-level, preserved gcd and full-orbit graph arguments are directly relevant. **High-concept, hard.**

## Priority 8 — Seven octahedral numbers, or an improved unconditional eight-term proof

Let O(n)=(2n³+n)/3 for integers n≥0. The archived result gives **every** N≥0 as a sum of eight O(n) **conditional on the named special three-square theorem GL3(1 mod4)**; its seven-summand target remains unresolved.

**Tasks:** Improve the conditional eight-term construction to a strictly elementary three-square input (or prove that input independently via the allowed tools); invent a seven-term covering identity with an exact finite-range return bound; or prove a nontrivial improved seven-term interval with a **proved** finite bound and lifting rule.

**Key starting identity:** paired terms O(c+d)+O(c−d) are a polynomial in c and d². Expand exactly and look for a three-squares/three-variables certificate, while preserving 0≤d≤c. The current eight-term large-range proof uses three pairs plus two singles; seven needs a genuine new compression.

**Status:** eight **WRITTEN-AUDITED** under GL3; seven **OPEN** in archive. [Source theorem and boundary](../sources/math-prime/research/proofs/theorems/octahedral_bounded_sums/README.md).

## Priority 9 — An independently elementary descent for Gauss's three triangular numbers

Known theorem: each n≥0 is a sum of three T_j=j(j+1)/2 for j≥0. The existing derived project uses an **Ankeny/lattice route with a Dirichlet import**; that is not a fresh elementary descent.

**Research target:** Starting from
\[
8n+3=(2a+1)^2+(2b+1)^2+(2c+1)^2,
\]
construct a **well-founded descent on integer triples** that both preserves oddness and eventually reaches a certified base case without importing three-squares or primes-in-progressions in the final proof. The obvious restatement of the three-squares theorem is **not** a proof.

**Why promising:** a known truth permits all methods to be audited against complete ground truth; a genuinely new descent could transfer to m-gonal representations and other ternary forms.

**Sources:** [existing non-descent route](../sources/math-treasure/research/three-triangular-proof/proof.md), [figurate reductions](../sources/math-treasure/research/figurate-polygonal/README.md).

## Priority 10 — Integer 2×2 matrices as sums of two cubes

The archive reports three cubes sufficient in M₂(Z) and supplies Cayley–Hamilton reduction:
\[
X^3=(\mathrm{tr}(X)^2-\det X)X-\mathrm{tr}(X)\det(X)I.
\]
Find exact arithmetic criteria for A=X³+Y³ with integral X,Y. Identify local obstructions, whether local-global sufficiency fails, and a complete classification in restricted trace/determinant strata.

**Why promising:** concrete polynomial Diophantine system, exact finite arithmetic certificates, and an already completed three-cube framework. Check outside literature before claiming any two-cube result novel.

**Source:** [three-cube card](../sources/math-prime/research/proofs/theorems/matrix_cubes/README.md).

## Priority 11 — A controlled reattack on the elementary Catalan gap (speculative)

Catalan's theorem has already been proved using advanced mathematics. The research challenge here is **an independent elementary proof**, not settling the known theorem anew.

The archival C36 source-count interface shows a missing bound on the actual source-primary-unit intersection J. Merely bounding Fourier deficiency, proving principality or adding independent ideal-field counts does not close it. C37 shows ordinary Padé auxiliaries fail to achieve enough cancellation in a wide exact class. C38 gives a second radical with joint prime-cofactor restrictions but no full source contradiction.

**Only productive next steps:** (a) an explicit **actual-source** arithmetic restriction forcing h≤H/4 or h=0, with no assumed unit splitting, or (b) a qualitatively different elementary descent that avoids cyclotomic prerequisites. Prove a genuinely new family or obstacle lemma first, and **do not claim** an all-band result from finite samples.

**Priority:** High conceptual impact, **low near-term solvability**. [C36](../sources/math-prime-pro-catalan/answers/C36.md), [C37](../sources/math-prime-pro-catalan/answers/C37.md), [C38](../sources/math-prime-pro-catalan/answers/C38.md).

## Priority 12 — New elementary cyclotomic / q-divisibility generalizations

One archived proof claims q-Euler–Gauss divisibility including repeated primes and zero factors; another proves a cyclovariety unit-torus intersection using finite exponent symmetry. Abstract the **finite divisibility and coefficient cancellation** parts (not the topological or categorical imports), then discover new arithmetic divisibility statements for polynomial sequences, q-binomial recurrences, and finite Möbius transforms.

**First theorems:** a new divisibility closure under composition or repeated prime powers with exact integer polynomial certificates, or a universal valuation congruence not already in the literature.

**Status:** method transfer, not a claim those related open questions are already solved. [q-Euler–Gauss checked note](../sources/math-treasure/research/connections/breadth-screening/q-euler-gauss/CHECKED.md), [cyclovariety](../sources/math-treasure/research/connections/cyclovariety-torus/README.md).

## Priority 13 — Long-shot exploratory tracks (do not displace tractable work)

- A new **purely elementary proof of Fermat's Last Theorem**: historically known theorem with advanced standard proof. Focus on a precise reusable hypothetical descent or a new obstruction class; avoid a circular primitive-coprimality lemma.
- ζ(5), ζ(7) and related irrationality/linear form questions: the archive lists strong partial and conditional work but no full ζ(7) result. Try only if a *new arithmetic positivity-and-denominator mechanism* can beat the proved barriers; [mechanism register](../sources/zeta5-zeta7-research/research/zeta7/MECHANISMS.md).
- Two independent perfect-power constraints on Pell/Lucas sequences, factorial or figurate recurrences, and simultaneous Thue equations not currently in this agenda; **screen each one independently**.

## The independent problem-discovery engine — mandatory

You are **not limited** to Priorities 1–13. Continuously run a second research loop:

1. **Harvest signatures** from techniques: AB=CD with adjacent factors; integer products nearly r-th powers; quadratic-in-one-variable involutions; cubic norm with small discriminant; a rank-r orbit with r vanishing coefficient markers; n²+d=y³; two independent recurrence indices; a perfect-power condition on a symmetric progression.
2. **Search eligible sources for unanswered questions:** the Erdős Problems database, published Diophantine/recurrence papers, OEIS comments as *leads only*, formal-conjectures, specialized journals and their latest arXiv versions. Prefer problem statements carrying an explicit missing lemma and modest enough structural complexity to match an exact elementary technique.
3. **Match each question to an explicit proof pipeline** from the playbook. Name the missing lemma and why it is stronger/weaker than known results.
4. **Discover at least three new candidates of your own** before committing to a single moonshot, not counting the agenda items above.
5. **Rank candidates** from 0–5 by (i) elementary fit, (ii) existing method match, (iii) likelihood of a first substantive lemma, (iv) potential mathematical significance, (v) novelty after literature search, (vi) checkability/formalization. A good project has both a plausible first lemma and a major expansion if it works.
6. **Choose a portfolio:** one tractable complete-family problem, one method-level generalization, and one longer-shot open problem. Replace stalled tracks when a countermodel shows the toolkit cannot yet close them.
7. After each successful lemma, **reverse-map it** to new problems whose structure now satisfies its hypotheses. Recompute the agenda. Prefer a theorem that closes several families over a bespoke isolated computation.
8. **Respect the hard exclusion on every new discovery.** Do not reinterpret an excluded special case as a “related new problem.”

## Proof milestones and output format

**Milestone A: exact formulation and current prior art.** A one-page statement with all integers/ranges/signs, and dated references. Confirm that the target is truly open or lacks a genuinely elementary route.

**Milestone B: proof graph.** Every arrow in normalization→orbit→height→seed→marker→exclusion→reconstruction has a named lemma. Tag any imported non-elementary theorem and any unresolved claim.

**Milestone C: exact certificates.** Supply polynomial identities symbolically, integer norm tables with strict margins, recurrence-period and valuation tests with finite algebra proofs, and an executable verifier that produces deterministic results.

**Milestone D: hostile independent review.** Try to break assumptions at p=2, zeros, negative indices, primes dividing discriminants, signed reflection boundaries, repeated roots, torsion points, p-adic false positives and rank-two missing markers.

**Milestone E: finished result.** A standalone elementary argument, reproduction script, proof status, novelty limitations and at least one newly identified family for transfer.

**No quiet success inflation:** "proved" means every infinite-range step is actually proved, not inferred from a large finite computation. If a route fails, report the sharpest true lemma, the counterexample and the next mathematically different idea.