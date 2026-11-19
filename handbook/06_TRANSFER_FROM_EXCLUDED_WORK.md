# 06 — Technique-only transfer from an explicitly excluded historical programme

## Scope rule — do not skip

This handbook **categorically prohibits** work on the Schäffer/Schaffer Conjecture, **including every subcase and parameter-specific variant**. The supporting original documents reside under [reference_only/excluded_problem_source/](../reference_only/excluded_problem_source/). They are **archives of mathematical techniques only**. Do not pursue, verify, amend, cite as an assignment, or repackage the historical problem. In particular, no target-specific index, exponent, endpoint, finite bound, or claimed completed subcase in those files belongs on the new model's to-do list.

Extract only mathematically meaningful, generalized moves below, and try them on wholly unrelated problems after proving the relevant premises independently.

## T1. Exact finite-difference reconstruction of integer polynomial sums

Given an integer-valued polynomial P(n), use Newton's binomial basis
\[
P(n)=\sum_{r\ge0}(\Delta^rP)(0)\binom nr,\quad\Delta P(n)=P(n+1)-P(n).
\]
The telescope \(\sum_{j=0}^{n-1}\binom jr=\binom n{r+1}\) gives an integer-valued primitive for polynomial sequences without destroying denominators early. This can replace Faulhaber/Bernoulli expansions when the goal is an exact divisibility or congruence certificate. The Bernoulli form remains useful for reflection symmetries, degree bounds, signs and factorized centred polynomials.

**New targets:** products of figurate numbers, Stirling/Faulhaber divisibilities, finite binomial sum identities and values of combinatorial counting polynomials. Do not reconstruct the historical excluded equation.

**Archive:** [notebook-methods index](../reference_only/excluded_problem_source/Ramanujans-Notebooks/research/schaffer/notebook-methods/README.md).

## T2. Reflection and centre before expansion

A polynomial sum over consecutive integers often has a hidden symmetry under replacing the endpoint by its negative offset. A centred variable N=2n+1 or T=n(n+1)/2 may turn a high-degree polynomial into a much lower-degree polynomial of N² or T. Work with the *canonical integer centre* and explicitly track denominators/parity.

**Transfer test:** compute P(−n−1), check exact sign, and derive the quotient in Z[T] or a rational ring with a proved integer-valued basis. Do not assume that a rational polynomial in T has integer coefficients.

**New targets:** polynomial representation identities, shifted figurate equations, families of consecutive-power products, integer-valued polynomial congruences.

## T3. Marked-factor graph and forced prime allocations

Factor a hypothetical perfect-power or norm equation into identifiable components. Attach labels describing which primes may divide which components and how valuations of overlaps are forced. Encode the constraints in a finite graph or hypergraph of factors. If a cycle would demand incompatible divisibility or multiplicative gains, eliminate it by an explicit integer invariant.

**What is established by the approach:** algebraic factor graphs, norm-cycle constraints and finite signed-coefficient covers; **what does not follow:** the graph formalism itself does not imply all possible integer solutions are eliminated.

**Transfer test:** prove a correspondence between *every actual solution* and a marked graph, not just examples; prove each graph transition preserves integrality and marker conditions.

**New targets:** simultaneous norm equations, products of consecutive terms equalling high powers, Markov/Hurwitz-type integer dynamics, cyclotomic cofactor systems.

**Archive:** [marked-graphs primer](../reference_only/excluded_problem_source/Ramanujans-Notebooks/research/schaffer/03_marked_graphs/README.md).

## T4. Exact residue-cover certificates rather than numerical searching

Generate a finite set of congruence conditions, each excluding a fully specified class; combine prime-power conditions by CRT only where moduli are coprime. To claim universality, prove that **all** permitted indices or bounded variables lie in the finite union of excluded classes, and give a separately proved analytic or descent reduction if a height bound was required.

A large finite sieve is valuable only for its certified bound. It cannot extend beyond it by an appeal to density or observed absence of counterexamples.

**New targets:** perfect powers in special recurrences, polynomial-value power obstructions, finite Thue reductions, simultaneous congruence and norm conditions.

**Archive:** [exact interval sieve](../reference_only/excluded_problem_source/Ramanujans-Notebooks/research/schaffer/31_endpoint_sieve/README.md).

## T5. Relative units, rational boxes and residue-disk uniqueness

In algebraic norm equations, reduce an object to finitely many unit orbits and use local expansions of unit powers. An integer exponent vector may be constrained modulo p^r by finite polynomial/binomial expansions. With as many independent coefficient conditions as unit directions, a nonzero local Jacobian isolates the known exponent class.

**Critical split:** A genuinely elementary variant derives the ring operations, unit generators, and local expansions by **integer arithmetic and exact finite certificates** (see norm.triples and lifting.binomial.multi). The historical method sometimes employs class-number, ideal or global unit theorems and numerical logarithm bounds; those imports must not be smuggled into a claimed elementary proof.

**New targets:** quartic Thue equations of signature (2,1), integral-point classifications reduced to fixed-norm triples, positive/negative generalized Pell norm constraints.

**Archive:** [relative units](../reference_only/excluded_problem_source/Ramanujans-Notebooks/research/schaffer/16_relative_units/README.md), [integer multi-unit analogue](../sources/math-prime/research/proofs/methods/ghost.unit.multi_unit.md).

## T6. Rational-approximation transfer and clearing denominators

A Ghost statement such as “a special integer expression must be a square” may imply that a rational approximant a/b is implausibly close to a fixed algebraic root. Clear denominators **exactly**:
\[
b^2T-a^2\in\mathbf Z
\]
when T∈Z. Construct an explicit integer polynomial identity forcing that difference to be either 0 or too small relative to a/b. A rational approximation bound must state both the numerator growth and the denominator dependence; missing a b² factor can invalidate it.

**New targets:** square values of triangular/polygonal sequences, integral Thue equations reducible to Pell sectors, near-square integer polynomials and binary cubic/quartic forms.

**Archive:** [Ghost rational-approximation transfer](../reference_only/excluded_problem_source/Ramanujans-Notebooks/research/schaffer/38_ghost_transfer/README.md), [norm-sign rigorous alternative](../sources/math-prime/research/proofs/methods/approx.normsign.md).

## T7. Prime-support allocation between two conjugate norms

When an equation produces two norm/cofactor factorizations, the primes dividing one side may be forced into a congruence class and constrained in the conjugate field. Instead of counting separate local ranks, seek a **joint restriction arising from the same actual integer solution**. Abstract modules that can be chosen independently often fail to provide any contradiction.

**New targets:** repunit perfect-power systems, two simultaneous quadratic-norm equations, catalan-like adjacent powers, multiple power-residue constraints and near-powerful triples.

**Archive:** [Catalan cofactor and two-field warning](../sources/math-prime-pro-catalan/answers/C38.md). This is a **partial method**, not a full Catalan proof.

## T8. Effective height bounds plus exact finite exclusion

The split “large heights ruled out by effective inequalities, small heights ruled out by finite exact computation” is a possible complete method, provided both halves are rigorous and **meet**. A huge bound can be useful only if the smaller-range certificate reaches it; otherwise the gap is mathematically real.

A log/log height bound obtained from algebraic units and linear forms in logarithms is generally **not strictly elementary** without reconstructing the imported theorem. Treat such work as a comparator, or replace its analytic argument with an exact rational-approximation ladder.

**New targets:** constrained Thue equations with exploitable small parameters, parametric perfect-power identities and explicit recurrence-index exclusions.

**Archive:** [height-method context](../reference_only/excluded_problem_source/Ramanujans-Notebooks/research/schaffer/50_odd_square_height/README.md).

## T9. Pell trace algebra and denominator monotonicity

A Pell orbit yields explicit recurrence polynomials, trace coordinates and denominators. One can prove divisibility of denominators under multiplication of indices and sometimes eliminate all indices once the first obstruction on each primitive odd chain is located. To handle a square/Pell marker, combine the recurrence with a *strict* integer gap, not just fixed congruence classes that admit infinitely many returns.

**New targets:** quartic Pell families from figurate equations, simultaneous Pell equations, fixed-power values of Lucas sequences, integer points on conics under two reflections.

**Archive:** [Pell trace research](../reference_only/excluded_problem_source/Ramanujans-Notebooks/research/schaffer/42_pell_trace/README.md), [elliptic component halving](../sources/math-prime/research/proofs/methods/lifting.halving.md).

## T10. Cyclotomic finite transforms, Bernoulli congruences and local rank

A polynomial source's cyclotomic components may be measured using discrete Fourier transforms over finite fields. Frobenius orbits, harmonic values, binomial coefficient recurrences and multiplicative periods provide exact rank lower bounds. The **number of source exponents surviving a modular filter is not by itself the rank of an actual arithmetic intersection**. Keep source classes, primary-unit classes and ideal classes separate.

**New targets:** new exact finite-field rank inequalities, Euler-system-like identities with explicitly elementary proof, cyclotomic divisibility of recurrence values and unit intersection barriers.

**Archive:** [Catalan C35](../sources/math-prime-pro-catalan/answers/C35.md), [C36](../sources/math-prime-pro-catalan/answers/C36.md), [C37](../sources/math-prime-pro-catalan/answers/C37.md).

## T11. Counterfactual models to test whether a proposed argument has enough information

Before attempting a grand statement, construct exact finite algebraic models that satisfy all abstracted axioms yet violate the desired conclusion. If such a countermodel exists, the missing premise must concern an **actual integer/cyclotomic source**, not merely an arbitrary object satisfying the simplified formalism. Similarly, use a p-adic zero or rational point not coming from an integer solution to test the purported local-to-global step.

**New targets:** any new extension of Ghost counting, class principalization, local rank packing or finite-ring sieve.

**Archive:** [source/unit countermodel](../sources/math-prime-pro-catalan/answers/C36.md), [multi-unit false zeros](../sources/math-prime/research/proofs/methods/ghost.unit.multi_unit.md).

## T12. Cross-method thesis

The most promising synthesis is:
1. **Normalization:** finite differences, reflection, recentering, sum-to-product.
2. **Exact decomposition:** Vierzahlensatz, gcd and prime-factor allocation.
3. **Arithmetic orbit:** elementary integer rings, unit transformations or Vieta involutions.
4. **Ghost lift:** replace a nonlinear condition with one or more coordinate markers.
5. **Closure:** multi-index exact lifting, archimedean rational approximation, invariant descent, or a proven finite height certificate.
6. **Independent verification:** exact scripts, full edge cases and hostile review.

The transfer is valuable only when the final target is genuinely **unrelated** to the excluded historical programme. A new proof of an already known result counts only if its *elementary character and novelty of proof technique* survive the full dependency and literature audit.

---

**Operating rule:** Do not read the archive's historical target lists as instructions. The only valid assignments are those in [04_RESEARCH_AGENDA](04_RESEARCH_AGENDA.md) plus new, independently discovered **eligible** problems.