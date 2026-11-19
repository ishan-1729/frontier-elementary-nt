# 01 — Elementary Number Theory method playbook

**How to use this file:** Every item specifies a *trigger*, a concrete transformation or lemma, proof obligations, transferable examples, and an internal source. Study the original copied cards and proofs rather than accepting this overview as a substitute. Unless tagged as a classical identity, a programme-specific result retains the audit status of its cited source.

## 1. The sum → product → coprime allocation cascade

**Trigger:** A Diophantine sum, polynomial power sum, difference of structured terms, or symmetric progression conceals factors that are nearly coprime.

**Moves:**
- Centre a progression: a−d, a, a+d or a−3d, a−d, a+d, a+3d. For even-length progressions use twice the centre to keep integer/parity bookkeeping exact.
- Use x²−y²=(x−y)(x+y), x³−y³=(x−y)(x²+xy+y²), x³+y³=(x+y)(x²−xy+y²); use aⁿ−bⁿ and cyclotomic factors only with fully justified gcd information.
- After F=UV with gcd(U,V)=1 and UV=s w^r where s is r-th-power-free, valuations force **each** prime power into a unique U or V and leave finitely many distributions of primes dividing s. For s=3 and r=2, a coprime product p q=3d² gives either p=u²,q=3v² or p=3u²,q=v², up to signs and signs of d, with u,v coprime.

**Proof obligations:** exact identities, primitive normalization, gcd formula, parity/zeros, allocation of *every* prime dividing the coefficient, and an actual descent or last contradiction after parameterizing. A product identity alone is not proof.

**Sources:** [MathAgent sum-to-product](../sources/MathAgent/knowledge/methods/sum_to_product.md), [factor splitting](../sources/MathAgent/knowledge/methods/euclid_lemma_factor_splitting.md), [arithmetic progressions](../sources/MathAgent/knowledge/methods/arithmetic_progressions.md), [derived BINGO products](../sources/math-treasure/research/sum-products/README.md).

## 2. Vierzahlensatz (“Four Numbers Theorem”) — project-specific terminology

**Exact statement:** Given integers A,B,C,D satisfying AB=CD, one can choose integer parameters p,q,r,s with

\[
A=pq,\quad B=rs,\quad C=pr,\quad D=qs.
\]

One standard construction for nonzero positive A,C: let p=gcd(A,C), q=A/p, r=C/p. Since gcd(q,r)=1 and AB=CD implies qB=rD, deduce r|B and q|D. Put s=B/r=D/q. Handle general signs by signed gcd/units and handle zeros separately.

**Why it helps:** Preserve equality while making two gcd systems visible. For example, if C,D are consecutive, gcd(C,D)=1 forces gcd(p,q)=1 as p|C,q|D; then coprime square allocation for A can determine p,q.

**Example entry point:**
\[
3d^2(4d^2+1)=c(c+1)
\]
becomes
\[
3d^2=pq,\quad 4d^2+1=rs,\quad c=pr,\quad c+1=qs.
\]
Then gcd(p,q)=1 and the above two square allocations are available; eliminate c from qs−pr=1. The upstream example is a **method outline**, not the missing full proof of its motivating Diophantine statement.

**Generalizations:** rectangular decompositions for AB=CD involving multiple product identities; repeated gcd graphs; coefficient allocation in higher perfect powers; cofactor squares on 2-coverings; expose an involution after parametrization.

**Sources:** [Vierzahlensatz method](../sources/MathAgent/knowledge/methods/vierzahlensatz.md), [worked outline](../sources/MathAgent/knowledge/examples/vierzahlensatz/kieren_mse_x2_plus_1_eq_y3.md), [primer and caveats](../sources/Ljunggren/Vierzahlensatz_Primer/README.md).

## 3. Bézout substitution and integer-coordinate changes

If gcd(m,n)=1, choose integers u,v with mu+nv=1. Use them to reparameterize linear divisibility conditions mx+ny=k, taking all solutions (x,y)=(ku+nt,kv−mt), t∈Z. Integrality and sign restrictions must be carried through both directions. A unimodular integer substitution is particularly useful after Vierzahlensatz, where one equation says a combination of two factored quantities is ±1. Rational substitutions with uncleared denominators are **not** equivalent.

**Source:** [Bézout method card](../sources/MathAgent/knowledge/methods/cauchy_bezout_substitution.md).

## 4. Discrete squeeze and adjacency barriers

**Idea:** A required integer perfect power is placed strictly between consecutive powers. Examples: n² < M < (n+1)² or n³ < M < (n+1)³, with n≥0. More generally, a recurrence candidate lies strictly between adjacent values of a monotone integer sequence.

**Procedure:** (i) derive exact integer expression M; (ii) choose the right adjacent values; (iii) establish positivity and strict bounds on the entire claimed range; (iv) check endpoints separately. Recurrences from Pell reduction often supply the adjacent bounds. If bounds are only non-strict, equality is a possible solution and must be treated.

**Source:** [squeeze-sequence card](../sources/MathAgent/knowledge/methods/squeeze_sequence.md), [Pell reduction](../sources/MathAgent/knowledge/methods/pell_equation.md).

## 5. Classical descent and Vieta/weighted involution descent

**General descent contract:** a solution S maps to S' with the same equation, integer positivity/domain preserved, and a natural number height h(S')<h(S). A mere involution or numerical reduction seen in examples is not enough.

**Vieta move:** For quadratic equation Az²+Bz+C=0 in z, companion root z'=−B/A−z. First prove it integral and nonnegative, then verify it reduces height. This extends to conics with two coordinate-wise reflections.

**Reusable weighted theorem template:** Let a,b,c,C positive integers, with m=c/a and ell=c/b integers. Put
\[
F(x,y)=ax^2+by^2-cxy.
\]
The integer involutions
\[
X(x,y)=(my-x,y),\qquad Y(x,y)=(x,\ell x-y)
\]
preserve F. For F(x,y)=C,
\[
a x(my-x)=by^2-C,\qquad b y(\ell x-y)=ax^2-C.
\]
Both jumps cannot fail to reduce the corresponding positive coordinate, for this would give −2C≥0. If in addition 0<C<a+b+c, a proposed jump crossing to opposite-sign nonzero coordinates is impossible because the transformed positive form a u²+cuv+b v² is ≥a+b+c. Thus positive points descend through nonnegative points to axis seeds using height x+y. The original dossier goes further: congruence, gcd and square-class restrictions, signed levels and parametric orbit classification.

**Key warnings:** For C outside the sign barrier, branch crossings change the argument. The all-quotient positive-offset classifier is an additional theorem with separate hypotheses; one cannot apply the axis theorem outside its range.

**Sources:** [original Vieta source dossier](../sources/math-treasure/research/vieta-descent/README.md), [signed levels](../sources/math-treasure/research/vieta-descent/signed-levels.md), [all positive offsets](../sources/math-treasure/research/positive-offset-classification/proof.md), [conic fundamental arc](../sources/math-prime/research/proofs/methods/involution.vieta.md).

## 6. The Ghost parametrization schema — the programme's central synthesis

A Ghost is an enlarged arithmetic object with an explicit relation to the original solutions; it is not an invented solution. **Each proof must supply eight objects:**

1. **G** — enlarged Ghost domain, including extra rational points, primitive Pythagorean nodes or integral tuples of a norm equation.
2. **P** — a Ghost Property invariant under the moves, such as a norm, primitivity, integrality or a sign/coprimality condition.
3. **M** — reversible or directed algebraic moves (Pythagorean parents, chord/tangent, explicit 2- or 3-isogeny formulas, unit multiplication, integer involutions).
4. **h** — positive integer descent height strictly lowered by a valid parent outside a finite set.
5. **S** — explicit finite seed set, **complete for all Ghosts** in the claimed domain.
6. **Gamma** — orbit/index group: Z for one signed index, Z^r for r independent units/moves, potentially a tree of involutions.
7. **mu** — an *if and only if* marker that recovers the original equation, e.g. denominator B=1 or coordinate w=0.
8. **E and rho** — all-index exclusion of spurious marker values and exact reconstruction of the original variables.

**Why enlarging works:** Exact integrality or perfect-power properties may be hard to test directly in the original variables but become the fixed value of one integer observable on a covering. First solve the bigger orbit geometry; then prove that the marker is rare. If the marker remains “is a square modulo p”, no fixed p-adic precision generally excludes it; construct a further covering to turn the square into a coordinate instead.

**Worked classes:** Ljunggren Ghost nodes r⁴+q²=2s⁴, with original slice r=1; reduced rational points on Y²=X³−7 with original denominator B=1; negative-sign Mordell fixed-norm integer triples with marker w=0; 2-descent quartic Thue integer triples with two vanishing markers.

**Sources:** [unified design](../sources/math-prime/research/proofs/GENERALIZATIONS.md), [rational ghosts](../sources/math-prime/research/proofs/methods/ghost.rational.md), [unit ghosts](../sources/math-prime/research/proofs/methods/ghost.unit.md), [audited Ljunggren Ghost route](../sources/math-treasure/research/ghost-classification-audit/proof-audit.md).

## 7. Exact finite lifting on recurrence orbits

**Trigger:** An integral polynomial transition T on a finite tuple (including a scaled Somos window) satisfies T^m(z)≡z modulo p^a on the relevant residue classes, and a marker f(T^{mk+c}z) has a primitive nonzero first increment.

**Core idea:** Let J=T^m and write locally J(z)=z+p^a Lambda(z), where Lambda has integer polynomial coefficients and a≥1 (for the p=2 doubling step usually require a≥2). Binomial finite differences of iterates acquire additional p-adic valuations. For a fixed seed residue class c, prove an exact drift law of shape
\[
v_p(f_{c+Mq}-f_c)=a_0+v_p(q)
\]
for **all nonzero integers q**. If f_c is the target value and differences never vanish, the marker is excluded on the entire infinite class except at q=0.

**Documented examples:**
- For x²+7=y³, v₅(s_{c+5q}−s_c)=2+v₅(q) on specified seed classes c=0,−1.
- In the quartic parameter t=7 case, v₂(a_n−a_k)=v₂(n−k)+3 when n≡k mod8.
- For x²+x+1=y³, the prime 3 is **misleading** (persistent high-precision residues); the prime 2 with period 8 yields v₂(s_{c+8q}−s_c)=4+v₂(q) on its seed classes.
- The Ljunggren proof uses an index-dependent modulus 3^{e+3}, never a fixed finite scan in lieu of an all-index argument.

**Proof obligations:** integral polynomial transition (possibly after scaling), primitive denominator/numerator lemma for the marker, period/residue partition, sign of q, exact higher-difference valuations, and proof covering every residue class. At p=2, a map ≡identity mod2 alone may fail the doubling estimate.

**Source:** [finite Somos lifting](../sources/math-prime/research/proofs/methods/lifting.somos.md).

## 8. A portable one-unit 3-adic coefficient theorem

Let R be a commutative ring free of finite rank over Z with distinguished basis 1,e₂,...,e_d. Let L extract one nonconstant basis coefficient. If eta=1+3w is a unit in R, L(w) is not divisible by 3, and eta^n exists integrally for every n∈Z, then for every n≠0:
\[
v_3(L(\eta^n))=1+v_3(n).
\]

**Proof for n>0:** the linear term in (1+3w)^n is 3nL(w). For j≥2,
\[
v_3\bigl(\binom nj 3^j\bigr)\ge v_3(n)+j-v_3(j)\ge v_3(n)+2,
\]
so higher terms cannot cancel the leading term. For n<0, eta^{-1}≡1 mod3, so eta^{-1}=1+3w' with w'≡−w mod3; apply the positive case to |n|. This excludes every nonzero exponent that would force L(eta^n)=0. It is a **general arithmetical lifting theorem** independent of any particular Mordell equation.

**Source, with complete norm and coefficient examples:** [small Mordell–Thue completion](../sources/math-treasure/research/mordell-thue-completion/README.md).

## 9. Multi-index finite polynomial lifting and the dimension test

For r independent unit exponents and t independent markers, the search becomes zeros of t functions on Z^r. The available zero-isolation strategy works most naturally when **t≥r**.

- Rational elliptic Ghosts generally have a single denominator-one marker: at Mordell–Weil rank 1 the one-index technique is feasible; at rank 2 a single marker cannot be isolated by a Jacobian full-rank condition.
- A degree-n Thue equation gives n−2 nonconstant coordinate-zero markers. A unit orbit has rank r_U=r₁+r₂−1. Therefore n−2−r_U = r₂−1. This dimension balance is favorable when the defining Thue field has at least one complex place; for quartic signature (2,1), **two markers and two unit indices** align.
- For two markers, two indices, the archive supplies exact finite-polynomial zero-count theorems if the stable reductions are coprime, a saturation certificate holds or a frame hypothesis holds. A **Jacobian determinant nonzero mod p** at a known solution isolates that p-adic zero on its exponent class.
- Use unit reductions modulo prime powers to construct kernel lattices Lambda_q=ker(Z²→(R/qR)^×). Enumerate exponent cosets; isolate known roots via nonzero Jacobians; exclude remaining cosets through other moduli or a proved global bound.
- **A sieve may never terminate**: there can be a nonintegral p-adic zero for every p-power. A second independent prime, an archimedean height barrier or a new marker is sometimes essential.
- Neither determinant nonvanishing nor finite termination follows automatically from the abstract setup. The source has explicit counterexamples.

**Source:** [two-index theorem and exact hypotheses](../sources/math-prime/research/proofs/methods/lifting.binomial.multi/proof.md), [method card](../sources/math-prime/research/proofs/methods/lifting.binomial.multi.md), [multi-unit orbit coverage](../sources/math-prime/research/proofs/methods/ghost.unit.multi_unit.md).

## 10. Integer-pair and integer-triple arithmetic in place of black-box fields

**Integer pairs:** encode a+bω by the pair (a,b) and a fixed quadratic multiplication table. Prove conjugation, norm multiplicativity, Euclidean division via exact rounding, gcd and Bézout, units, and coprime cube/square splitting **from scratch**. This can replace statements about algebraic integer rings for small Mordell equations. **Counterexample:** a pair-ring with discriminant −23 does *not* admit naive cube extraction even when a conjugate pair is coprime in a simplistic sense; do not infer UFD from a norm factorization alone.

**Integer triples:** encode a+bθ+cθ² by (a,b,c), with θ satisfying an explicit cubic. Reduce products modulo that cubic. Let M_α be the 3×3 integer matrix of multiplication by α; define N(α)=det(M_α), and α* from the adjugate so α α*=N(α). This constructs norms and inverse denominators by finite integer arithmetic, no field invariant tables.

**Euclidean certificate:** tile the fundamental rational box of residue triples with finitely many rational subboxes, and for each provide an integer translate with a verified upper bound |N|<1 (via exact polynomial Taylor bounds). For θ³=4θ−2, the archived certificate uses 24 boxes with maximum verified bound 2023/2048. From exact division derive Bézout, Euclid's lemma, coprime square allocation and units.

**Beyond norm Euclidean:** a Dedekind–Hasse-type finite-multiplier certificate finds s in a finite set S and kappa in the integer triple ring with 0<|N(s xi−kappa)|<1 for every relevant residue xi, including half-lattice exceptional cosets. The archive uses S={1,2} for cases not norm Euclidean, with a margin 8191/8192 in one case. The multiplier condition has to be proved on **all** residues; arbitrary integer multipliers do not always suffice.

**Units by elementary descent:** seek explicit units shrinking an elementary positive height (e.g. trace-square) outside a bounded box. Enumerate every unit in the box exactly, yielding finite generators without importing a unit theorem. For multi-unit fields use an explicit reduced cone in logarithmic/archimedean coordinates with rational interval certification.

**Sources:** [norm pairs](../sources/math-prime/research/proofs/methods/norm.pairs.md), [norm triples](../sources/math-prime/research/proofs/methods/norm.triples.md), [finite multiplier method](../sources/math-prime/research/proofs/methods/norm.triples.multiplier.md), [unit Ghosts](../sources/math-prime/research/proofs/methods/ghost.unit.multi_unit.md).

## 11. Signed orbit coverage, chord/tangent moves and real-component halving

Rather than assuming an elliptic-curve group law, derive chord/tangent coordinate maps explicitly and verify their polynomial identities. Show every rational Ghost has an integral-height-decreasing parent; prove the signed index covers **all** nodes by induction in both directions.

When the real cubic has three roots, it has bounded and unbounded real components. If odd-index points lie on the bounded component ('egg') and a doubling map preserves denominator divisibility B(P)|B(2P), then integral even indices halve back to integral points. Integer abscissas on the bounded egg occupy a finite interval, so enumerate all odd base cases; one nonintegral double blocks every higher double in that chain. This is a genuinely **archimedean** all-index exclusion, not a p-adic sieve.

**Source:** [Ghost rational](../sources/math-prime/research/proofs/methods/ghost.rational.md), [component halving](../sources/math-prime/research/proofs/methods/lifting.halving.md), [Mordell 1963 complete proof](../sources/math-prime/research/proofs/theorems/mordell_1963/proof.md).

## 12. Rational approximation with a polynomial Wronskian sign

**Trigger:** A Thue-type form can be reduced to F(U,V)=Re((U+iV)^m)−a·Im((U+iV)^m), whose small root is rho_m(a)=tan(arctan(1/a)/m).

Construct rational ladders T_n(a)/S_n(a) through
\[
S_0=1,\ S_1=2ma,\quad T_0=0,\ T_1=2,
\]
and, for Y∈{S,T},
\[
Y_{n+1}=\frac{2m(2n+1)}{n+1}aY_n+
\frac{4(m^2n^2-1)}{n(n+1)}Y_{n-1}.
\]
The source's polynomial Wronskian identity controls signs, alternating brackets, denominator growth and an explicit rational approximation lower bound. For degree m and n-dependent denominators, multiply by d_m(n)=∏_{p|m}p^{v_p(n!)} before claiming integer separation. For m=3 unnormalized T₃ has coefficient 140/3: naive integrality is false. Once the rational approximation exponent beats the Diophantine target’s degree, finitely many small U remain; finish by exact congruences.

**Sources:** [norm-sign ladder](../sources/math-prime/research/proofs/methods/approx.normsign.md), [all-degree extension](../sources/math-prime/research/proofs/methods/approx.normsign.all_degree.md), [Ljunggren full independent proof](../sources/math-prime/research/proofs/theorems/ljunggren/norm_sign/elementary_proof.md).

## 13. Binomial/Padé gap and the finite separation problem

For fixed x²+d=y³, the archive's positive-d case reduces by pigeonhole and Bézout to finitely many cubic equations, then compares large integer solutions with a fixed one using explicit binomial approximants. If two required nonzero integer pair vectors would each have squared length strictly below 1, contradiction. The method proves finiteness under the stated source hypotheses; it does **not automatically list all solutions** and may be ineffective if the constants are not explicit.

The Catalan programme also builds an instructive **negative theorem**: ordinary Padé denominators for a particular cyclotomic source polynomial do not achieve enough norm-degree cancellation, even at large degree. For a new problem, compute the exact required cancellation order and compare it to a coefficient-count lower bound **before** assuming a Padé device closes the argument.

**Sources:** [binomial approximation](../sources/math-prime/research/proofs/methods/approx.binomial.md), [Catalan C37 obstruction and calculations](../sources/math-prime-pro-catalan/answers/C37.md).

## 14. Pell recurrence, Vieta orbit and index descent

A quadratic square-discriminant condition leads to X²−D Y²=N after clearing denominators and completing a square. Develop an elementary Pell recurrence and its integer invariants. Check whether:
- a power/property at a composite index forces a smaller index (index descent),
- a residue class has exact valuation drift,
- a recurrence lies between consecutive squares,
- a Vieta companion root forces a smaller positive solution,
- or a larger cover provides a **linear** marker.

**Critical limitation:** classifying square terms of a Pell sequence by finite search is generally inadequate; a supposed “only these squares” lemma may be exactly the hard theorem being proved. The archive explicitly records a prior failed soft-proof shortcut of this kind.

**Sources:** [Pell method](../sources/MathAgent/knowledge/methods/pell_equation.md), [squeeze](../sources/MathAgent/knowledge/methods/squeeze_sequence.md), [Vierzahlensatz primer warning](../sources/Ljunggren/Vierzahlensatz_Primer/README.md).

## 15. Finite differences, Faulhaber/Bernoulli, binomial bases and denominators

For any integer polynomial P(n), expand in the binomial basis:
\[
P(n)=\sum_{j\ge0}\Delta^jP(0)\binom nj.
\]
Since \(\sum_{k=0}^{n-1}\binom{k}{j}=\binom{n}{j+1}\), exact integer-valued power-sum expressions follow without unnecessarily introducing Bernoulli denominators. The method combines with centered terms, reflection n↦−n−1, factorial valuations, high-precision congruence identities and marked factor gcds.

**Use as transfer, not a target:** the richest development arose during an off-limits conjecture programme. See [06_TRANSFER_FROM_EXCLUDED_WORK](06_TRANSFER_FROM_EXCLUDED_WORK.md). No question from that programme or its subcases is permitted.

## 16. Marked-factor graphs, modular supports and whole-residue covers

If a family has forced factors or norm cofactors, form a graph of prime allocations (factor nodes, residue-class labels, local valuations). Elementary reciprocal/period restrictions can prohibit cycles or force a “fresh prime” in each hypothetical solution. Test not only local coverage but whether a global integer lift exists. CRT provides compatible congruences only for a **finite selected set of moduli**; simultaneous all-modulus compatibility still may not imply an integral solution.

A useful certificate for an **infinite** exclusion is either (a) a proved bound reducing to finitely many classes, (b) an index-drift valuation, (c) a finite-state recurrence whose classes all terminate, or (d) an analytic gap proving a finite height ceiling. Avoid extrapolating density statements to exceptions unless the statement explicitly implies none.

**Source:** [reference-only transfer](06_TRANSFER_FROM_EXCLUDED_WORK.md), [finite lifting](../sources/math-prime/research/proofs/methods/lifting.somos.md).

## 17. Integer matrices and Cayley–Hamilton linearization

For X∈M₂(Z), with t=tr X and d=det X,
\[
X^2=tX-dI,\qquad X^3=(t^2-d)X-tdI.
\]
A matrix Waring problem for cubes can be reduced to linear forms in X with arithmetic coefficients. Choose idempotents or involutions whose cubes are themselves; use shifts with prescribed trace to control a third cube. The source reports an exact three-cube classification (LEAN) and local two-cube criteria, but **integer two-cube classification is open** in its programme. Candidate next: higher powers and low-rank matrices under exact trace/conjugacy control.

**Source:** [matrix Cayley-Hamilton method](../sources/math-prime/research/proofs/methods/matrix.cayley_hamilton.md), [result card](../sources/math-prime/research/proofs/theorems/matrix_cubes/README.md).

## 18. Indexed generating products and finite mixed-radix arithmetic

In a commutative ring, indexed subsets (not subsets of distinct *values*) satisfy
\[
\prod_{i=1}^n(u_i+v_i y x^{a_i})
=\sum_{I\subseteq\{1,\ldots,n\}}
\Bigl(\prod_{i\notin I}u_i\Bigr)
\Bigl(\prod_{i\in I}v_i\Bigr)y^{|I|}x^{\sum_{i\in I}a_i}.
\]
For distinct-index subset sums, a nondecreasing positive sequence covers all nonnegative integers iff a₁=1 and a_j≤1+∑_{i<j}a_i. Adding uniqueness forces a_j=2^{j−1}. This is a classical criterion **reproved**, not claimed new. The technique can turn additive representation questions into finite product identities, recurrences and residue constraints.

For radices r_i≥2 and Q_i=∏_{j≤i}r_j, Q=Q_N,
\[
\prod_{i=1}^N\sum_{d=0}^{r_i-1}t^{dQ/Q_i}=\sum_{k=0}^{Q-1}t^k.
\]
This is a finite mixed-radix sum-product collapse.

**Source:** [derived BINGO sum-products](../sources/math-treasure/research/sum-products/README.md).

## 19. Finite-character cancellation on invariant exponent sets (optional extension)

For a finite lattice exponent set A invariant under integral U with U−I invertible over Q and A∩2L empty, the archived cyclovariety proof uses a winding/character sum contradiction to establish a zero on the unit torus, then continuity of the absolute value to attain 1. This may inspire **purely finite** modular-character analogues, but the actual proof uses real topology; it is **not automatically elementary number theory in the strictest sense**. The source's priority is unestablished.

**Source:** [full argument](../sources/math-treasure/research/connections/breadth-screening/cyclovariety/R281-candidate.md).

## 20. Critical synthesis patterns to try next

- **P1:** sum-to-product → Vierzahlensatz → gcd allocation → discriminant square → Pell recurrence → strict squeeze/descent.
- **P2:** quartic Diophantine → Ghost covering → descent to finite seeds → integer Somos recurrence → exact period/drift → marker classification.
- **P3:** norm factorization → construct integer pair/triple multiplication → finite-box Euclidean/Dedekind–Hasse certificate → unit orbit coverage → coefficient lifting → complete Thue solutions.
- **P4:** rank-two Ghost → identify a second independent marker via covering, conjugate condition or transport → Jacobian full-rank lifting on Z² → kernel-lattice sieve → archimedean gap for leftover cosets.
- **P5:** involution dynamics → orbit fundamental arc → finite seeds → invariant gcd label → explicit family/classification rather than bounded search.
- **P6:** finite binomial transform / Bernoulli normal form → forced factor graph → primitive prime allocation → uniform residue obstruction → whole-range completion.
- **P7:** high-degree Thue root → explicit Wronskian-signed Padé ladder → integer denominator normalization → irrationality separation → small-case congruences.
- **P8:** exact finite computations → construct *proof-generating* certificate (rational grid cover, interval inequality, polynomial identity), rather than an incomplete exhaustive search.

**No pattern grants a theorem without its premises and all-index bridge.** The best new work will identify a strong common lemma that makes several of these pipelines provably complete at once.