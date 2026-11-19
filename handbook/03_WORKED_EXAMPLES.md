# 03 — Worked examples, tricks, small proofs and reconstruction tests

This file supplies **independent small proofs** where feasible and labels **pipeline sketches** where an underlying theorem needs deeper machinery. For complete large proofs, follow the local archived primary sources. Do not silently upgrade a sketch to a solution.

## Example 1 — Explicit Vierzahlensatz parameterization

Take 6×10=3×20. One valid four-number parametrization is
\[
p=3,\ q=2,\ r=1,\ s=10,
\]
since
\[
6=pq,\quad 10=rs,\quad 3=pr,\quad20=qs.
\]
The general positive-integer construction starts p=gcd(A,C), q=A/p, r=C/p for AB=CD; gcd(q,r)=1 follows. The identity qB=rD gives r|B and q|D, so s=B/r=D/q. This proves the positive form of Vierzahlensatz from Euclid's lemma.

**Why this matters:** It transfers *two interacting factorisations at once*. If gcd(C,D)=1 then gcd(p,q)=1; this gives prime-power allocation for A. If B and D are adjacent, gcd(r,q) or gcd(s,...) can give another split. Always check the pair that is actually coprime, rather than assuming all four parameters are pairwise coprime.

**Read:** [Vierzahlensatz](../sources/MathAgent/knowledge/methods/vierzahlensatz.md).

## Example 2 — A research pipeline, explicitly not a complete proof

Suppose you have derived
\[
3d^2(4d^2+1)=c(c+1).
\]
Choose p,q,r,s such that
\[
pq=3d^2,\quad rs=4d^2+1,\quad pr=c,\quad qs=c+1.
\]
The last two relations show gcd(p,q)=1. Thus nonnegative d, after zero/sign handling, has p=u²,q=3v² or p=3u²,q=v², with d=uv up to sign. The adjacent relation reduces to qs−pr=1. One can now combine the other equation rs=4u²v²+1 with these parametric forms, seek a strict squeeze of a discriminant, or use a Bézout substitution.

**Status:** This is a *worked reduction only*. The original archive associates it with a route toward X²+1=Y³; this page **does not assert** that the remaining elimination is completed.

**Read:** [source outline](../sources/MathAgent/knowledge/examples/vierzahlensatz/kieren_mse_x2_plus_1_eq_y3.md).

## Example 3 — A complete elementary discrete-square squeeze

For every integer n≥1,
\[
n^2<n(n+1)<(n+1)^2,
\]
because n(n+1)−n²=n>0, and (n+1)²−n(n+1)=n+1>0. Since no square of an integer can lie strictly between consecutive integer squares, **n(n+1) cannot be a perfect square** for n≥1.

This illustrates the exact three items to export to harder problems: prove that the middle expression is an *integer*, identify *consecutive* feasible perfect powers, and handle equality endpoints. A Pell recurrence can supply nontrivial consecutive candidates where n,n+1 no longer appear visibly.

**Read:** [squeeze method](../sources/MathAgent/knowledge/methods/squeeze_sequence.md).

## Example 4 — A complete weighted Vieta descent step

Consider
\[
x^2+y^2-3xy=1,\quad x,y\ge0.
\]
The maps X(x,y)=(3y−x,y) and Y(x,y)=(x,3x−y) are integral involutions preserving the equation. For instance, (1,3) is a positive solution, because 1+9−9=1. Its Y-jump is (1,0), a nonnegative solution of lower height 1 rather than height 4. Conversely, applying Y to (1,0) recovers (1,3), then applying X creates (8,3), also a solution.

More generally, for a,b,c,C>0 with a|c,b|c, the involutions from the playbook satisfy axx'=by²−C, byy'=ax²−C. At least one coordinate is strictly decreased for every positive solution. If 0<C<a+b+c, a jump cannot leave the nonnegative quadrant through two coordinates of opposite nonzero sign. Therefore every positive solution descends to an axis point by decreasing x+y.

**Note:** This establishes the descent mechanism and seed coverage in that stated range, **not** an unqualified classification for arbitrary coefficients or signs. A reversal can produce many positive points from one seed, often governed by a second-order linear recurrence.

**Read:** [weighted Vieta dossier](../sources/math-treasure/research/vieta-descent/README.md).

## Example 5 — Exact finite subset-product coefficient

Consider
\[
P(x)=\prod_{j=0}^{10}(x^{2^j}+2^j).
\]
The maximum exponent is 2^{11}−1=2047, attained by selecting each x^{2^j}. The coefficient of x^{2012} corresponds to omitted exponents whose sum is 2047−2012=35=1+2+32. Uniqueness of binary expansion implies the omitted indices are exactly j=0,1,5. At those factors select constants 2^0,2^1,2^5. Therefore
\[
[x^{2012}]P(x)=2^{0+1+5}=\boxed{64}.
\]

**Portability:** converting a count/sum to an indexed generating product avoids overcounting equal values; the missing-index complement can transform a large exponent into a very small integer and hence a one-line arithmetic certificate.

**Read:** [sum-products full dossier](../sources/math-treasure/research/sum-products/README.md).

## Example 6 — A complete unit coefficient-lifting instance in Z[sqrt(2)]

Work in the integer pair ring R=Z[√2] with norm N(a+b√2)=a²−2b². The element
\[
\eta=577+408\sqrt2
\]
has norm 577²−2·408²=1 and inverse 577−408√2 in R. Rewrite
\[
\eta=1+3w,\quad w=192+136\sqrt2.
\]
Let L(a+b√2)=b. Since L(w)=136 is not divisible by 3, the elementary finite-binomial argument gives, for **all nonzero integers n**,
\[
\boxed{v_3(L(\eta^n))=1+v_3(n).}
\]
Check n=1: coefficient 408=3·136 has valuation exactly 1. For n>0 the binomial linear term 3n L(w) has valuation 1+v₃(n), while every term j≥2 has at least one extra factor 3; they cannot cancel it. For n<0, \eta^{-1}=1+3w' with L(w')≡−136 (mod3), so the same proof applies to −n.

**Result:** The √2-coefficient of \eta^n is **never zero** for n≠0, without scanning exponent bounds. This is a tiny complete model for the coefficient-exclusion step that finishes some unit-orbit Diophantine problems.

**Read:** [general proof and actual Mordell applications](../sources/math-treasure/research/mordell-thue-completion/README.md).

## Example 7 — Integer triples with norm from a determinant

Let θ³=4θ−2, and represent each element of the integer polynomial quotient ring as (a,b,c)=a+bθ+cθ². Multiplication by θ has columns determined by
\[
\theta\cdot 1=\theta,\quad\theta\cdot\theta=\theta^2,\quad\theta\cdot\theta^2=-2+4\theta.
\]
Thus
\[
M_\theta=\begin{pmatrix}
0&0&-2\\1&0&4\\0&1&0
\end{pmatrix},\qquad N(\theta)=\det M_\theta=-2.
\]
The multiplication map for any α=a+bθ+cθ² is M_α=aI+bM_θ+cM_θ²; define N(α)=det M_α. Since M_{αβ}=M_α M_β by polynomial-reduction algebra, the determinant proves
\[
N(\alpha\beta)=N(\alpha)N(\beta)
\]
by pure finite integer matrix arithmetic. The adjugate gives the explicit formal inverse/norm relation when N(α)≠0.

**Caution:** This does **not** yet prove Euclidean division, UFD, or a complete solution classification. Those require the archive's *independent* exact rational-box norm cover, minimal-norm Bézout and unit descent.

**Read:** [integer triple method](../sources/math-prime/research/proofs/methods/norm.triples.md), [Mordell proof](../sources/math-prime/research/proofs/theorems/mordell_1963/proof.md).

## Example 8 — Cayley–Hamilton turns a cube into a linear expression

Let X be a 2×2 integer matrix, t=tr X and d=det X. Cayley–Hamilton gives X²=tX−dI. Multiply by X:
\[
X^3=tX^2-dX=(t^2-d)X-tdI.
\]
For X= \(\left(\begin{smallmatrix}1&2\\0&0\end{smallmatrix}\right)\), t=1,d=0 and X²=X, so X³=X. This gives infinitely many explicit integer matrix cubes that are **linear** objects. Use them as building blocks rather than solving for an arbitrary cubic polynomial in four matrix entries.

**Read:** [matrix method](../sources/math-prime/research/proofs/methods/matrix.cayley_hamilton.md).

## Example 9 — How a rank-two ghost obstruction arises

Suppose Ghosts have complete orbit coordinates (m,n)∈Z² and the original integer solutions satisfy one condition μ(m,n)=0. A single congruence at p may reduce μ to a scalar function of m,n; its differential is a 1×2 row, which cannot have full column rank 2. Thus the elementary inverse-function/Jacobian zero-isolation used for two markers **cannot** directly certify isolated points on the entire two-index lattice from only that one marker.

This is a **dimension warning**, not a theorem forbidding all possible clever one-equation exclusions. There could be other structure (growth, algebraic dependence, a second covering marker, a global height bound). Proposed repair: add a second independent integer marker ν(m,n)=0 that is *forced by the original problem*, then seek a nonzero mod-p determinant
\[
\det\begin{pmatrix}
\partial_m\mu&\partial_n\mu\\
\partial_m\nu&\partial_n\nu
\end{pmatrix}\not\equiv0\pmod p
\]
on each seed class. The archive proves conditional finite-polynomial zero counts, not a universal existence of ν.

**Read:** [generalization design](../sources/math-prime/research/proofs/GENERALIZATIONS.md), [two-index full theorem](../sources/math-prime/research/proofs/methods/lifting.binomial.multi/proof.md).

## Example 10 — Quotient classification with a genuinely unbounded offset

From the derived positive-offset programme, consider integer
\[
\frac{x^2+y^2}{xy+3}
\]
for positive x,y. The archived complete theorem classifies the **positive integer quotients** as
\[
\{3t^2:t\ge1\}\ \cup\ \{1,5\}.
\]
It records the finite mixed seed tuples (k,u,v)=(1,1,1),(5,1,2),(5,2,1), while the square-multiple family comes from axis and reflection orbits. The full proof handles both orientations, nonprimitive labels and exact rational offsets in a wider parameter range.

**Check one example:** x=y=1 gives 2/(1+3)=1/2, **not** an integer; this illustrates why **seed coordinates in the source's transformed u,v system are not the original x,y values**. Verify the reconstruction map rather than treating a seed label as a solution tuple. The proof and its statement must be read with exact variable names.

**Read:** [classification](../sources/math-treasure/research/positive-offset-classification/proof.md), [README with examples](../sources/math-treasure/research/positive-offset-classification/README.md).

## Example 11 — Complete small Mordell equations, but not a one-line proof

The derived proof package proves, for *all* integers X,Y,
\[
X^3=Y^2-4 \Longleftrightarrow (X,Y)=(0,\pm2),
\]
\[
X^3=Y^2-3 \Longleftrightarrow (X,Y)=(1,\pm2),
\]
\[
X^3=Y^2-2 \Longleftrightarrow (X,Y)=(-1,\pm1).
\]
Its mechanism is coprime cube allocation in an explicitly constructed quadratic pair ring → finitely many unit exponents → coefficient equations → the valuation theorem from Example 6. This is a **complete documented proof** in the copied dossier, not a claim that finite enumeration of X,Y is enough.

**Read:** [complete written argument](../sources/math-treasure/research/mordell-thue-completion/README.md).

## Example 12 — A local-global false implication to avoid

Suppose an integer polynomial f(n) has a root modulo p^k for each k≥1. This alone does **not** imply an integer root, or even a rational root. Example: f(n)=n²−7 and p=3. We have f(1)=−6 divisible by3 and f'(1)=2 not divisible by3, so elementary Hensel lifting constructs a unique compatible 3-adic root at every precision. But no integer n has n²=7, since 2²<7<3².

In higher-dimensional Ghost problems, the same phenomenon appears as a profinite “false zero” of a marker, which a single-prime sieve cannot dismiss. A **global height argument** or a separate modulus/marker is necessary.

**Read:** [multi-unit limits and false zeros](../sources/math-prime/research/proofs/methods/ghost.unit.multi_unit.md).

## Example 13 — Full proof versus proof-pattern use in Ljunggren

The exact classified theorem is
\[
x^2+1=2y^4,\quad x,y>0
\quad\Longleftrightarrow\quad
(x,y)=(1,1)\ \text{or}\ (239,13).
\]
The **Ghost route** starts from primitive nodes r⁴+q²=2s⁴ and proves every node is on a signed orbit. The original equation occurs on a fixed-coordinate slice; a Somos-4 marker and an **index-dependent exact 3-adic** law rule out all unwanted orbit positions. The **independent norm-sign route** reduces to a quartic Thue equation, constructs an integer rational-approximation ladder with an alternating Wronskian sign, obtains a finite bound, then finishes by congruence.

**Warning:** A bare “negative-Pell square term list” with no proof is merely the original hard problem relocated. The authoritative primitive-orbit coverage and all-index exclusion must be checked in the long manuscripts and Lean audit.

**Read:** [theorem and audit](../sources/math-prime/research/proofs/theorems/ljunggren/README.md), [independent norm-sign](../sources/math-prime/research/proofs/theorems/ljunggren/norm_sign/elementary_proof.md), [Ghost proof arithmetic audit](../sources/math-treasure/research/ghost-classification-audit/proof-audit.md).

---

**Instruction to the frontier model:** For any new problem, first derive its **smallest analogous exact worked example**. If you cannot exhibit the defining transformation with verified bidirectional reconstruction on a small case, you do not yet have a proof route. Do not fill the missing step with unbounded search or an imported rank/class theorem.