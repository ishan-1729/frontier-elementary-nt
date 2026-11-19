# Weighted Vieta descent: signed levels and finite seeds

**Status: checked written proofs, 25 October 2026.** An independent reviewer
read the complete positive and signed-level arguments, checked their critical
steps and exceptional cases, and reran the exact diagnostics. See the
[review receipt](review.md). No formal verification or novelty claim is made.

This dossier recovers and extends Ishan Singh's generalized IMO 1988 Problem 6
manuscripts. The integral involutions, descent, axis families and Chebyshev
parametrization are already in his source, including the weighted `M,N` version.
The present contribution is a uniform proof with a larger sign barrier,
complete finite seeds for every positive right-hand side, and explicit domain
and exceptional-case handling. The coefficient formulation also makes a
normalization already compatible with the owner’s framework explicit.
See [source audit](source-audit.md) for provenance, corrections and antecedents.
The continuation [negative and zero levels](signed-levels.md) completes the
signed-level family under the same integral-reflection hypotheses.

## 1. Exact problem and notation

Fix positive integers `a,b,c,C` such that `a|c` and `b|c`. Put

\[
 m=c/a,\quad \ell=c/b,\quad
 F(x,y)=ax^2+by^2-cxy,\qquad R(u,v)=au^2+cuv+bv^2.
\]

We classify solutions of `F(x,y)=C` with **positive integer coordinates**.
Nonnegative coordinates are used only as intermediate states and seeds.
Neither coprimality of `a,b` nor nonsquareness of a discriminant is assumed.
In particular, `ab|c` is unnecessary: `(a,b,c)=(2,2,2)` is allowed.
This is a normalization convenience, not a substantive novelty claim. If
`d=gcd(a,b)`, any solution requires `d|C`; after dividing all coefficients
and `C` by `d`, the coprime weights `a/d,b/d` have product dividing `c/d`.
Thus the normalized coefficients fit the owner’s product-coefficient shape
(with quotient parameter 1), though his bound on the right-hand side remains
a separate restriction.

The weighted form from the prompt is obtained with `a=1,b=n,c=nh`.
The owner's more general manuscript corresponds to
`a=M,b=N,c=kMNp,C=kq`. Its assumption `q≤MNp` implies `C≤c`.
The distinction between fixed `c,C` and a varying quotient `k` matters: the
finite algorithm below is finite **for each fixed coefficient tuple**, not a
uniform finite search through every possible `k` in this positive-level
setting. The separate [negative-q extension](negative-quotients.md) does
prove such a finite all-positive-quotient search for denominator `abpxy−q`;
that extension has its own [independent review](negative-quotients-review.md).

## 2. Integral reflections and a decreasing height

Define

\[
 X(x,y)=(my-x,y),\qquad Y(x,y)=(x,\ell x-y).
\]

Both maps have integer coefficients, are involutions, and preserve `F`:
substitute `my-x` and use `am=c`, or substitute `ℓx-y` and use `bℓ=c`.
For a solution, writing `x'=my-x,y'=ℓx-y`, we have

\[
 axx'=by^2-C,\qquad byy'=ax^2-C. \tag{1}
\]

They also preserve `gcd(x,y)` because
`gcd(my-x,y)=gcd(x,y)` and `gcd(x,ℓx-y)=gcd(x,y)`.

**Lemma 1.** At every positive solution at least one of `x'<x` or `y'<y`
holds. Indeed, if both inequalities failed, (1) would imply
`by²-C≥ax²` and `ax²-C≥by²`. Adding gives `−2C≥0`, contrary to `C>0`.
A flip that keeps coordinates nonnegative and strictly reduces its changed
coordinate strictly decreases the nonnegative integer height `H=x+y`.
Thus a succession of such flips terminates in at most the original `H` steps.
This is a well-founded descent on the whole pair, not just an uncontrolled
comparison of one coordinate while the other may grow.

## 3. Sharp axis barrier and gcd classification

**Theorem 2.** Suppose `0<C<a+b+c`. Every positive integral solution descends,
through nonnegative integral solutions, to an axis point. Consequently, if
`g=gcd(x,y)`, then

\[
 C=ag^2\quad\hbox{or}\quad C=bg^2. \tag{2}
\]

Conversely, either `C=ar²` or `C=bs²` for a positive integer `r` or `s`
gives a positive solution. Therefore, in this range, positive solvability is
exactly this union of two square classes.

**Proof.** If a flip at a positive point produced a negative coordinate, its
new point would have opposite, nonzero signs. Its unchanged value would be
`R(u,v)` for positive integers `u,v`, and hence at least `a+b+c`. This
contradicts the strict bound. Both flips therefore remain nonnegative.
Lemma 1 gives a strictly decreasing choice whenever both coordinates are
positive. Termination must occur on an axis, since there is a decreasing
choice at every positive point. The endpoint cannot be `(0,0)` because `C>0`.
The preserved gcd makes it `(g,0)` or `(0,g)`, proving (2).
Conversely, `Y(r,0)=(r,ℓr)` and `X(0,s)=(ms,s)` are positive solutions. ∎

The barrier is sharp for this **universal gcd/axis-descent conclusion**.
At `C=a+b+c`, the point `(m+1,1)` is positive, has gcd 1, and satisfies
`F=C`, while `X(m+1,1)=(-1,1)`. Its `C` is neither `a` nor `b`.
It cannot be connected to an axis by these gcd-preserving reflections.
This does not assert that the same `C` has no other, nonprimitive orbit.

For `a=1,b=n,c=nh`, the integer bound becomes `C≤nh+n`, and
`C/g²∈{1,n}`. For the owner’s equation

\[
 x^2+ny^2=k(npxy+q),\qquad q\le np,
\]

put `h=kp,C=kq`. Since `kq≤knp≤knp+n`, all its positive solutions have
`kq=g²` or `kq=ng²`. For `q=1`, `k=g²` or `k=ng²`; at `n=p=q=1`
this recovers the known IMO assertion and the gcd refinement in Ishan's text.
More generally, for the owner’s weighted `M,N` equation the conclusion is
`kq=Mg²` or `kq=Ng²`. A larger sufficient bound, for a fixed candidate `k`, is
`kq<kMNp+M+N`. None of these statements claims this bound is necessary for
an individual solution to descend to an axis.

## 4. All positive C: a genuinely finite boundary reduction

Define the finite set

\[
 \mathcal R_C=\{(u,v)\in\mathbb Z_{\ge0}^2:
                R(u,v)=C\}.
\]

The origin is automatically excluded. Every member satisfies
`u≤⌊√(C/a)⌋` and `v≤⌊√(C/b)⌋`; these bounds follow from nonnegativity
of all three terms **on this quadrant**. We do not call `R` positive definite
on the whole plane: it need not be when `c²≥4ab`.

Form the following seed set `B_C`:

* `(u,0)` for `(u,0)∈R_C`, and `(0,v)` for `(0,v)∈R_C`;
* `(mv+u,v)` and `(u,ℓu+v)` for each `(u,v)∈R_C` with `u,v>0`.

**Theorem 3 (finite seeds).** All seeds solve `F=C`. Every positive solution
is obtained from a seed by a finite word in `X,Y`. Equivalently, the positive
solutions are exactly the positive-coordinate points in the union of these
orbits. In particular,

\[
 F(x,y)=C\text{ has }x,y>0
 \quad\Longleftrightarrow\quad\mathcal R_C\ne\varnothing. \tag{3}
\]

**Proof of completeness.** Start at any positive solution. Choose a strictly
decreasing flip supplied by Lemma 1. If it becomes negative, stop before
executing it. If it remains nonnegative, execute it; if an axis is reached,
stop. Otherwise repeat. The height decreases at each executed step, so this
terminates. Axis endpoints are the first two kinds of seed. If an `X`-flip
would produce `(-u,v)` with `u,v>0`, then `C=R(u,v)` and the preceding
point is `(mv+u,v)`. If a `Y`-flip would produce `(u,-v)`, the preceding
point is `(u,ℓu+v)`. These are exactly the other seed types. Reversing the
executed involutions reconstructs the original solution. No uniqueness of
seeds or disjointness of orbits is needed or asserted.

**Proof of the converse.** The axis seeds visibly solve `F=C`. The mixed
seeds are `X(-u,v)` and `Y(u,-v)`; invariance proves their equation, and
all their coordinates are positive. Axis seeds also yield positive points
by the single flip used in Theorem 2. Every orbit point solves the equation;
retain only its positive-coordinate points. This proves (3). ∎

This is not a circular restatement of the original unbounded positive search.
It replaces that search by the explicitly bounded rectangle above, with
only integer operations. One can scan `0≤v≤⌊√(C/b)⌋` and solve
`au²+cvu+bv²=C` using the exact discriminant

\[
 D_v=(c^2-4ab)v^2+4aC,
 \quad u=(\sqrt{D_v}-cv)/(2a).
\]

Require `D_v≥0`, an integer square root, exact divisibility by `2a`, and
`u≥0`. Check the resulting equation. This is a one-dimensional finite seed
algorithm; it imports no unproved Pell solvability assertion.

**Example.** For `(a,b,c,C)=(1,2,4,7)`, the barrier ends at 6. The mixed
representation is `R(1,1)=7`, yielding positive seeds `(5,1)` and `(1,3)`.
The first has gcd 1 and so cannot lie in an axis orbit (`C∉{1,2}`). This
exhibits the exact new phenomenon beyond the barrier. For `(1,2,4,1)`,
the axis `(1,0)` instead produces `(1,2)`, then `(7,2)`, then `(7,12)`.

## 5. Recurrences, complete parametrization, and exceptional regimes

Use column vectors and compose right to left. The matrices are

\[
 X=\begin{pmatrix}-1&m\\0&1\end{pmatrix},\quad
 Y=\begin{pmatrix}1&0\\\ell&-1\end{pmatrix},\quad
 A=XY=\begin{pmatrix}m\ell-1&-m\\\ell&-1\end{pmatrix}.
\]

Here `det A=1` and `t=tr A=mℓ−2` is an integer. Direct multiplication
shows `A²−tA+I=0`. Define integers `S_0=0,S_1=1` and
`S_{j+1}=tS_j−S_{j-1}` for `j≥1`. Induction, starting with `j=1`, gives

\[
 A^j=S_j A-S_{j-1}I\qquad(j\ge1). \tag{4}
\]

These are `S_j=U_{j-1}(t/2)` if `U_{−1}=0,U_0=1` and
`U_{r+1}=2zU_r−U_{r-1}` define the Chebyshev polynomials. No eigenvalue
division or unproved closed-form identity is used. Negative powers are
integral too: `A^{-1}=tI−A`; use ordinary integer powers of this inverse.
Every coordinate of `A^j z` obeys `w_{j+2}=tw_{j+1}−w_j`, with generating
function `(w_0+(w_1−tw_0)z)/(1−tz+z²)` for nonnegative indices as a formal
power series. Multiplying by the denominator proves that assertion directly.

**Complete parameter domain.** For each `s∈B_C`, take

\[
 A^j s\quad\hbox{and}\quad X A^j s\qquad(j\in\mathbb Z), \tag{5}
\]

and retain precisely the vectors whose **two coordinates are strictly
positive**. Formula (5) is exhaustive: cancellation of adjacent `XX` or `YY`
reduces any word to an alternating word. Even alternating words are `A^j`;
odd ones have the form `XA^j` (use `Y=XA` and `AX=XA^{-1}`). Conversely
every displayed vector lies in a reflection orbit. Duplicate parameter
values may represent the same point; uniqueness is not claimed.
Intermediate points in (5) need not remain positive. This causes no false
solutions because both reflections preserve the equation everywhere, and
Theorem 3 already supplies a nonnegative reconstruction for every target.

Since `m,ℓ≥1`, we have `t≥−1`. For `t=−1,0,1`, respectively, the matrix
identities give `A³=I`, `A²=−I`, and `A³=−I`. Thus it suffices in (5)
to use `0≤j<3,4,6`, respectively, followed by the positivity filter. These
are finite orbits; no infinite family of distinct positive points is implied.
For `t=2`, `(A−I)²=0` and `A^j=I+j(A−I)` for every integer `j`, so use
this linear formula, again with the positivity filter. For `t>2`, (4) and
its inverse remain valid; the two distinct real eigenvalues are an optional
interpretation, not a premise of completeness. The algorithm never divides
by `√(t²−4)`, so the zero-discriminant regime is included.

To recover the owner's first axis family, start with `(r,0)`. For `N≥0`,

\[
 A^N(r,0)=
 \bigl(r[U_N(t/2)+U_{N-1}(t/2)],\;\ell r U_{N-1}(t/2)\bigr).
\]

The second family starts at `X(0,s)=(ms,s)` and gives

\[
 A^N(ms,s)=
 \bigl(ms U_N(t/2),\;s[U_N(t/2)+U_{N-1}(t/2)]\bigr).
\]

These reproduce his unexpanded formulas. By themselves, unrestricted
nonnegative `N` do not guarantee positive coordinates: at `m=ℓ=1`, the
first family starts `(r,0),(0,r),(-r,-r),…`. Formula (5), its reflected
partners, and the explicit filter supply a safe complete specification.

## 6. Reversible Pell-form map and its arithmetic cost

Set `Δ=c²−4ab`. Completing the square gives the exact identity

\[
 (2ax-cy)^2-\Delta y^2=4aF(x,y).
\]

Hence the map `T(x,y)=(u,v)=(2ax−cy,y)` bijects the desired solutions
with integer solutions of

\[
 u^2-\Delta v^2=4aC,\qquad
 v>0,\quad u+cv>0,\quad u+cv\equiv0\pmod{2a}. \tag{6}
\]

The inverse is `x=(u+cv)/(2a),y=v`. The congruence ensures integrality,
the inequalities ensure positivity, and substitution into the identity
returns `F=C`. Thus a solution of the norm equation without the congruence
or positivity requirements is not automatically an original solution.
For `a=1,b=n,c=nh`, `Δ=n²h²−4n` and the congruence is
`u+nhv≡0 (mod 2)`.

The sign of `Δ=ab(t−2)` distinguishes ellipse (`t<2`), degenerate conic
(`t=2`), and hyperbola (`t>2`). The reflection matrices preserve the exact
integer lattice and gcd before this change of coordinates. We do not claim
that an arbitrary unit of the norm equation preserves the congruence in (6),
or that one Pell seed generates every orbit. The finite seeds of Theorem 3
prove the required coverage without either assumption.

For the preceding example `(a,b,c,C)=(1,2,4,1)`, `(7,2)` maps to `(6,2)`:
`6²−8·2²=4`, and the inverse recovers `x=(6+8)/2=7`. The same map with
`C=7` sends `(5,1)` to `(6,1)`, for which `36−8=28`.

## 7. Evidence, scope, and next obligation

Run `python3 experiments/vieta-descent/verify.py`; its deterministic record
is [verification.json](../../experiments/vieta-descent/verification.json).
The finite diagnostics compare seed enumeration to a separate bounded
representation scan, execute/reverse descent on every qualifying point,
check matrix formulas and exceptional regimes, round-trip the Pell map, and
retain sharpness and source-typo witnesses. Their exact ranges are in the
record and script. They test the implementation and boundary handling;
**the universal conclusions depend on the written proofs above**.

The positive-level classification is a useful bounded reduction. The
[continuation](signed-levels.md) proves the corresponding negative-level
finite reduction with unique reduced seeds and classifies the zero level.
Neither part covers nonintegral reflections (`a∤c` or `b∤c`) or establishes
novelty. Canonical representatives without duplicate orbits for the
**positive-level** mixed seeds remain open here. Independent scrutiny of
Theorems 2–5 and the domain of (5) remains the immediate verification
obligation.
