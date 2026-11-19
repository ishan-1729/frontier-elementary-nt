# Negative and zero levels of the same integral reflections

**Status: checked written proofs, 25 October 2026; see [review](review.md).** This
continues the positive-level [dossier](README.md), retaining its notation
`F(x,y)=ax²+by²−cxy`, positive integers `a,b,c` with `a|c,b|c`, and positive
integer coordinates. Put `m=c/a,ℓ=c/b,Δ=c²−4ab,d=mℓ=c²/(ab)`.
The involutions are `X(x,y)=(my−x,y)` and `Y(x,y)=(x,ℓx−y)`.
The same matrix `A=XY` and recurrence formulas remain valid.

The negative-level investigation is already present in Ishan Singh's
*Descent Trees in the Negative q Case* and its modified version, which end
with an unresolved descent obligation in Case 2. The finite reduction below
is a repository completion of that obligation for each **fixed** quotient
and coefficient tuple. The later [varying-quotient draft](negative-quotients.md)
adds a uniform finite search for all positive quotients in the owner’s
negative-q family, including q>abp. Classical quadratic-form reduction is
a relevant antecedent; no novelty claim is made.

## 8. Negative levels: positivity and reduced points

Fix a positive integer `N` and ask for `F(x,y)=−N`.
If `Δ≤0`, there are no positive solutions, since

\[
 F(x,y)=a\left(x-\frac{c}{2a}y\right)^2
          +\frac{4ab-c^2}{4a}y^2\ge0.
\]

Suppose henceforth `Δ>0`. Both complementary roots are strictly positive
integers at every positive solution, because

\[
 ax(my-x)=by^2+N>0,\qquad
 by(\ell x-y)=ax^2+N>0. \tag{7}
\]

Call a solution **reduced** when

\[
 2ax\le cy,\qquad 2by\le cx. \tag{8}
\]

These inequalities say precisely that neither reflection decreases
`H=x+y`. If a point is not reduced, apply a decreasing reflection.
Positivity is preserved by (7), and the positive integer height strictly
decreases, so this process terminates at a reduced point. It takes at most
the initial height minus 2 steps. Reflections preserve the gcd as before.
There is no need to force the alternating comparisons `x>y>x'` used in
the source diagrams.

## 9. Explicit seed bounds

For a reduced point let `r=x/y`. Inequalities (8) give

\[
 \frac{2b}{c}\le r\le\frac{c}{2a},\qquad
 \frac{N}{y^2}=cr-ar^2-b=:f(r).
\]

The quadratic `f` is concave. On this closed interval its minimum is at
an endpoint: this follows, for example, by comparing its value with the
linear interpolation between the two endpoints. The endpoint values are

\[
 f(2b/c)=\frac{b\Delta}{c^2},\qquad
 f(c/(2a))=\frac{\Delta}{4a}.
\]

The first is the smaller because `c²>4ab`. Therefore
`N/y²≥bΔ/c²`. Exchanging the roles of `(a,x)` and `(b,y)` proves

\[
 x^2\le\frac{Nc^2}{a\Delta},\qquad
 y^2\le\frac{Nc^2}{b\Delta}. \tag{9}
\]

These are exact rational bounds. No approximation of a square root is
needed for enumeration: use `isqrt((N*c*c)//(a*Δ))` and its `b` analogue.
Because `d=c²/(ab)` is an integer and `d>4`, it is at least 5. Hence

\[
 \frac{c^2}{\Delta}=\frac{d}{d-4}\le5,
 \qquad ax^2\le5N,\quad by^2\le5N. \tag{10}
\]

The uniform constant 5 cannot be decreased over this coefficient family:
`(a,b,c,N,x,y)=(1,5,5,1,2,1)` is reduced and has `by²=5N`.
Exchanging the two coordinates and weights gives equality in `ax²≤5N`.

The integrality of the reflection ratios is essential to this uniform
constant. Without it, a real or rational value `d>4` could be arbitrarily
close to 4; the same argument would retain (9), not imply (10).

Define `S_N` to be all positive integer points in the rectangle (9) which
satisfy `F=−N` and (8). This is a finite set computable by exact arithmetic.
Every positive solution descends to a member of `S_N`; conversely every
seed is already a solution. Thus negative-level positive solvability is
exactly the nonemptiness of this bounded set. This is a finite reduction,
not a test of the original unbounded coordinate range.

## 10. Exactly one reduced seed per orbit

**Theorem 4.** For `Δ>0`, each positive reflection orbit on `F=−N`
contains exactly one member of `S_N`. Every positive solution is obtained
from that seed by `A^j` or `XA^j` for an integer `j`. No positivity filter
is needed at this level: (7) ensures that every word in the reflections
stays positive.

**Proof.** Existence follows from height descent. At most one reflection
can strictly decrease at a positive point. Otherwise `2ax>cy` and
`2by>cx`; multiplying and cancelling the positive `xy` would give
`4ab>c²`, contradicting `Δ>0`.

Consider the undirected graph whose vertices are the positive solutions
in an orbit and whose edges are distinct points related by a reflection.
Ignore a reflection that fixes its point. Each edge changes only one
coordinate, so its endpoints have different heights. Suppose two distinct
reduced vertices were connected. A finite simple path between them exists
by the definition of an orbit (erase repetitions from any connecting word).
Its endpoints cannot have maximal height along this path: their adjacent
vertices have strictly larger height by reducedness and absence of loops.
A vertex of maximal height therefore occurs internally, and both its
neighbors have smaller height. This supplies two distinct decreasing
reflection edges at one vertex, contrary to the preceding paragraph.
The reduced vertices must coincide. The word normal form proved in the
positive-level dossier yields the stated orbit parametrization. ∎

The same argument rules out nontrivial cycles: a vertex of maximal height
on a simple cycle would have two decreasing neighbors. Thus “tree” in the
source can be made precise as the reflection graph with loops omitted,
although its degree is at most two. Different reduced seeds cannot describe
the same positive orbit. Parameters within one orbit may still repeat a
point through a stabilizing reflection; the theorem is about the seed,
not a claim of unique words.

Every nonempty negative-level orbit is infinite. At a reduced point at
least one reflection increases height, since two fixed reflections would
force `Δ=0`. A fixed reflection at **any** positive point implies that
point is reduced: for example, `2ax=cy` gives
`cx=c²y/(2a)>2by`, and the other case is symmetric. After one increasing
step the inverse edge decreases. The other reflection cannot decrease
(at most one can), and cannot fix a nonreduced point, so it increases.
This continues indefinitely with strictly increasing height. This statement
is for `C<0`; it does not extend to the finite-orbit `C>0` regimes.

**Source example.** Set the source's parameters
`(M,N,p,q,k)=(1,2,1,1,3)`. Its point `(x,y)=(5,1)` satisfies
`(x²+2y²)/(2xy−1)=27/9=3`. In the present notation
`(a,b,c,N)=(1,2,6,3)` and `Δ=28`. The reflection gives
`X(5,1)=(1,1)`, which is reduced: `2≤6` and `4≤6`.
The bounds (9) allow only `x=y=1`, so this is the sole reduced seed.
The source's strict assertion `x'<y` is equality in this example; height
nevertheless decreases from 6 to 2. This repairs the descent target without
requiring that strict interlacing claim.

## 11. Zero level

**Theorem 5.** Under the same positive coefficient and integral-reflection
hypotheses, positive integer solutions of `F(x,y)=0` exist if and only if
`Δ=0`. In that case every one is

\[
 (x,y)=t\left(\frac{c}{g},\frac{2a}{g}\right),
 \qquad t\in\mathbb Z_{>0},\quad g=\gcd(c,2a). \tag{11}
\]

**Proof.** If `Δ<0`, the completed-square identity in Section 8 is
strictly positive for `y>0`. If `Δ>0`, a positive solution would make
`x/y` a rational root of `ar²−cr+b=0`, so `√Δ/c` would be rational.
Consequently `√(d(d−4))=d√Δ/c` would be rational, where `d=mℓ≥5`
is an integer. A rational square root of an integer is an integer: if
`p/q` is in lowest terms and its square is an integer, `q²|p²` forces
`q=1`. But for every integer `d≥5`,

\[
 (d-3)^2<d(d-4)<(d-2)^2.
\]

The first difference is `2d−9>0`; the second difference is 4. These
are consecutive nonnegative integer squares, so `d(d−4)` cannot be a
square. This contradiction excludes `Δ>0`.

For `Δ=0`, the completed square vanishes exactly on `2ax=cy`.
Dividing by `g=gcd(c,2a)` makes the two coefficients coprime, giving
precisely the positive multiples (11). ∎

Thus the zero level is a primitive ray and all its positive integer
scalings, not a finite set of fixed-level axis or reduced seeds. Both
reflections fix every point of this ray. The matrix formula still holds,
but cannot generate a different scaling. The integrality hypotheses cannot
be dropped: `(a,b,c)=(1,2,3)` has `Δ=1` and positive zero-level points
such as `(1,1)`, while `b∤c`.

## 12. Signed classification and remaining scope

For fixed positive integer coefficients with `a|c,b|c`, the combined dossier
now gives:

* `C>0`: finitely many axis and mixed boundary seeds, an exact bounded
  existence criterion, and a sharp small-level gcd classification;
* `C<0`: no solutions for `Δ≤0`; otherwise finitely many uniquely reduced
  seeds, explicit bounds, and their infinite positive reflection orbits;
* `C=0`: no positive solutions unless `Δ=0`, when (11) is the complete ray.

The reversible map `(x,y)↦(2ax−cy,y)` and its congruence/positivity conditions
apply to every sign of `C` unchanged. These fixed-level results alone do
not classify varying quotients. The subsequent [checked negative-q extension](negative-quotients.md)
closes that obligation for positive `(ax²+by²)/(abpxy−q)` by cancelling k
from a reduced-product bound. It does not supply the same conclusion for
the positive-offset family. Canonical orbit representatives for positive-level
mixed seeds, nonintegral-reflection lattices, and priority of individual
bounds remain outside that extension. Exact diagnostics supplement these proofs; they do not replace
them or constitute formal verification.
