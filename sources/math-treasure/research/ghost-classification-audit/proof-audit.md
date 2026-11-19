# Arithmetic proof and audit obligations

Notation follows the recovered manuscript where useful. Every rational-map
assertion below is on its stated nonzero domain; identities alone cannot
supply that domain. The independent exact certificate implementation is
[verify.py](../../experiments/ghost-classification-audit/verify.py).

## 1. Signed descent is exhaustive

A node is (r,s,q) with r,s>0, gcd(r,s)=1,
r⁴+q²=2s⁴, and q≡1 (mod4). Original solutions give (1,y,±x), with the
unique sign satisfying this congruence. The original equation implies x and
y are odd: y even makes x²≡−1 mod16 impossible, and y odd makes x odd.
Thus the sign choice exists.

For any node (x,y,z), parity and coprimality imply all three coordinates
are odd and pairwise coprime. An odd common prime in any pair divides the
third; if x were even and y odd, z²≡2 mod16 is impossible. Also z²≡1
mod16, and z≡1 mod4 therefore implies z≡1 mod8. Set

    A=(x²+z)/2, B=(x²−z)/2.

They are coprime signed integers with A²+B²=y⁴, A≡1 mod4 and 4|B.
A common odd prime would divide x and z; parity excludes 2. A zero leg
forces x=y=1 by coprimality, and then z=1. For y>1 neither leg is zero.

For a primitive triple a²+b²=h² with a odd, b even nonzero and h>0,
the positive integers (h+a)/2 and (h−a)/2 are coprime with square product.
Taking square roots and assigning the sign of b to the second root proves

    a=u²−v², b=2uv, h=u²+v²,

with u>0 initially, v signed, gcd(u,v)=1 and opposite parity. This works
for negative a. Simultaneous negation of u,v is also permitted.
Applying it gives A=u²−v², B=2uv and y²=u²+v². Since A≡1 mod4, u is odd
and v even; reduction of y²=u²+v² modulo8 gives 4|v. Moreover

    x²+2v²=(u+v)².

Negate both u,v if needed to have u+v>0. Apply the same parametrization
to y²=u²+v², now allowing signed u, to obtain

    u=k²−ell², v=2k ell, y=k²+ell²,
    k>0, ell≠0, gcd(k,ell)=1.

Separately gcd(x,v)=1, since a common prime would divide u. Factor
((u+v)−x)/2 and ((u+v)+x)/2. They are positive coprime integers, of
product v²/2, one odd. The odd one is t² and the even one 2w²; choose
w's sign to match v. Thus

    x=|t²−2w²|, v=2tw, u+v=t²+2w²,
    t>0 odd, w even, gcd(t,w)=1.

Here w is even because 4|v. This yields u≡1 mod4, so k is odd and ell
even, even if the earlier u was negative. No r≤s or sign restriction on
ell,w has been smuggled in.

Since k ell=tw, let λ=gcd(k,t)>0, k=sλ,t=rλ. Euclid's lemma gives
ell=rμ,w=sμ, where r,s positive odd and coprime, λ odd, μ even nonzero,
and gcd(λ,μ)=1. All four factors are in fact pairwise coprime by the two
previous primitive factorizations. Equating the expressions for u+v gives

    D μ²−2rsλμ+(r²−s²)λ²=0, D=r²+2s².

Modulo λ this implies λ|D. Completing the square yields

    (Dμ−rsλ)²+λ²(r⁴−2s⁴)
       =D[Dμ²−2rsλμ+(r²−s²)λ²]=0.

Consequently q₀=(D/λ)μ−rs is an odd integer with q₀²+r⁴=2s⁴.
Choose q=±q₀≡1 mod4, writing q₀=δq, δ=±1. The parent is a node and

    Dμ=λ(rs+δq), s≤k<k²+ell²=y.

The strict inequality needs only k≥1 and ell≠0. Reconstruction is exact:

    F=r²λ²−2s²μ², H=s²λ²+r²μ²,
    U=s²λ²−r²μ², V=2rsλμ,
    (x,y,z)=(|F|,H,U²−V²−2UV).

The root is uniquely (1,1,1). Strong induction on height therefore covers
all nodes by actual reversed descent steps. It does not claim arbitrary
ascent parameters give primitive nodes or preserve r=1.

## 2. The integer sequence and reduced arithmetic marker

Set b_−2=13,b_−1=b_0=1,b_1=13 and use

    b_(l+2)b_(l−2)=1764b_(l+1)b_(l−1)−3107b_l².       (S)

For four consecutive terms define

    I(a,b,c,d)=a²d²+1764(ac³+b³d)−3107b²c²−422abcd.

Its seed value is zero. With e=(1764bd−3107c²)/a, the exact identities are
I(b,c,d,e)=(e/a)I(a,b,c,d) and, when I=0,

    ae=1764bd−3107c²,
    b²e=422bcd−ad²−1764c³.                          (8)

The verifier also checks the off-invariant identity: the difference
between the second left and right sides is I/a.
Assume a,b,c,d are nonzero integers, individually coprime to 1764, with
adjacent and distance-two pairs coprime. Both ae and b²e are integers;
gcd(a,b²)=1 makes a Bezout combination prove e integral. Modulo each of
2,3,7, equation (8) gives ae≡c², so e is prime to 1764 and nonzero.
If a prime divides e and d, it divides neither b nor c; the two equations
force it to divide both 3107 and 1764, whose gcd is 1. If it divides e and
c, it divides neither b nor d and would divide 1764. Both are impossible.
These are exactly the new adjacent/distance-two conditions. Induction
constructs every forward term and the same reversed rule every backward
term. The symmetric seeds give b_(−l−1)=b_l.

Modulo3 all terms equal1. Modulo8 the period-four residue pattern, starting
at l=−2, is 5,1,1,5. Four successive entries have ae≡4bd+5c²≡1 mod8;
odd residues are self-inverse, so e≡a and the pattern persists in both
directions. Define

    f_l=(9b_l²−b_(l−1)b_(l+1))/4, X_l=f_l/b_l².

The numerator is 4 mod8, so f_l is an odd integer. A common prime of f_l
and b_l would divide an adjacent term, impossible. Reflection preserves
f_l as well as b_l. If X_l=−1/y² for y>0, then (−f_l)y²=b_l² and f_l<0.
A prime dividing −f_l would divide b_l, contradicting coprimality. Hence
f_l=−1 and |b_l|=y. This is the necessary primitive return to the original
slice; unnormalized coordinates would not suffice.

## 3. Rational maps, coefficients and exhaustive orbit coverage

Put P=(−1,1), Q=(9/4,−21/8), N(X,Y)=(X,−Y), on Y²=X³−2X. Define A by
the reflected chord with Q:

    m=(Y+21/8)/(X−9/4),
    A(X,Y)=(m²−X−9/4, m(X−(m²−X−9/4))−Y).

Define D by the reflected tangent followed by the reflected chord with P:

    d=(3X²−2)/(2Y), X₂=d²−2X, Y₂=d(X−X₂)−Y,
    e=(Y₂−1)/(X₂+1),
    D(X,Y)=(e²−X₂+1, e(X₂−(e²−X₂+1))−Y₂).

These are explicit rational maps; no elliptic group, rank or generator
theorem is used. Let L consist of rational points on the cubic with X
3-integral and X≡−1 mod3. Then Y is 3-integral and Y²≡1 mod3. Every node
maps into L via E(r,s,q)=(−r²/s²,rq/s³): reducing the node equation mod3
excludes divisibility of r,s,q by3 using gcd(r,s)=1.

The A denominator 4X−9≡2 mod3 is nonzero. D has nonzero 2Y, and

    X₂+1=(X⁴+4X³+4X²−8X+4)/(4Y²)

has numerator1 mod3. Direct reduction gives A(X,Y)≡(−1,−Y) and
D(X,Y)≡(−1,−1); N plainly preserves L. Thus every composition used below
is defined and stays in L.

The manuscript Appendix A supplies weighted polynomial triples A,B,D
normalizing as (X/Z²,Y/Z³), with B representing A². Our independent
implementation transcribes those formulas and computes in
Z[t,w]/(w²−t³+2t), storing each element as two univariate polynomials.
Multiplication is (a,b)(c,d)=(ac+(t³−2t)bd,ad+bc). Unique monic division
by w²−t³+2t proves that a zero pair is a polynomial identity on the curve;
this is not a numerical test. The eight appendix certificates verify curve
preservation of A,D, both coordinates of A²=B, NANA=id, and DA=BD.
Additional certificates check the correspondence of A with its chord and
D with its tangent/chord formula in both coordinates, and all three
printed explicit curve-factor certificates. Weighted homogeneity extends
the Z=1 checks to nonzero denominators. B's denominator is −Z² mod3
for a 3-unit representative, hence never vanishes on L. Thus

    A⁻¹=NAN, DA=A²D, D(P)=A(P), A(NP)=P.           (9)

The two base substitutions are checked exactly, with
D(P)=A(P)=(−1/169,239/2197). Conjugating NANA=id by N gives the other
inverse composition, so A really is bijective on L.

The descent adapter also requires signs. Write t=r/s,w=q/s² and
h=(t+δw)/(t²+2)=μ/λ; set

    R=t²−2h², H=1+t²h²,
    U₀=1−t²h², V₀=2th, Z₀=U₀²−V₀²−2U₀V₀.

Then w²+t⁴=2 and the child's point is (−R²/H²,|R|Z₀/H³).
Exact reduction modulo w²+t⁴−2, separately for δ=±1, proves

    D(−t²,δtw)=(−R²/H²,−RZ₀/H³).

Here t²+2,H are positive and R≠0 since its scaled integer F is odd.
The input sign is δ and the output sign adjustment is −sign(R). Therefore
every actual child is N^epsilon D N^eta of its parent. No branch is lost.

It remains to identify the actual sequence, not just some orbit. Define
q_l=b_(l−1)b_(l+1)/b_l² and
Y_l=q_l²(q_(l−1)−q_(l+1))/672. From (S) and I=0,

    q_(l−1)q_(l+1)=(1764q_l−3107)/q_l²,
    J(u,v)=u²v²−422uv+1764(u+v)−3107=0,

where u=q_l,v=q_(l+1), both nonzero. In these coordinates

    X=(9−u)/4, Y=(1764u−3107−u²v²)/(672v),
    X'=(9−v)/4, Y'=(u²v²−1764v+3107)/(672u).

The independently checked polynomial identities in the script show
Y²=X³−2X, A's first coordinate equals X' modulo J, and its second equals
Y' exactly. Denominators u,v,X−9/4=−u/4 are nonzero. The seeds give
P_0=(X_0,Y_0)=P; reflection gives P_(−l−1)=NP_l. Since b_l≡1 and f_l≡−1
mod3, X_l is 3-integral with residue−1. The curve equation then proves
Y_l 3-integral, despite the displayed factor672 in its denominator.
Thus P_l∈L and P_(l+1)=A(P_l) for every signed integer l.

Identity (9) yields D(P_l)=P_(2l+1): start at l=0, propagate forward by
DA=A²D, and backward by injectivity of A². Strong induction on node height
now proves every node equals P_l for some integer l, because each actual
ascent is generated by D and N, acting on indices by l↦2l+1 and l↦−l−1.
This does not require an inverse of D or uniqueness of the index.

## 4. Universal exclusion, not a fixed-index sieve

Let G(j)=b_(2j)b_(2j−2)−9b_(2j−1)²=−4f_(2j−1). We prove for every
integer j that G(j)=4 exactly at j=0,1. Because b_l≡1 mod3, the scaled
integers z_l=(b_l−1)/3 exist before any modular reduction.
For four consecutive scaled values a,b,c,d put

    H=1764db+588d+a²−2ac−2a+588b−3107c²−2072c−150.

An exact integer polynomial identity (independently reproduced) rewrites
the recurrence as

    (1+3a)(z_new−(2−a+2c))=3H.                    (10)

It follows by first comparing unscaled expressions, whose difference is
9H, and then cancelling 3 over the integers. No cancellation in a residue
ring is being used.

For an arbitrary positive precision N work in the full polynomial ring
R_N=(Z/3^N Z)[Z₀,Z₁,Z₂,Z₃]. The polynomial
I_N(a)=Σ_(k=0)^(N−1)(−3a)^k satisfies (1+3a)I_N(a)=1 coefficientwise in
this ring. Thus the actual scaled recurrence modulo3^N is the polynomial
map

    V_N(a,b,c,d)=(b,c,d,2−a+2c+3H I_N(a)).

Modulo3 it is the affine map L(a,b,c,d)=(b,c,d,2−a+2c). Direct symbolic
composition gives

    L⁶=(−2a+3c+6,−2b+3d+6,−3a+4c+12,−3b+4d+12).

Hence U_N=V_N⁶ is the identity modulo3 coefficientwise, for arbitrary
variables, not merely along the computed orbit. Let E act on polynomials
by substitution of U_N and Δ=E−id. Expansion of each monomial proves
Δh∈3R_N for every h, and Δ(3^k h)∈3^(k+1)R_N. These statements remain
meaningful even when those ideals are zero; they do not require division
inside R_N.

The observable

    g=(1+3Z₀)(1+3Z₂)−9(1+3Z₁)²
     =−8+3(Z₀+Z₂+3Z₀Z₂−18Z₁−27Z₁²)

has constant term plus a multiple of3. Since Δ kills constants,
Δ^r g∈3^(r+1)R_N for r≥1. At the actual integer state
Z(j)=(z_(2j−2),z_(2j−1),z_(2j),z_(2j+1)), g evaluates to G(j), and U_N
advances this state to Z(j+3) modulo3^N, including negative starting j.

For m>0 with 3^e|m choose N=e+3. The finite binomial theorem for the
commuting additive operators id and Δ gives

    E^m g−g−mΔg=Σ_(r=2)^m binom(m,r)Δ^r g.

For 2≤r≤m, r binom(m,r)=m binom(m−1,r−1) implies
v₃(binom(m,r))≥max(0,e−v₃(r)). Also v₃(r)≤r−2: r=2,3 are direct, and
3^(r−2)>r for r≥4 follows by induction. Therefore

    r+1+max(0,e−v₃(r)) ≥ e+3.

For e≥v₃(r) use r−v₃(r)≥2; otherwise e≤v₃(r)≤r−2 already suffices.
Every summand vanishes modulo3^(e+3), so

    G(j+3m)−G(j) ≡ m(G(j+3)−G(j)) mod 3^(e+3).    (11)

Write a(j)=G(j+3)−G(j). The same divisibility with r=1,N=2 and r=2,N=3
gives 9|a(j) and a(j+3)≡a(j) mod27. All statements here are finite
polynomial statements at an arbitrarily chosen finite precision. Neither
analytic interpolation nor an infinite 3-adic construction is assumed.

The necessary finite seed certificate is small. Exact recurrence modulo27
uses only invertible denominators and yields, for l=−2,...,8,

    b_l ≡ 13,1,1,13,13,22,22,1,1,4,4.

Consequently G(j) for j=0,...,4 is 4,4,7,13,22 modulo27, with G(0)=G(1)=4
as exact integers. Thus a(0)≡9 and a(1)≡18 mod27.

If j≡2 mod3, telescope increments by3 in either direction: G(j)−4≡3
mod9, so it is nonzero. Otherwise j=c+3m with c=0 or1. If m=0 the exact
seeds apply. For m>0, take e=v₃(m) in (11). Then G(j)−4≡m a(c)
mod3^(e+3), whose right side has exact valuation e+2 and is nonzero.
For m<0, put h=−m and e=v₃(h), and apply (11) from j forward to c:

    4−G(j)≡h a(j) mod3^(e+3).

Telescoping a along its residue class gives a(j)≡a(c) mod27. Multiplying
the difference by h makes it divisible by3^(e+3), so replacement is
legitimate at the required precision. The same nonzero residue excludes
j. This proves the all-integer exclusion without a bounded-index premise.

## 5. Return to Ljunggren and the precise index set

Given any positive solution, its signed node is on the orbit. Hence
X_l=−1/y² for some integer l. Section2 forces f_l=−1 and |b_l|=y.
Exactly one of l and −l−1 is odd. Reflection preserves f and b, so choose
that representative and write l=2j−1. Then G(j)=4; Section4 forces j=0
or1, hence l=−1 or1 and y=1 or13. Substitution gives x=1 or239.
Both positive pairs satisfy the equation directly.

In the complete signed sequence the marker indices are precisely
{−2,−1,0,1}: reflecting the two odd representatives adds0 and−2, and the
same argument excludes every other index. This explains the correction
to the sibling README without changing the manuscript's proof.

## 6. What the checks establish

The universal parts are the descent argument, arithmetic induction,
polynomial identities, two-sided orbit induction and finite-difference
proof just given. The retained polynomial checks are exact identities in
specified quotient rings, and verify every Appendix A table row plus
normalization/adaptation formulas and the recurrence/scaling interfaces.
The script separately checks actual integer terms, signed nodes, rational
coordinates and modular jumps over explicit finite ranges. Those checks
would not replace any universal part if the written proof failed.

Only elementary prime factorization, Euclid/Bezout, the proved primitive
triple parametrizations, polynomial arithmetic, and finite induction and
binomial identities are needed. No elliptic group theorem, Laurent
phenomenon theorem, analytic interpolation theorem, approximation result
or reciprocity theorem is imported. The earlier source's attribution and
formal-verification claims remain distinct historical questions.
