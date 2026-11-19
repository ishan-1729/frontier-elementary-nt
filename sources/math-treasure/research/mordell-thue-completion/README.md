# Completing the three remaining small Mordell cases

**T026/R65 — independently approved by R73 and R74, 26 October 2026.**
See [review.md](review.md) for the precise audit scope.

For integer X,Y, the remaining equations in the source's nine-case request
have exactly these solutions:

| Equation | All integer solutions (X,Y) |
| --- | --- |
| X³=Y²−4 | (0,±2) |
| X³=Y²−3 | (1,±2) |
| X³=Y²−2 | (−1,±1) |

Thus only the middle equation has strictly positive solutions: (1,2).
These new exclusions combine with the six independently reviewed cases in
[small-mordell](../small-mordell/README.md) to classify the requested nine
members. They do not finish the original
note's comparisons, simultaneous-square question, quartic investigations,
Cauchy substitution or method-transfer requests.

The mechanism is classical: an algebraic norm equation becomes a unit power,
then coefficient congruences exclude every nonzero exponent. All needed unit
generation and exponent arguments are proved below. We do not assume the
maximal order, unique factorization, Dirichlet's unit theorem, a general
Skolem algorithm, or a conditional computer-algebra result for the cubic
orders. The only finite exhaustive computations used in the proof are tiny
*fundamental unit interval* certificates with analytically proved bounds,
not searches for Mordell or Thue solutions. Their complete reproducible
integer certificates are in [the verifier](../../experiments/mordell-thue-completion/verify.py).

## 1. The exact reductions being completed

The R58-reviewed [earlier dossier](../small-mordell/README.md), §§7–8,
proved the following reductions, including all signs. They are repeated to
make the dependency and the exclusions explicit.

For Y²=X³+4, X≥−1. X=−1 is impossible and X=0 gives ±2. If Y is odd,
Y−2 and Y+2 are coprime cubes differing by 4, impossible modulo 9. Thus
Y=2v,X=2u; v²=2u³+1 forces u even. Put X=4w. Then v²=16w³+1,
v≡±1 mod8; change Y's sign to have v≡1 mod8. The coprime integers
(v+1)/2 and (v−1)/8 have product w³, so they equal r³,s³ respectively.
Consequently

    r³−4s³=1,       X=4rs,       Y=16s³+2.                 (1)

Conversely these formulas produce a solution, and both signs of Y cover
the earlier sign choice. Section 5 proves that (r,s)=(1,0) is the sole
integer root of (1).

For d=2,3, Y²=X³+d factors in Z[√d]. The earlier proof established
Euclidean division, unit generation, and coprimality of Y±√d. In detail,
X is odd; for d=2 Y is odd, and for d=3 Y is even with 3∤X. A common
factor's absolute norm divides 4d and |X|³, hence is one. The units are
±e_d^n with e_2=1+√2, e_3=2+√3. Thus

    Y+√d=e_d^r(a+b√d)³,       r=0,1,2.                    (2)

The coefficient equation for r=0 is b(3a²+db²)=1, impossible for both d.
For r=1 the coefficient equations are

    d=2: a³+3a²b+6ab²+2b³=1;
    d=3: a³+6a²b+9ab²+6b³=1.                              (3)

The substitutions t=a+b and t=a+2b, respectively, give

    t³+3tb²−2b³=1;        t³−3tb²+4b³=1.                 (4)

Each equation in (4) has only (t,b)=(1,0), by §§3–4 below.
There is no need to separately solve the two r=2 forms: conjugate (2)
and negate it. Since the square of the norm of e_d is 1,

    −conjugate(e_d² γ³)=e_d(−e_d^(−1) conjugate(γ))³.

The new factor has real part −Y and coefficient of √d equal to 1; its
cube argument remains in Z[√d]. Thus every r=2 solution maps to an r=1
solution with the same X and the opposite Y, and this operation is reversible.
The r=1 root a=1,b=0 gives Y=1,X=−1 for d=2, and Y=2,X=1 for d=3.
Taking the opposite signs of Y proves the table, once (4) is excluded.

## 2. A finite-congruence lemma for every integer exponent

Let R be a commutative ring free of finite rank as a Z-module with a fixed
basis containing 1. Congruences mean coefficientwise divisibility. Let L
be extraction of a basis coefficient other than the constant one. Suppose
η=1+3w is a unit and L(w) is not divisible by 3. Then for every nonzero
integer n,

    v_3(L(η^n))=1+v_3(n).                                 (5)

For positive n put v=v_3(n). The linear term in the binomial expansion has
coefficient 3nL(w), of valuation v+1. For j≥2,

    v_3(binomial(n,j) 3^j) ≥ v+j−v_3(j) ≥ v+2,

because binomial(n,j)=(n/j)binomial(n−1,j−1) and j−v_3(j)≥2.
Coefficients of w^j are integers, so higher terms cannot cancel that linear
term. This proves (5) for positive n without an infinite expansion.

For negative n, first note by induction that
η^(3^h)≡1 mod3^(h+1), h≥0: cubing 1+3^(h+1)z increases the divisibility
by at least one. Inverses preserve these congruences because η is a unit.
Choose h≥v_3(n)+1 and a positive integer m congruent to n modulo 3^h.
Then v_3(m)=v_3(n) and η^m≡η^n mod3^(h+1). The positive-exponent result
therefore gives exactly (5) also for n<0. No convergence or unproved
p-adic zero-counting theorem is hidden here.

If ε³=η, L(ε) and L(ε²) are nonzero modulo3, and L(1)=0, then

    L(ε^n)=0 for n∈Z  if and only if  n=0.                 (6)

Indeed ε³≡1 mod3 forces n≡0 mod3, after which (5) excludes n/3≠0.

## 3. Fundamental unit certificates in two cubic orders

For (C,D)=(3,−2) or (−3,4), let θ be the real root of
f(T)=T³+CT+D and work in R=Z[θ]. Rational-root testing at the divisors
of D shows f is irreducible. Its discriminant −4C³−27D² is −216 or −324,
so it has exactly one real root ρ and two nonreal conjugates z,bar(z).
Every element has unique integer coefficients (a,b,c), and its norm is

    N(a,b,c)=a³−2Ca²c+Cab²+C²ac²+3Dabc−Db³−CDbc²+D²c³.    (7)

This is the determinant of multiplication by a+bθ+cθ². A norm-one
element is a unit: its multiplication matrix has integer inverse, whose
application to 1 supplies the inverse in R. Also its real value U is
positive, since N=U|V|²=1, where V is its nonreal embedding.

The proposed generators and their inverses are

| (C,D) | ε | E=ε^(-1) | ε³ |
| --- | --- | --- | --- |
| (3,−2) | 1−θ−θ² | 17+3θ+5θ² | 25−33θ−15θ² |
| (−3,4) | −7−θ+θ² | 9−11θ+5θ² | −455−3θ+93θ² |

Multiplication modulo f verifies εE=1 and (7) verifies N(ε)=N(E)=1.
The root intervals

    59/100 < ρ < 3/5;
    −2196/1000 < ρ < −2195/1000

follow by opposite signs at the endpoints (uniqueness of the real root
was already established). They imply 0<ε<1 and 1<E<21 or 1<E<58,
respectively. These inequalities are rational interval checks, not decimal
approximations.

Let any norm-one unit have real value U>0. Multiplication by an integer
power of E puts its real value in the half-open interval [1,E). This uses
only E>1 and the Archimedean property. We now exhaust *all* integer
coefficient triples of norm one in that interval.

For the normalized element q(T)=a+bT+cT², |V|=1/√U≤1. Lagrange
interpolation at ρ,z,bar(z), using f(T)/(T−r)=T²+rT+(r²+C), gives

    c= Σ q(r)/f'(r),
    b= Σ r q(r)/f'(r),
    a= Σ (r²+C) q(r)/f'(r).                              (8)

Here the sums run over the three distinct roots. In particular
f'(ρ)=|ρ−z|²>0, |z|²=C+ρ², and
|f'(z)|²=|disc(f)|/f'(ρ). These follow by multiplying root differences.
The following conservative bounds follow from the root intervals:

| Case | |ρ| upper | |z| upper | f'(ρ) lower, upper | |f'(z)| lower | |ρ²+C| upper |
| --- | ---: | ---: | --- | ---: | ---: |
| C=3 | 3/5 | 2 | 4, 102/25 | 7 | 84/25 |
| C=−3 | 11/5 | 7/5 | 11, 288/25 | 5 | 46/25 |

All bounds are strict. Applying the triangle inequality in (8), also using
|z²+C|≤|z|²+|C|, yields

    C=3:  |c| < 21/4+2/7 <6,
          |b| < (3/5)21/4+4/7 <4,
          |a| < (84/25)21/4+2(4+3)/7 <20;

    C=−3: |c| < 58/11+2/5 <6,
          |b| < (11/5)58/11+2(7/5)/5 <13,
          |a| < (46/25)58/11+2(49/25+3)/5 <12.

Hence the entire normalized unit interval lies in these finite integer boxes:

    C=3:  |a|≤19, |b|≤3,  |c|≤5  (3003 triples);
    C=−3: |a|≤11, |b|≤12, |c|≤5  (6325 triples).

Evaluating the explicit integer polynomial (7) on every triple gives:

| Case | Complete list with N=1 inside the stated box |
| --- | --- |
| C=3 | (1,−1,−1), (1,0,0), (17,3,5) |
| C=−3 | (−7,−1,1), (1,0,0), (9,−11,5) |

These lists are ε,1,E in each case. Since ε<1 and E is the excluded
upper endpoint, the only normalized unit is 1. Reversing the normalization
proves that **every norm-one element of R is ε^n for some n∈Z**.
This is a direct generation proof for these orders; no class-number or
maximal-order theorem is being presumed. The certificate's finite loops
are displayed and retained in the verifier; they enumerate every triple
in the mathematically proved boxes, using (7) and an independent determinant.

## 4. Closing both real-quadratic Thue branches

Take L to extract the θ² coefficient. In the first cubic order,

    ε²=5−6θ−4θ²,    (ε³−1)/3=8−11θ−5θ².

In the second,

    ε²=57+4θ−10θ²,  (ε³−1)/3=−152−θ+31θ².

In both cases ε³≡1 mod3, the θ² coefficients of ε and ε² are nonzero
mod3, and L((ε³−1)/3) is a unit modulo3. Sections 2–3 therefore show
that the only norm-one element with θ² coefficient zero is 1.

But N(t−bθ)=t³+Ctb²+Db³. Thus both equations (4) force t−bθ=1,
so t=1,b=0. Section 1 then gives the full classifications for k=−2,−3,
including the opposite-sign r=2 branches. This also proves the complete
integer-root sets of the four cubic forms previously left unresolved:
(1,0), (−1,1), (1,0), (−2,1), in the order of the earlier table.
For r=2 these coordinates follow by using γ=−conjugate(e_d) in the reverse branch transformation; substitution verifies each listed pair. Cubing in
the real quadratic field is injective, so a fixed returned Y determines γ
uniquely and no additional representations of the same point are lost.

## 5. Closing r³−4s³=1 in Z[cuberoot(2)]

Let α=∛2, ε=α−1 and E=1+α+α²=ε^(-1). The polynomial T³−2 is
irreducible by rational-root testing. In R=Z[α],

    N(a+bα+cα²)=a³+2b³+4c³−6abc.

Norm-one elements are units with positive real value, as before. Normalize
any such element to real value U∈[1,E). Since 1<α<13/10, E<4; its
nonreal embedding V has |V|≤1. With ω a primitive cube root of unity,
the elementary three-term Fourier inversion gives

    3a=U+V+bar(V),
    3αb=U+ω²V+ω bar(V),
    3α²c=U+ωV+ω² bar(V).

Each numerator lies strictly between −1 and 6 (the lower bound may be
replaced by ≥−1, which gives the same integer conclusion). Therefore
0≤a,b,c≤1. Evaluating the norm on these eight triples gives only 1 and
E as norm-one elements. E is the excluded endpoint, so every norm-one
element is ε^n, n∈Z. This argument again needs no UFD or maximal-order claim.

Now take L to extract the α coefficient. We have

    ε²=1−2α+α²,       ε³=1+3α−3α².

The hypotheses of (6) hold, with L((ε³−1)/3)=1. Thus the only norm-one
element whose α coefficient vanishes is 1. From r³−4s³=1,
N(r−sα²)=1 and that coefficient is zero, so r−sα²=1, giving r=1,s=0.
Equation (1) and its sign symmetry now give precisely (X,Y)=(0,±2).

## 6. Evidence, failed routes and limits

The [source audit](source-audit.md) distinguishes the owner's request,
the reviewed quadratic reductions, and classical unit/congruence antecedents.
This is a self-contained proof reconstruction with explicit finite unit
certificates, not a claim of a new classification or new general method.

Several routes do not supply the needed exclusions. Real-quadratic unit
reduction alone leaves the four cubic equations, and repeating a bounded
Mordell scan does not resolve them. The available MSE answer for Y²=X³+2
mentions rank and heights but omits the height exclusion; it is not a proof
dependency. The tempting Gaussian/√−2 proof on that page explicitly solves
the opposite-sign equation. General certified Thue software would require
its unconditional guarantees and dependencies to be audited; none is used.
A direct FLT-three citation does not, without a correct bridge, solve these
specific forms; no such unsupported bridge is asserted.

Here the missing step is genuinely discharged: normalized unit intervals
have universal proved bounds, followed by exact exhaustive certificates;
then (5) excludes *all* nonzero integer exponents. The accompanying finite
power tests only guard the formulas and are not substitutes for that lemma.
R73 and R74 independently approved the written proofs and finite certificates
without mathematical corrections. This is not formal verification; see
[review.md](review.md).
