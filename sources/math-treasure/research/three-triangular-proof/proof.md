# Three triangular numbers through Ankeny's geometric argument

T026/R255. Independently approved without corrections in R258, relative to
the explicit standard Dirichlet import below; see [review.md](review.md).
The mathematical proof body is unchanged. This supplies the theorem G
required by the Dolan–Macmillan reduction. It is **not** the owner's requested
descent proof, not a novel theorem, and not a proof of every case of Legendre's
three-square classification.

## Statement, reduction and dependencies

**Theorem.** For every integer n>=0 there are a,b,c>=0 integers with
n=a(a+1)/2+b(b+1)/2+c(c+1)/2.

We prove that every positive integer M=3 modulo 8 is a sum of three integer
squares. Write M=d²m with m squarefree, using the even and odd exponents
in its unique prime factorization. M and d are odd; hence d²=1 modulo 8
and m=3 modulo 8. In particular m>=3 and m is not a square. A representation
of m gives one of M by multiplying each root by d. It therefore suffices
to prove the assertion for squarefree m=3 modulo 8.

**Imported theorem D (Dirichlet).** If positive integers A,B have gcd(A,B)=1,
there are infinitely many primes q with q=A modulo B. Only existence beyond
a finite bound is used. D was actually checked against Kedlaya's author notes
Theorems 1 and 3 and their derivation, using the existing R232 source snapshot;
its analytic L-function nonvanishing foundation is an explicit standard import,
not reproved here. It is not a disguised three-square assertion. No prime-size
estimate or effective complexity assertion is used. [sources.md](sources.md)
records exact reading and limits. All quadratic reciprocity, Minkowski and
two-square ingredients used below are proved in [arithmetic.md](arithmetic.md).

## An explicit auxiliary progression

Set A=(m-1)/2 and B=4m. Since m=3 modulo 8, A=1 modulo 4. Also
2A=m-1 shows gcd(A,m)=1, so gcd(A,4m)=1. The two conditions
q=A modulo m and q=1 modulo 4 are therefore consistent and identify exactly
A modulo 4m: this is the elementary CRT, and in this case the representative
is already explicit. By D choose a prime q>m in this progression.

We have q odd, q=1 modulo 4, q not dividing m, and -2q=1 modulo m.
Thus (-2q/p)=1 at every p dividing m and (-2q/m)=1. Reciprocity and the
supplements from arithmetic.md give

    1=(-2q/m)=(-2/m)(q/m)=(m/q)=(-m/q).

The last equality uses q=1 modulo 4. Choose b with b²=-m modulo q.
A root exists by the last character equality. Choose its representative
0<=b<q, replacing it by b+q if necessary to make it odd; then 0<b<2q.
Since b²+m=0 modulo q and b²+m=1+3=0 modulo 4, there is an integer
h=(b²+m)/(4q)>0. No simultaneous square-root choice modulo prime powers
is needed. In Ankeny's notation t²=-1/(2q) modulo m; our progression makes
t=1 valid, since 2q=-1 modulo m.

## The real lattice and the exact integer length

For x,y,z in Z define

    U=2qx+by,
    R=U+mz,
    S=U/sqrt(2q),
    T=y sqrt(m/(2q)).

The three column vectors for (R,S,T) are
(2q,sqrt(2q),0), (b,b/sqrt(2q),sqrt(m/(2q))), (m,0,0).
Expanding the determinant along the third column gives absolute determinant
m*sqrt(2q)*sqrt(m/(2q))=m^(3/2)>0. Thus this is a genuine full-rank
lattice Lambda with that covolume; nonzero (x,y,z) maps to a nonzero vector.

The open ball R²+S²+T²<2m has volume
(4pi/3)(2m)^(3/2). Dividing by the covolume gives 8sqrt(2)pi/3>8,
using pi>3. Minkowski's proved strict-volume lemma yields a nonzero integer
(x,y,z) whose image lies in the ball.

Put

    V=S²+T²=(U²+my²)/(2q)=2qx²+2bxy+2hy².

The last equality uses b²+m=4qh, so V is an even integer, although S and T
need not be integers. Also R is an integer. Modulo m the denominator 2q is
invertible and

    R²+V = U²+(2q)^(-1)U² = 0 modulo m.

Since this is a positive integer strictly below 2m, it equals m. This is
the first arithmetic return: a small real lattice norm becomes exactly the
specified integer. It is NOT yet three integer squares. Because m is not
a square, V=m-R² cannot be zero; thus 0<V<=m<q. In particular no prime
factor of V equals q. The inequality q>m was chosen to remove that special
prime case from the norm argument; it is legitimate by D.

## The binary norm return: every prime 3 modulo 4 has even valuation

Let p be an odd prime occurring in V to odd exponent. It cannot equal q.
We show p=1 modulo 4.

First suppose p does not divide m. From m=R²+V we have m=R² modulo p,
so (m/p)=1 (R is nonzero modulo p). If (-m/p)=-1, the congruence
U²+my²=0 modulo p forces U=y=0 modulo p: otherwise y would be nonzero
and (U/y)²=-m. More precisely for a nonzero pair (U,y), let e be the
minimum of their p-adic valuations (with valuation of zero taken as infinity).
After removing p^e, at least one coordinate is nonzero modulo p, so the
remaining U²+my² is a p-adic unit. Its valuation is exactly 2e, even.
But U²+my²=2qV has the same odd p-adic valuation as V, contradiction.
Hence (-m/p)=1, and division by (m/p)=1 yields (-1/p)=1, so p=1 modulo 4.
This argument covers one coordinate zero as well as both nonzero; both zero
would make V=0, already excluded.

Next suppose p divides m. Squarefreeness is essential here. The equation
m=R²+V implies p divides R. The identity 2qV=U²+my² then implies p divides
U. Multiply the first equation by 2q and substitute to obtain

    2qm=2qR²+U²+my².

Divide by p, reduce modulo p, and cancel the nonzero m/p; this gives
y²=2q modulo p. In particular (2q/p)=1. But our prime selection also gives
(-2q/p)=1. Multiplicativity again forces (-1/p)=1. This step fails without
squarefreeness, which is why the initial reduction was explicit.

Every odd prime of V with odd exponent is therefore 1 modulo 4. The proved
two-square criterion gives integers v,w with V=v²+w². Thus m=R²+v²+w².
Scaling by d proves the three-square assertion for every original M=3 modulo 8.

## Return to triangular indices, including zero

Take M=8n+3. Modulo 8 a square is 0,1 or 4; the only triple of these
residues summing to 3 is (1,1,1). All three roots are odd and nonzero.
Replacing them by absolute values and setting a=(|dR|-1)/2,
b=(|dv|-1)/2, c=(|dw|-1)/2 gives nonnegative integer indices. Subtracting
three from the sum of squares and dividing by eight proves the theorem.
For n=0, M=3 and the same argument gives three roots of absolute value 1,
so the indices are (0,0,0). There is no exceptional division by zero.

The reverse identity 8[a(a+1)/2]+1=(2a+1)² proves the familiar equivalence.
Existence above is geometric plus Dirichlet, not a well-founded descent in
n. The owner's distinct request for a descent proof remains unresolved.

## A terminating exact algorithm

Factor M by trial division, producing d,m. Test the progression A+4mk in
increasing k until a prime q>m is found, using trial division for primality.
D proves termination, without supplying a running-time bound. Search the
q residue classes for b²=-m; the proved character condition ensures success,
and parity adjustment produces 0<b<2q.

The ball search can be made explicitly finite using integers only. For any
point of its inverse image, T²<2m gives y²<4q. For each such integer y,
U²+my²<4qm bounds U=2qx+by, so there are finitely many x. For each pair,
compute the integer V=(U²+my²)/(2q); its integrality follows from the identity
above. With V<2m, enumerate integers R=U modulo m satisfying
R²<2m-V, and set z=(R-U)/m. Test 0<R²+V<2m. Minkowski proves at least
one such triple exists; congruence proves its squared norm equals m.

Factor V, construct two-square factors as in arithmetic.md and multiply them
using the two-square identity, or simply test 0<=v<=floor(sqrt V) for
V-v² a square: existence was proved before this finite search. Scale and
convert to triangular indices. Thus every search terminates and every output
has an exact arithmetic certificate. An unbounded blind enumeration of
three-square representations is not used to justify the theorem.
