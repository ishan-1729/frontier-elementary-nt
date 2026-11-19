# The thirteenth power sum has no nontrivial perfect powers

**Theorem.** For positive integers \(x,y\) and any integer \(m\ge2\),
\[
                 1^{13}+2^{13}+\cdots+x^{13}=y^m
                 \quad\Longleftrightarrow\quad x=y=1.
\]
The new part is the square case, including unrestricted endpoint height. The [odd-prime argument](../13_cyclotomic_cofactor/README.md) handles every \(m\) divisible by an odd prime. The proof below uses a relative norm equation, exact ideal factorizations, and finite approximations to the logarithms of global units. Either one of the primes 19 and 79 proves the local and uniqueness steps.

## From a square to an integral norm equation

Write \(T=x(x+1)/2\), \(n=2x+1\), \(R(t)=960t^5-2800t^4+4592t^3-4720t^2+2764t-691\). One has \(105S_{13}(x)=T^2R(T)\). Let \(w\) be a root of the monic polynomial \(g(W)=60^5R(W/60)/960\) and put \(K=\mathbb Q(w)\). The [quintic integral descent](../14_square_descent/README.md) proves that any such square gives \(z\in\mathcal O_K\) with \((60-w)z^2=60T-w\). That descent includes an unconditional class-number-one certificate for \(K\). A [second descent on two genus-two quotients](../14_square_descent/paired.md) gives the same square class even for rational triangular points.

Let \(d^2=30(60-w)\), \(L=K(d)\), \(\sigma(d)=-d\) and \(C=2025-d^2=225+30w\). The polynomial \(-30^5g(60-X^2/30)\) for \(d\) is monic irreducible of degree ten. Since \(8T+1=n^2\), our equation becomes
\[
                A=15n+dz\in\mathcal O_L,\qquad A\sigma(A)=C .
\]
At \(n=3,z=1\) put \(A_*=45+d\). We determine all integral \(n\) and integral \(z\in K\) satisfying this norm equation. For a related quadratic norm calculation in the notebook commentary, the equation \(n^2-3m^2=6\) is treated by multiplying one solution by powers of a unit. Here the known solution supplies a base ideal, and the relative units have rank three. [Lost Notebook IV, printed p. 174](../../../notebooks/lost-4/pages/0189.webp) · [Skolem's method for norm-form equations](https://doi.org/10.1016/0022-314X%2888%2990098-4)

## Sixteen possible ideals

Factor \((C)\) in the maximal order of \(L\). A prime fixed by \(\sigma\) has \(v_{\mathfrak p}(A)=v_{\mathfrak p}(C)/2\). For a pair exchanged by \(\sigma\), put \(e=v_{\mathfrak p}(C)\). The two valuations of \(A\) are \(a,e-a\), with \(0\le a\le e\). If unequal, their minimum is the valuation of \(A+\sigma A=30n\), hence is at least \(v_{\mathfrak p}(30)\).

| Rational prime | Valuations of \(C\) at fixed primes | Valuations within exchanged pairs | Possible valuations of \(A\) in each pair |
|---:|---|---|---|
| 3 | 6, 2 | (2,2), (2,2) | (1,1) in both |
| 5 | 6, 8 | (3,3) | (1,2) in either order |
| 7 | | (1,1) | (0,1) in either order |
| 43 | | (1,1) | (0,1) in either order |
| 127 | | (1,1) | (0,1) in either order |

At every other prime the valuation is zero. Starting from \((A_*)\) there are four independent possible interchanges. At each pair choose the prime \(P\) where \(A_*\) has larger valuation. Interchanging the valuations multiplies the ideal by \(\sigma(P)/P\). The [verifier](verify.gp) constructs four factored elements \(b_i\) with exactly these fractional ideals. It factors their constituents and sums the exponents at every prime to verify each equality. For a subset \(I\) of the four interchanges put \(b_I=A_*\prod_{i\in I} b_i\). Then any solution has \(A/b_I\in\mathcal O_L^*\) for at least one of these sixteen choices.

## Three unit directions in a five-dimensional local space

The signatures of \(K,L\) are \((1,2),(2,4)\). Dirichlet's theorem gives absolute unit ranks two and five. The image of the norm on units contains the squares of the units of \(K\), so the norm-one unit group \(U\) has rank three. Since \(L\) has a real embedding, its only torsion units are \(\pm1\).

Fix \(p=19\) or \(79\). The integral polynomial giving the reduced presentation of \(L\) in the verifier is irreducible modulo \(p\). Thus its power basis at \(p\) gives the unramified extension of degree ten, with residue field \(\mathbb F_{p^{10}}\). The fixed subfield has residue field \(\mathbb F_{p^5}\). Put \(E=p^{10}-1\) and for local units define
\[
                  D(a)=p^{-1}\log(a^E).
\]
The logarithm is additive on products and well defined since \(a^E\equiv1\bmod p\). If \(a^E=1+pu\), its reduction modulo \(p^{h-1}\) for \(2\le h\le4\) is the exact finite calculation
\[
              D(a)\equiv\sum_{j=1}^{h-1}(-p)^{j-1}u^j/j\pmod{p^{h-1}}.
\]
Every omitted term has valuation at least \(h-1\), including terms whose denominator is divisible by \(p\).

The verifier constructs units \(e_i\in\mathcal O_L^*\), checks that all their prime ideal valuations vanish, and computes the columns \(D(e_i/\sigma e_i)\). Three are independent modulo \(p\). They belong to \(D(U)\), and their preimages generate a subgroup of index prime to \(p\) in \(U\) modulo torsion. Indeed any dependence in \(U/U^p\) would give a dependence after applying \(D\bmod p\), and \(U/U^p\) has dimension three. It follows that **every** norm-one unit logarithm belongs to the \(\mathbb Z_p\)-span of the three columns. That span, denoted \(W\), is saturated in the local additive lattice because it has a unit three-by-three minor. In particular, a vector of \(W\) divisible by \(p^j\) can be divided by \(p^j\) inside \(W\). The proof needs neither the full class group of \(L\) nor a certified basis of its global units. The factor and rank checks establish all properties used from the proposed generators. [PARI documentation on compact number-field generators](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html)

For an integer candidate \(n\) and a sign of the local square root, let \(u=dz\), so \(u^2=(15n)^2-C\). If \((A)=(b_I)\), the quotient \((A/\sigma A)/(b_I/\sigma b_I)\) is a norm-one unit. Consequently
\[
                  D(A/\sigma A)-D(b_I/\sigma b_I)\in W .
\]
This condition uses the global ideal and unit equations, in addition to local squareness of the power sum.

## Enumerating residue disks and proving uniqueness

At either prime the verifier enumerates every \(n\bmod p^2\). It checks the square of \(z\) in the degree-five subfield; all candidate values are nonzero modulo \(p\). Since \(p^5\equiv3\bmod4\), exponentiation by \((p^5+1)/4\) and a squaring check give the two possible signs when roots exist. Hensel lifting preserves each sign uniquely. At precision \(p^h\), an invertible minor of the three unit logarithms gives a linear projection annihilating their span modulo \(p^{h-1}\). The test above is evaluated for each ideal and sign. Only lifts of survivors need testing at the next precision.

| Prime | Triples (residue, root sign, ideal) at depth two | Depth three | Depth four |
|---:|---:|---:|---:|
| 19 | 12 | 8 | 4 |
| 79 | 28 | 4 | 4 |

At depth four these are exactly the reductions of \((n_0,z_0)=(3,\pm1),(-3,\pm1)\) with their actual ideals. To prove that each disk contains only its known integral point, suppose there is another solution in the same disk. If \(n=n_0\), the given sign fixes \(z\). Otherwise write \(r=v_p(n-n_0)\ge4\), \(u_0=dz_0\), \(A_0=15n_0+u_0\). The identity \((u-u_0)(u+u_0)=225(n-n_0)(n+n_0)\), with \(u+u_0\) a unit, gives
\[
 \frac{A}{A_0}\equiv 1+\frac{15}{u_0}(n-n_0)\pmod{p^{r+1}}.
\]
The matching ideal makes \(A/A_0\) a norm-one unit. Thus \(D((A/A_0)^2)\in W\). Dividing by \(p^{r-1}\) and reducing modulo \(p\) gives
\[
       \frac{30E}{u_0}\,\frac{n-n_0}{p^r}\pmod p\quad\in W\bmod p .
\]
The scalar on the right of the first factor is nonzero modulo \(p\). At each of the four known branches the verifier checks \(30E/u_0\notin W\bmod p\). This contradiction excludes the entire punctured disk at every subsequent depth. It follows that \(n=\pm3\). Since \(n=2x+1\) with \(x>0\), one obtains \(x=1\), as claimed.

Run the [quintic verifier](../14_square_descent/verify.gp) and the [relative-unit verifier](verify.gp) with PARI/GP, for example \(\mathtt{gp\ -q\ -s\ 512M\ verify.gp}\) in the corresponding directory. The [paired descent and endpoint sieve](../14_square_descent/paired.md) supply a separate check through \(10^{20}\). The computer-assisted theorem classifies index 13; Schäffer's conjecture still requires control of unbounded source indices.

A [prime-index extension](prime_indices.md) proves squarefreeness, a nonreal factor and a unit-rank deficit for every prime source index at least eleven.
