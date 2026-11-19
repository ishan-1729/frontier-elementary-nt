# Two principal norms attached to a seventeenth-power sum

A square at index 17 supplies two primitive norms in \(\mathbb Q(\sqrt{-6404582})\). Their prime ideals above the midpoint have opposite orientations. At \(x=86758\), every prime in both norms has the required splitting and each norm passes an order-five class test separately. Their required choices in that same quotient are incompatible. At \(x=141646\), the original norm even has a primitive representation by the principal quadratic form. Its trace must be divisible by a factor depending on \(x\). The corresponding binary form has a primitive solution over every completion and none over the integers.

## The common root and the degree-twelve cofactor

Let \(n=2x+1,\ T=(n^2-1)/8,\ d=2\cdot73\cdot43867,\ s=\sqrt{-d}\). The ring of integers is \(\mathbb Z[s]\). A square \(S_{17}(x)=y^2\) gives a positive integer \(Z=y/T\). Write \(U=128Z,\ Q=S_{17}(x)/T^2\). The identity below defines \(Q\) for all positive endpoints. The integrality of \(Z\) for square sources is [PRO-8, Lemma 1](../pro/PRO8_PRIME_POWER_SEED_RETURN.md).

\[
 45\cdot8192 Q=R(n^2),\qquad
 R(u)=5u^7-245u^6+6645u^5-123485u^4+
 1598807u^3-13618551u^2+67348527u-144103095 .
\]

Set \(r_0=R(0)=-45d/2,\ J(u)=(R(u)-r_0)/u,\ H(u)=((u-1)^2R(u)-r_0)/u\).
The algebraic integers and their positive norms are
\[
 \alpha_0=U+s,\quad N_0=2^{14}Q+d=2n^2J(n^2)/45,\qquad
 \alpha_1=(n^2-1)U+s,\quad N_1=2^{20}S_{17}(x)+d=2n^2H(n^2)/45 .
\]
Each imaginary coordinate equals one. The first formula supplies the degree-twelve cofactor \(J(n^2)\). The second recovers the [degree-sixteen cofactor](../notebook-methods/deep/lost-3-5/cofactor.md). At \(x=193\) every prime in the midpoint and every unramified prime in \(N_1\) splits, yet \(71\mid N_0\) and \((-d/71)=-1\).

In fact \(3\mid n\). Otherwise \(n^2\equiv1,4,7\bmod9\), and \(R\) is 6 modulo 9 at each of these residues, contradicting integral \(Q\). Write \(a=v_3(n)\). Both \(J(0)\) and \(H(0)\) have 3-adic valuation one, giving \(v_3(N_i)=2a-1\). Also \(5\nmid n\) by the [midpoint calculation](../notebook-methods/deep/lost-3-5/midpoint.md). At other primes dividing \(n\), each norm has valuation at least \(2v_p(n)\). Both norms have 2-adic valuation one. Consequently
\[
 B_0=J(n^2)/15,\quad B_1=H(n^2)/15,\qquad N_i=6(n/3)^2 B_i
\]
with \(B_i\) positive odd integers prime to 3.

## A simultaneous class condition

Let \(G\) be the proper form class group of discriminant \(-4d\), identified with the ideal class group. For every split prime \(p\mid N_0N_1\) choose either prime ideal above \(p\), with class \(c_p\), and put \(e_i=v_p(N_i)\). No rational prime divides \(\alpha_i\) in \(\mathbb Z[s]\), since its second coordinate is one. Thus an inert prime cannot divide \(N_i\), and at a split prime the entire exponent belongs to one ideal above \(p\). A ramified prime dividing \(N_i\) has exponent one since the square of its unique prime ideal is \((p)\).

At a shared split prime one has
\[
 N_1-(n^2-1)^2N_0=d n^2(2-n^2).
\]
Hence \(n^2-1\equiv-1\) if \(p\mid n\), and \(n^2-1\equiv+1\) otherwise. Reduction of \(s\) at the ideal of \(\alpha_0\) gives \(-U\bmod p\). The ideal of \(\alpha_1\) is its conjugate in the first case and the same ideal in the second. This determines a single sign \(\epsilon_p\) for the pair. In \(G^2\) the possible contributions are
\[
       (c_p^{e_0},c_p^{\eta_p e_1})^{\epsilon_p},\qquad
       \eta_p=\begin{cases}-1&p\mid n,\\1&p\nmid n.\end{cases}
\]
The same rule allows both choices independently of convention for a prime dividing just one norm. At ramified primes insert the fixed class to the appropriate exponents. A square source requires that products of these contributions reach \((1,1)\), since both ideals are principal. Any homomorphism out of \(G\) preserves this necessary condition.

Part V's treatment of the Sun–Williams coefficients uses a character of the form class group to retain the class represented by a prime. In a single class variable the two orientations are inverse. Here the character lies on \(G\times G\). At a shared split prime its two orientations give one pair of inverse terms, with the relative sign determined by the displayed norm identity. Multiplying two independent single-norm factors would introduce two additional choices at that prime. Equivalently one can keep the signed sums directly in the integral group ring. The projected simultaneous condition says its identity coefficient is positive. This applies the class weighting of Part V to these two actual integer norms. [Part V, printed pp. 338–339, \(F(M,n)\)](../../../notebooks/lost-5/pages/0344.webp) · [Special class-character combinations](../../../notebooks/lost-5/pages/0345.webp)

The reduced form inequalities give \(|G|=1800\). The class \(A=[3,2,2134861]\) has order 150. Its powers with exponents \(150,75,50,30\) are represented respectively by \([1,0,d],[2,0,3202291],[1894,364,3399],[1574,1088,4257]\). These computations require no analytic class-number bound. Since \(1800=72\cdot25\), the 5-primary group is cyclic of order 25. The map \(c\mapsto c^{360}\) has cyclic image of order five generated by \(A^{360}\).

At \(x=86758\), \(n=3\cdot57839\), \(T\equiv1\bmod4\) and \(Q\) is integral. Every odd norm prime splits. In additive coordinates on the order-five quotient, the nonzero local vectors, each with one free sign, are
\[
               (1,-1),\ (-1,1),\ (2,0),\ (2,0),\ (0,2).
\]
The first two come from \(3,57839\). Their sum is \((b,-b)\) with \(b\in\{0,2,-2\}\). To cancel its second coordinate requires \(b=\pm2\), while cancelling its first requires \(b\in\{0,1,-1\}\). Each coordinate can vanish separately. The joint zero is impossible. The verifier proves primality of every factor, composes and reduces the prime forms, and obtains 16 different pairs. This endpoint also passes the ordinary square-residue tests on \(Q\) at \(2,3,5,7,11,13,73,43867\). There are other individual prime tests that detect its nonsquareness. The asserted incompatibility concerns the two complete norm supports and their order-five classes.

One can compute the number of compatible orientations with the fifth-root dissection used in Part III. There the proof of Entry 2.1.1 writes the summation index as \(5m+\nu\), with \(\zeta_5^\nu\) retaining its residue. Here, for the five nonzero vectors \(v_i\) above, let \(L(r,s)=\prod_i(\zeta^{(r,s)\cdot v_i}+\zeta^{-(r,s)\cdot v_i})\). Expanding the product and using \(\sum_{r=0}^4\zeta^{rk}=5[5\mid k]\) counts oriented sums. The two one-coordinate counts \(5^{-1}\sum_rL(r,0)\) and \(5^{-1}\sum_sL(0,s)\) are both 8. The simultaneous count \(25^{-1}\sum_{r,s}L(r,s)\) is zero. This gives a second exact check in \(\mathbb Q(\zeta_5)\). [Part III, printed p. 17, proof of (2.2.1)](../../../notebooks/lost-3/pages/0026.webp)

The square-divisor coefficient from Part IV also explains which local factor belongs in this product. At a split prime and for a character value \(z\), all ideals of norm \(p^e\) contribute \(L_e=\sum_{j=0}^e z^{e-2j}\). An ideal divisible by the rational ideal \((p)\) comes from one of norm \(p^{e-2}\), with the same class. Subtracting this square divisor gives \(L_e-L_{e-2}=z^e+z^{-e}\) for \(e\ge1\), with negative-index terms set to zero. This is the character-weighted primitive count for a single norm. Requiring a *shared* orientation explains why only two terms also occur at a prime common to both norms. [Part IV, printed p. 192, square-divisor coefficient of (8.3.2)](../../../notebooks/lost-4/pages/0206.webp)

## Forced primes of class order divisible by 25

For **every** square source with \(v_3(n)=a\not\equiv3\bmod5\), at least one of the following holds.

* A prime dividing \(n/3^a\) has an ideal above it of class order divisible by 25.
* Each of \(B_0,B_1\) has a split prime divisor of class order divisible by 25 whose exponent in that cofactor is prime to 5.

Indeed, the oriented contribution at 3 has exponent \(2a-1\), which is nonzero modulo 5. If every other midpoint prime satisfies \(25\nmid\mathrm{ord}(c_p)\), its 360th power is trivial. The ramified primes have order at most two. In each norm a nonzero contribution from a prime of \(B_i\) must cancel the contribution at 3. The decomposition \(N_i=6(n/3)^2B_i\) applies even if \(B_i\) and \(n\) have a common factor.

The two forced cofactor primes are **distinct** whenever \(n^2\not\equiv2\bmod89\) and \(n^2\not\equiv2\bmod503\). In fact
\[
 H(u)=(u-1)^2J(u)+r_0(u-2),\qquad J(2)=89\cdot503\cdot1019 .
\]
A split cofactor prime of order divisible by 25 is coprime to \(r_0\), since the only prime factors of \(r_0\) are \(3,5,73,43867\), and 3 does not divide either \(B_i\). If the prime divides both cofactors it must divide \(n^2-2\) and \(J(2)\). The prime 1019 has \((2/1019)=-1\). For example, throughout the unbounded family \(n=3^{1+5020t}\), \(t\ge0\), any square source requires two *different* such split primes. Indeed the residues \(9^a=2\) occur at \(a=8\bmod44\) modulo 89 and \(a=28\bmod251\) modulo 503. The hypothesis concerns the prime classes in the midpoint. The cofactors have unrestricted size and number of prime factors.

## An integral local-to-global failure in the trace condition

Part V, Entry 15.4.1, recovers an actual trace from a representation of a split prime. Its degree-one case gives \(2A\) after choosing the normalized first coordinate of \(A^2+3B^2=p\). Here the real coordinate of \(\alpha_1\) equals \(128(n^2-1)Z\). Keeping that divisibility gives a primitive binary representation with a discriminant growing with \(n\). [Part V, printed p. 343, Entry 15.4.1](../../../notebooks/lost-5/pages/0349.webp) · [Printed p. 344, proof by pairing \(B\) and \(-B\)](../../../notebooks/lost-5/pages/0350.webp)

Put \(m=n^2-1\) and \(D_0=d/2\). Every square source requires
\[
            D_0 V^2+8192 m^2 W^2=N_1/2
\]
with coprime \(V,W\), since \((V,W)=(1,Z)\) is such a pair. At \(x=141646\), \(n=283293\) and \(\gcd(m,D_0)=1\), so this is a primitive quadratic form. The norm \(N_1\) is primitively represented by the maximal-order principal form. Enumeration of **all** primitive representations of \(N_1=A^2+dB^2\) gives one absolute pair
\[
 \begin{aligned}
 |A|&=5048394168549874380469867764962276249996106336220,\\
 |B|&=901606224507814536760687515570528494467249907 .
 \end{aligned}
\]
Here \(128m=10272630252544\) and \(|A|\bmod(128m)=8590850489308\).

Any primitive integral solution \((V,W)\) of the restricted equation would give the maximal-order representation \((128mW,V)\). It is primitive, since a common prime would have to divide \(128m\) and \(V\). Yet \(N_1/2=D_0+8192m^2Q\) is coprime to \(2m\), and reduction of the restricted equation forbids such a prime. The listed real coordinates exclude every global primitive solution. The verifier also solves the restricted form directly and obtains the empty list.

At every odd prime \(p\nmid m\), the displayed maximal-order solution gives a primitive local solution by taking \(W=A/(128m)\) and \(V=B\). At any odd \(p\mid m\), the unit \((N_1/2)/D_0=1+8192m^2Q/D_0\) has a square root over \(\mathbb Z_p\). At 2 this unit is 1 modulo 8 and also has a square root. Use that root for \(V\) and take \(W=0\). A real solution exists by positivity. Thus all local primitive representation tests succeed for this specific form. This is a failure of integral local-to-global representation at an actual endpoint. The unit-coordinate requirement \(V=1\) was already weakened to arbitrary coprime coordinates.

These necessary conditions use the rational quadratic field. The [septic principal branch](../25_seventeenth_norm/README.md) and [normalized septic 2-class obstruction](../pro/PRO6_PRINCIPAL_PRIMES_AND_GLOBAL_KERNEL.md) remain separate from the order-five calculation.

Run `gp -q -f verify.gp` in this directory. The verifier compares both polynomial identities with 1,000 direct sums, counts reduced forms, proves the class order and prime factorizations, and recomputes the signed class products and integral representations. [Verifier](verify.gp) · [Recorded output](verification.txt) · [PARI primitive representations](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Arithmetic_functions.html#qfbsolve)
