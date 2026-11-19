# Rational approximations to the square root of the triangular number

Let \(T=x(x+1)/2\). If \(x>1\) and \(S_{17}(x)\) is a square, then for **every** pair of positive integers \(a,b\),
\[
              a<10^6\bigl(b^2+|b^2T-a^2|\bigr)^4.                    \tag{1}
\]
In particular such a source satisfies
\[
        |b^2T-a^2|>100000 \qquad(1\le b\le800,\ a\ge1).             \tag{2}
\]
The finite window (2) uses the project's [index-17 interval theorem](../31_endpoint_sieve/README.md). Inequality (1) holds at all heights including \(x=1\), requires no class-group assumption, and has no coprimality hypothesis. It extends the [integer square approximation](../21_shifted_roots/seventeenth.md) to rational approximants \(a/b\).

There is a concrete branch with a separate short certificate. The only positive endpoint with **both** \(S_{17}(x)\) and \(4T-3\) square is \(x=1\). The equation \(4T-3=a^2\) itself has infinitely many positive endpoints in \(x\equiv1\bmod M\) for every fixed \(M\), at growing distance in \(T\) from the old integer-square windows.

## Clearing the denominators in the square approximation

Use the [septic identity](../18_seventeenth/README.md) \(45S_{17}(x)=T^2P(T)\), with
\[
 P(t)=1280t^7-6720t^6+21120t^5-46880t^4+72912t^3-74220t^2+43404t-10851 .
\]
Write \(u=a/b,\ c=b^2T-a^2,\ v=c/b^2,\ H=1+|v|,\ D=b^2+|c|=b^2H\). Thus \(T=u^2+v\). If the power sum is square, \(z=\sqrt{P(T)/45}\) is a positive rational whose reduced denominator divides 3, since the square of that denominator divides 45.

Consider the following polynomial part of the positive square root at infinity:
\[
\begin{aligned}
 G(u,v)={}&1024u^7+(3584v-2688)u^5+(4480v^2-6720v+4920)u^3\\
          &+(2240v^3-5040v^2+7380v-5837)u .
\end{aligned}
\]
The approximation is \(G/192\). Every monomial \(u^iv^j\) in \(G\) has \(i+2j\le7\), so \(b^7G\) is integral. The difference \(z-G/192\), if nonzero, has absolute value at least \(1/(192b^7)\). Put
\[
 E=4096P(u^2+v)-5G^2=184320\bigl(z^2-(G/192)^2\bigr).
\]
The coefficient of \(u^{2j}\) in \(E\), for \(0\le j\le3\), has \(v\)-degree at most \(7-j\). Its sum of absolute coefficients is, in order, \(1136177152,983571035,436006736,80964032\). All other coefficients vanish. The corresponding coefficient sum of \(G-1024u^7\) is 42889. The [exact polynomial check](verify.py) also verifies the power-sum identity independently from \(B_{18}(x+1)\).

Suppose \(a\ge10^6D^4\). Then \(u\ge10^6 b^7H^4\), so \(u^2\ge256H\). The three lower powers of \(u\) in \(G\) give \(|G-1024u^7|\le42889Hu^5\), whence \(G>768u^7\). Similarly
\[
    |E|\le2636718955\,H^4u^6 .
\]
If \(E\ne0\), division by \(184320(z+G/192)\) gives
\[
 |z-G/192|
 <\frac{2636718955\,H^4}{4\cdot184320\,u}
 <\frac{1}{192b^7}.
\]
This contradicts the common denominator.

If \(E=0\), multiply by \(b^{14}\) and reduce modulo \(a^2\). Every positive power of \(u\) in \(E\) is at least quadratic, and every monomial has weighted degree at most 14. It follows that
\[
                 a^2\mid4096b^{14}P(c/b^2).
\]
The integer on the right is nonzero, since \(P\) is irreducible modulo 23 and has no rational root. Its absolute value is at most \(4096\cdot277387D^7\), using the coefficient sum of \(P\) and \(D\ge b^2,|c|\). This is smaller than \(10^{12}D^8\le a^2\). This proves (1).

For (2) suppose \(|c|\le C=100000\). Since \(x<2\sqrt T\),
\[
        x < 2\left(10^6\frac{(b^2+C)^4}{b}+317\right).
\]
The derivative of \((b^2+C)^4/b\) has the sign of \(7b^2-C\) for positive \(b\), hence its maximum on \([1,800]\) occurs at an endpoint. At both endpoints the displayed upper bound is below \(10^{27}\). The [existing certificate](../31_endpoint_sieve/README.md) leaves only \(x=1\). The proof of (1) is an explicit instance of approximation at infinity on a curve with a rational leading square. [Walsh, *A quantitative version of Runge's theorem*](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/62/2/107516/a-quantitative-version-of-runge-s-theorem-on-diophantine-equations)

## An infinite Pell branch that returns to every seed congruence

Set \(b=2,c=3\). Then \(n=2x+1\) and \(a\) obey \(n^2-2a^2=7\), and (1) gives \(a<2401000000\). The matrix
\[
 B=\begin{pmatrix}3&4\\2&3\end{pmatrix}
\]
preserves this equation and increases each positive coordinate. If \(a\ge6\), then \(4/3<n/a=\sqrt{2+7/a^2}<3/2\). Therefore the inverse step \((3n-4a,3a-2n)\) is positive and reduces \(a\). At \(0\le a\le5\) the only nonnegative solutions with \(n>0\) are \((3,1)\) and \((5,3)\). Iterating these two seeds to the bound tests 25 endpoints. Exact integer square roots of \(T^2P(T)/45\) leave only \(x=1\). This check is independent of the \(10^{27}\) sieve. Since all higher targets at [index 17](../18_seventeenth/README.md) are settled, this also classifies every perfect power on this branch.

The determinant of \(B\) is 1. Its reduction modulo \(2M\) has finite order \(h\), so \(B^{jh}(3,1)\) gives unbounded endpoints with \(x\equiv1\bmod M\). For each odd \(a\ge1\), \(T=(a^2+3)/4\) has distances \((a+1)/2\) and \((a-1)/2\) to the neighboring squares \(((a-1)/2)^2\) and \(((a+1)/2)^2\), respectively. Thus a bounded window around integer squares misses all sufficiently large members. Taking \(M\) large enough retains the same finitely many local tests at which the seed is a unit, as in the [marked local lifting theorem](../06_ramanujan_audit/Ramanujan_connection_and_proof_audit.md).

## Triangular and pentagonal exponents in the notebook

Write \(G_5(j)=j(3j-1)/2\). Berndt's proof of Entry 12(i) in Part III, Chapter 20 selects terms of a theta expansion using the two exponents \(T(n)\) and \(G_5(j)\), via the congruence \(T(n)-G_5(j)\equiv7\bmod17\). The underlying identity is credited there to Rogers. Here the independent algebraic identity
\[
                 4\bigl(T(x)-6G_5(j)\bigr)-1
                        =4T(x)-(6j-1)^2
\]
puts those same two exponents into (2). Hence at any nontrivial square source,
\[
                 T(x)-6G_5(j) \notin [-24999,25000]
                       \qquad\text{for every }j\in\mathbb Z .
\]
The scale 6 makes the square completion integral. The source equation does not establish that some \(j\) attains a value in the interval. [Part III, printed p. 398, physical p. 413, equation (12.6)](../../../notebooks/notebook-3/pages/0413.webp) · [Prior Parts III–V audit](../notebook-methods/deep/regular-3-5/README.md)

At every source index \(k\ge9,\ k\equiv1\bmod4\), the [general inequality](all_indices.md) uses the squarefree part of \((k+1)/2\) to normalize the leading coefficient. It also proves a limit on the method. At index 17, inequality (1) by itself can exclude at most \(O_\epsilon(X^{2/7+\epsilon})\) endpoints through \(X\), even when every denominator is allowed.

The [generic division calculation](division.md) gives a second test of the Ghost Methods transfer, for maps defined on whole marked covers. To reproduce the finite checks, run `python3 research/schaffer/38_ghost_transfer/verify.py` with SymPy installed.
