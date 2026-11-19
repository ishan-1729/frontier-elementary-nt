# Primes forced into a variable base at every possible index

Let \(S_k(x)=y^m\), with \(x,y,m\ge2,\ k\ge1\). For every diverging positive function \(F(x)\), the set of endpoints admitting **any** solution with
\[
                            m\ge (k+1)F(x)
\]
has natural density zero. There is no restriction on the parity or size of \(k\), or on the prime factors of \(y\). Quantitatively, for large fixed \(H\),
\[
 \limsup_{X\to\infty}\frac1X\#\{X\le x<2X:\exists k,y,m,\
               S_k(x)=y^m,\ m\ge H(k+1)\}
                        \ \ll\ 1/\log\log H.                         \tag{1}
\]
There is also a restriction on the number of prime factors. If \(h(x)=o(\log\log x)\), outside a density-zero set every perfect-power base has more than \(h(x)\) distinct prime factors **in each of \(x\) and \(x+1\)**. For even \(k\) the same conclusion holds simultaneously for \(2x+1\). These three linear polynomials are pairwise coprime. The conclusions restrict the possible bases even at indices larger than any function prescribed in advance.

The [Carmichael-function argument](carmichael.md) sharpens the ratio bound to \((1+o(1))/2\) for every source and target at all but \(O(X\log\log X/\log X)\) endpoints. There is a [two-sided gap near the diagonal](diagonal.md) at even indices and a [101-digit endpoint certificate](endpoint.md) across both parities. The argument below separately forces large primes of each mark into the base.

## The local implication

Write \(C_p(x)=x-\lfloor x/p\rfloor\). If \(v_p(C_p(x))=1\), put \(d_p=p(p-1)\). Then
\[
                  d_p\nmid k,\qquad
 \prod_{\substack{\ell\mid x(x+1)\\\ell\ {\rm prime},\
                      \ell\equiv1\bmod d_p}}\ell\ \mid y.             \tag{2}
\]
For even \(k\) one may include the primes dividing \(2x+1\) in the product. Indeed if \(d_p\mid k\), each unit contributes 1 modulo \(p^2\), and each multiple of \(p\) contributes zero, since \(k\ge2\). Hence \(S_k(x)\equiv C_p(x)\bmod p^2\) would have valuation one, impossible for a perfect power. For an indicated \(\ell\), the divisibility \(\ell-1\mid k\) would imply \(d_p\mid k\), so a complete block sums to zero modulo \(\ell\). Adding or removing an endpoint multiple leaves the same residue. For \(\ell\mid2x+1\), reflection of the two half blocks gives \(2S_k(x)\equiv0\bmod\ell\) at even \(k\); all such \(\ell\) are odd. This proves (2). It uses the zero-versus-nonzero exponent classes of the [simultaneous residue sieve](../23_endpoint_projection/README.md) and divisibility of the actual base.

For \(x=ap+b,\ 0\le b<p\), \(C_p=a(p-1)+b\). With \(b\) fixed, this is equidistributed modulo \(p^2\) as \(a\) varies modulo \(p^2\). The witness thus occupies the exact fraction \((p-1)/p^2\) of classes modulo \(p^3\). The conditions for different primes are independent in natural density. The density of endpoints with no witness \(p\le z\) is
\[
       \prod_{p\le z}\left(1-\frac{p-1}{p^2}\right)\ll1/\log z ,
\]
by the ordinary Mertens product (the correction to \(\prod(1-1/p)\) converges). In particular almost every endpoint has at least one such witness. A witness alone excludes a class of source indices; the second prime in (2) supplies a divisor of \(y\).

## A large common prime at almost every endpoint

Fix \(z\) and \(H>8\), before taking \(X\to\infty\). For each \(p\le z\) count primes \(\ell\equiv1\bmod d_p\) in \(X^{2/H}<\ell\le X^{1/4}\) dividing \(x\). On \([X,2X)\) the mean of this count tends to
\[
                         \mu_p=\frac{\log(H/8)}{\varphi(d_p)} .
\]
This is Mertens' theorem in a **fixed** progression; no uniformity in its modulus is needed. The variance is at most \(\mu_p+o(1)\). Each single congruence has frequency \(1/\ell+O(1/X)\) and each distinct pair \(1/(\ell_1\ell_2)+O(1/X)\). At most \(X^{1/2}\) ordered pairs occur, so the total rounding error tends to zero. For the mean and variance at \(x+1\) or \(2x+1\), the same congruence counts apply. [Mertens' theorem for Chebotarev sets, specialized to a cyclotomic field](https://arxiv.org/abs/2103.14747)

Chebyshev's inequality bounds the limiting upper fraction having no such prime by \(1/\mu_p\). Take \(z=(\log H)^{1/4}\). The upper density of endpoints for which any of the required prime sets for \(p\le z\) is empty is at most \(O(z^3/\log H)\), since \(\varphi(d_p)\le p^2\). Outside these and the endpoints with no witness, apply (2). In particular \(y\) has a prime factor of \(x\) exceeding \(X^{2/H}\ge x^{1/H}\) on this interval for \(X\ge2\). Since \(y^m<x^{k+1}\), this excludes \(m\ge H(k+1)\). The sum of the two exceptional proportions is \(O(1/\log\log H)\). This proves (1). For diverging \(F\), bound the dyadic fractions first by any fixed \(H\), then let \(H\) tend to infinity. Dyadic summation yields the natural density claim.

## Many distinct common primes

For fixed \(d\), count all odd \(\ell\le X^{1/4}\) with \(\ell\equiv1\bmod d\) dividing one of the three specified marks, considered individually. The same variance argument gives mean \((\log\log X)/\varphi(d)+O_d(1)\) and variance at most the mean plus \(o(1)\). Outside a set of relative density tending to zero each count is at least \((\log\log X)/(2\varphi(d))\). This holds simultaneously for all \(d=d_p\) with \(p\le z\) and all three marks when \(z\) is fixed.

Outside the no-witness set a perfect power must therefore contain, at each mark covered by (2), at least \((\log\log X)/(2z^2)\) distinct primes of that mark in its base. For any \(h(x)=o(\log\log x)\) this exceeds \(h(x)\) eventually. Let \(z\) grow after taking the upper density to obtain the stated density-zero exceptional set. Equivalently, for fixed small \(c>0\), the limiting upper proportion admitting a perfect power with fewer than \(c\log\log x\) common prime divisors at one of the applicable marks is \(O(1/\log(1/c))\), by taking \(z\) a sufficiently small constant times \(c^{-1/2}\). This argument supplies no fixed universal exponent cutoff at every endpoint.
