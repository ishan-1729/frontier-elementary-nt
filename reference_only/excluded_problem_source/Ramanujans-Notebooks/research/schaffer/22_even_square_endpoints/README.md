# Almost every endpoint excludes all even source indices simultaneously

Let \(A\) be the set of positive integers \(x\) for which \(S_k(x)\) is a square for **at least one** even integer \(k\ge2\). Then, uniformly in the starting point of every interval \(I\) of \(X\ge3\) consecutive positive integers,
\[
                |A\cap I|\ \ll\ \frac{X}{(\log X)^{1/2}} .
\]
In particular almost every endpoint excludes every even source index and every even target exponent at once. This differs in quantifier order from density zero for each fixed \(k\). The proof also describes a restriction on every exceptional endpoint. If \(S_k(x)=y^2\) with \(k\) even, every prime \(p\equiv3\bmod4\) dividing \((x+2)(2x+3)\) occurs there exactly once and satisfies \(p-1\mid k\). Thus the entire part of \((x+2)(2x+3)\) supported on these primes divides \(\prod_{p\equiv3(4),\ p-1\mid k}p\).

## Two rational sources and their reflections

Take a prime \(p\equiv3\bmod4\) dividing \(x+2\), and write \(x=lp-2\). In the first \(lp\) terms the sum in a complete block, modulo \(p\), is 0 if \(p-1\nmid k\) and \(-1\) otherwise. Removing the last two terms subtracts 1 modulo \(p\) because \(k\) is positive and even. Thus
\[
 S_k(x)\equiv
 \begin{cases}
 -1,&p-1\nmid k,\\
 -1-l,&p-1\mid k
 \end{cases}\pmod p.
\]
The first alternative is a nonsquare. If \(p^2\mid x+2\), the second is a nonsquare too. More generally, define the forbidden residue set modulo \(p^2\) by
\[
 \Omega_p=\{pl-2\bmod p^2:\ l\bmod p,
                  (\tfrac{-1-l}{p})=-1\}.
\]
Every member rules out a square for every even \(k\), regardless of its size. There are \((p-1)/2\) such residues. For example \(p=3\) supplies \(x\equiv7\bmod9\), and \(p=7\) supplies three residues modulo 49. The true cannonball endpoint \(24\) and the endpoint 1 satisfy the necessary conditions.

The second mark is \(2x+3=pt\), where \(t\) is odd. If \(p-1\nmid k\), the half-block sum through \((p-1)/2\) vanishes modulo \(p\), by pairing each unit with its negative. Stopping one term earlier gives \(-2^{-k}\). The preceding complete blocks vanish. If \(p-1\mid k\), all units contribute 1. There are \((t-1)/2\) complete blocks followed by \((p-3)/2\) terms. Hence
\[
 S_k(x)\equiv
 \begin{cases}
 -2^{-k},&p-1\nmid k,\\
 -1-t/2,&p-1\mid k
 \end{cases}\pmod p.
\]
The first is a nonsquare because \(k\) is even. Thus \(p-1\mid k\) and \(p^2\nmid 2x+3\) at any square source. Add to \(\Omega_p\) the \((p-1)/2\) classes for which \((\tfrac{-1-t/2}{p})=-1\). The value of \(t\bmod p\) runs bijectively over the field as \(x\) varies over the lifts of \(2x+3=0\bmod p\). These classes are disjoint from the first family already modulo \(p\), since \(2(x+2)-(2x+3)=1\). The combined forbidden set has **\(p-1\)** residues modulo \(p^2\).

The two signs have an exact interpretation over \(\mathbb Q\). For even \(k\), \(S_k(1)=1\), \(S_k(-1/2)=0\), and \(S_k(1/2)=2^{-k}\) by the difference identity. The latter is itself a rational square. Reflection \(x\mapsto-1-x\) changes the sign of \(S_k\), giving \(-1\) and \(-2^{-k}\) at the two reflected rational sources. The block computations retain the correction even at the primes where these polynomial evaluations cannot simply be reduced.

The complete orbit in the block is the same finite power sum used in the [moving-prime argument](../08_moving_prime/README.md). The [Ghost factor cover](../03_marked_graphs/Schaffer_marked_graphs.md) tracks the endpoints \(x,x+1\); reflection of the known square \(S_k(1)=1\) contributes the additional shift \(x+2\), and the rational source at one half gives the coprime mark \(2x+3\).

## Density and quantitative sieve

Write \(g_p=|\Omega_p|/p^2\) and \(h_p=g_p/(1-g_p)\), using only primes \(3\bmod4\). For finitely many primes the forbidden events are independent by the Chinese remainder theorem. Since \(\sum g_p=\infty\), the product of the surviving proportions tends to zero. This already proves density zero for \(A\), including relative to any fixed progression after omitting primes dividing its modulus. For the logarithmic rate we use the analytic large sieve and the Mertens product in a fixed progression. [Montgomery and Vaughan, *The large sieve*](https://doi.org/10.1112/S0025579300004708) · [Williams, *Mertens' Theorem for Arithmetic Progressions*](https://people.math.carleton.ca/~williams/papers/pdf/057.pdf)

Here are the details for the square moduli. Let \(E=A\cap I\). The function \(f_p=\mathbf1_{\Omega_p^c}/(1-g_p)-1\) has mean zero, mean square \(h_p\), and value \(h_p\) on \(E\). For each squarefree \(d\) composed of our primes take \(f_d=\prod_{p\mid d} f_p,\ h_d=\prod_{p\mid d}h_p\). Its Fourier coefficients on \(\mathbb Z/d^2\) have squared norm \(h_d\) with the expectation normalization. Each nonzero coefficient has a frequency whose reduced denominator has prime support exactly \(d\), by the local mean-zero property. Summing \(f_d\) over \(E\), then applying Cauchy–Schwarz to its Fourier expansion, gives
\[
 |E|^2h_d\ \le\ \sum_{\theta\in\mathrm{supp}(\widehat f_d)}
        \left|\sum_{n\in E}e^{2\pi i n\theta}\right|^2 .
\]
For \(d=1\) use the constant function and the zero frequency. The supports for distinct \(d\) are disjoint. For \(d\le Q\) all denominators are at most \(Q^2\), with spacing at least \(Q^{-4}\) on the circle. The analytic large sieve now gives \(|E|\sum_{d\le Q}h_d\le X+Q^4\), where only the specified squarefree \(d\)'s are summed.

To bound this sum, temporarily form all products of primes up to \(z=Q^c\), weighted by \(h_d/\prod_{p\le z}(1+h_p)\). The expected \(\log d\) is \(\sum_{p\le z}g_p\log p=O(\log z)\). For small fixed \(c>0\), at least half the mass has \(d\le Q\) by Markov's inequality. Moreover
\[
 \log\prod_{p\le z}(1+h_p)
 =\sum_{p\le z}g_p+O(1)=\tfrac12\log\log z+O(1).
\]
Here \(g_p=1/p+O(p^{-2})\) and the sum is restricted to \(p\equiv3\bmod4\). Taking \(Q=\lfloor X^{1/4}\rfloor\) proves the stated interval bound, with a larger constant for small \(X\).

## How many individual witnesses occur

For a uniformly chosen \(1\le x\le X\), let \(W(x)\) be the number of primes whose forbidden set contains \(x\). Its mean is \(\tfrac12\log\log X+O(1)\), and its variance is \(O(\log\log X)\). To check this, first count the primes \(p\le X^{1/10}\). Their mean differs from \(\sum g_p\) by at most \(O(X^{-1}\sum p)\). For distinct \(p,r\) the joint probability is \(g_p g_r+O(pr/X)\) by counting classes modulo \(p^2r^2\). Summing these errors gives \(o(1)\), even before restricting to primes. The variance for the truncated count is therefore at most its mean plus \(o(1)\). Every omitted witness divides \((x+2)(2x+3)\) and is greater than \(X^{1/10}\), so the tail count is bounded by an absolute constant. This also gives the stated variance for the full count, by Cauchy–Schwarz. Consequently \(W(x)/\log\log X\) converges in probability to \(1/2\). Almost every endpoint has arbitrarily many distinct single-prime certificates in this density sense.

The [finite checker](verify.py) independently compares the reflection formula to direct sums over a full set of even exponent residues and all endpoint residues for small primes. Positive-density statements about a fixed finite sieve and the asymptotic estimates above have different scopes. For example, at \(x=49\) the marks are 51 and 101; the sole prime \(3\bmod4\) to test is 3, and the first mark gives the allowed residue 0 at that prime. This proof leaves those endpoints available for other methods.

The [odd-target result](odd_targets.md) uses the two half-integer marks to exclude, for almost every endpoint, all even source exponents that are not multiples of a prescribed odd target exponent.

Combining both theorems with the two stationary valuations gives a [simultaneous result over every target exponent](all_targets.md): more than 63.38 percent of endpoints in lower density exclude all perfect powers at all even source indices.

Counting complete endpoint residue sets at primes with a shared source-exponent period improves the interval bound to [\(O_\epsilon(X^{7/8+\epsilon})\)](../23_endpoint_projection/even_squares.md). A [finite projection of all exponents](../23_endpoint_projection/README.md) also raises the excluded lower density for all targets at even indices to over two thirds.
