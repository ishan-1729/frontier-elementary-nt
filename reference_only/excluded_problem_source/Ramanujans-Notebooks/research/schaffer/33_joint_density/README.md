# Uniform joint density in the endpoint and source index

Let \(I,J\) be intervals of positive integers of lengths \(H,K\). For \(x>1\), write \(S_k(x)=\sum_{j=1}^x j^k\). Then
\[
 \sup_{|I|=H,\ |J|=K}
 \frac{\#\{(x,k)\in I\times J:x>1,\ S_k(x)=y^m,\ y,m\ge2\}}{HK}
      \longrightarrow0\qquad(\min(H,K)\longrightarrow\infty).
\]
The two side lengths may grow at unrelated rates and the intervals may be translated independently. There is no cap on the target exponent. The result counts ordered pairs. It leaves open whether one exceptional index might occur at each of many endpoints.

There is a finite congruence form of the result. For every \(\epsilon>0\) a finite collection of modular tests, on the pair \((x,k)\), eliminates all but an \(\epsilon\)-fraction of residue pairs, for all sufficiently large \(k\). No effective rate for the rectangle limit is claimed here. We first exclude individual prime targets using just the residue \(x=2\bmod p\). A separate valuation argument prevents arbitrarily large target primes from escaping the finite tests.

## A shifted multiplicative orbit for every fixed target

Use independent uniform \(X,B\in\widehat{\mathbb Z}=\prod_p\mathbb Z_p\), modeling the endpoint and source exponent. Fix a prime \(q\). At any odd \(p\equiv1\bmod q\) with \(p-1\nmid B\), the complete power-sum block vanishes modulo \(p\). At \(X\equiv2\bmod p\) the source residue is \(1+2^B\), with exponent read modulo \(p-1\). The [moving-prime argument](../08_moving_prime/README.md) also singles out this two-summand endpoint. The complete-orbit projection is analogous to the root-of-unity twists in [Lost Notebook I, printed p. 18](../../../notebooks/lost-1/pages/0029.webp).

For a profinite exponent \(B\) let \(E_p(B)\) indicate that \(p-1\nmid B\) and \(1+2^B\) is neither zero nor a \(q\)-th power mod \(p\). Conditional on any fixed cylinder \(B\equiv r\bmod M\), \(2^B\) is uniform on a coset of \(\langle2^M\rangle\), of size at least \(\mathrm{ord}_p(2)/M\). For a multiplicative coset \(V\) in \(\mathbb F_p^\times\) and a nontrivial character \(\chi\),
\[
                    \left|\sum_{v\in V}\chi(1+v)\right|\le\sqrt p .
\]
Indeed expand the indicator of \(V\) in the characters trivial on its subgroup and use the Jacobi sum bound for each term, extending characters by zero at zero. The bound follows from the Gauss sum identity when the product character is nontrivial and direct summation in the exceptional cases. Expanding the \(q\)-power indicator therefore gives conditional frequency \(1-1/q+O(\sqrt p/|V|+M/(p-1))\) for \(E_p\). This includes the zero residue and the excluded exponent class.

For all but a set of primes of convergent reciprocal sum, \(\mathrm{ord}_p(2)\ge\sqrt p\exp((\log p)^\delta)\) for some \(\delta>0\), by Erdős–Murty, Theorem 3. Thus conditional on any fixed cylinder this expectation tends to \(1-1/q\) away from those primes. The reciprocal prime sum \(T(z)=\sum_{p\le z,\ p\equiv1\bmod q}1/p\) diverges. The bounded functions
\[
                  V_z(B)=T(z)^{-1}\sum_{p\le z,\ p\equiv1\bmod q} E_p(B)/p
\]
have average tending to \(1-1/q\) on every cylinder. Approximation by finite unions of cylinders in Haar measure gives \(\int_A V_z\to(1-1/q)\mu(A)\) for every measurable \(A\). Where the unnormalized sum converges, \(V_z\to0\). Dominated convergence gives that this convergence set has measure zero. [Erdős and Murty, *On the Order of a (mod p)*, Theorem 3](https://mast.queensu.ca/~murty/erdos-ram.pdf)

For almost every \(B\), the primes flagged by \(E_p(B)\) consequently have divergent reciprocal sum. Conditional on \(B\), the independent events \(X\equiv2\bmod p\) at these primes occur infinitely often almost surely, by the second Borel–Cantelli lemma. Thus almost every pair has witnesses excluding the specified target \(q\). Taking the intersection of these full-measure sets covers each fixed prime target. To obtain a finite sieve for actual integer pairs we must also bound which targets need checking.

## Almost every pair has a finite positive source valuation

Fix an odd prime \(p\) and a residue \(r\) of the exponent modulo \(p-1\). Let \(x_p,z\in\mathbb Z_p\) be the endpoint and exponent coordinates. On positive integer endpoints \(n\) put
\[
 F_r(n,z)=\sum_{\substack{1\le j\le n\\p\nmid j}}\omega(j)^r
                         \exp(z\log\langle j\rangle),
\]
where \(\omega(j)\) is the Teichmüller unit and \(\langle j\rangle=j/\omega(j)\). This extends continuously in \(n\) to every \(x_p\). For details, for any continuous \(g:\mathbb Z_p\to\mathbb Z_p\), its indefinite sum \(\sum_{j=1}^n g(j)\) is uniformly continuous in \(n\): modulo \(p^A\), \(g\) is constant on residues modulo some \(p^b\), and a full block of length \(p^{A+b}\) sums to zero. Indefinite summation does not increase the sup norm. Apply this to the coefficients of \(z^i\) in the displayed sum, setting each summand to zero at nonunits. Their uniform valuations are at least \(i-v_p(i!)\), tending to infinity. Thus \(F_r(x_p,z)\) is a restricted analytic series in \(z\), continuous jointly in both coordinates.

Choose a positive integer \(j_0\equiv r\bmod p-1\). On \(x_p=pu+a\), \(0\le a<p\), evaluation at \(z=j_0\) is exactly
\[
                 F_r(pu+a,j_0)=S_{j_0}(pu+a)-p^{j_0}S_{j_0}(u).
\]
This follows on nonnegative integer \(u\) by removing multiples of \(p\), then by continuity. The right side is a polynomial with leading coefficient \((p^{j_0+1}-p^{j_0})/(j_0+1)\ne0\), by the finite power-sum formula in [Notebook I, printed p. 151](../../../notebooks/notebook-1/pages/0161.webp). Outside its finitely many zeros on each disk, the analytic function of \(z\) is nonzero as a series, so its zeros have measure zero by Strassmann's theorem. Fubini shows \(F_r(X_p,B_p)\ne0\) almost surely. The other exponent coordinate \(B\bmod(p-1)\) is independent of \(B_p\). The conclusion holds simultaneously at every odd \(p\), by countability.

The expected sum \(\mathbb E\sum_{p:\,p-1\mid B}1/p=\sum_p1/(p(p-1))\) is finite. For almost every \(B\), the remaining primes therefore have divergent reciprocal sum. Conditional on \(B\), infinitely many of these have \(p\mid X\). At each such prime, \(F_r(X_p,B_p)\) is zero modulo \(p\), since full blocks sum to zero at nonzero exponent residue \(r\). Its valuation is finite almost surely by the preceding argument. This supplies a prime at which the modeled source has a positive, finite valuation \(a\). Every possible prime target must divide \(a\). In particular only the finitely many targets \(q\le a\) need witnesses.

## Transfer the certificates to integers

For almost every pair, choose one such valuation and one congruence witness for each prime \(q\le a\). These are finite certificates with open-and-closed conditions on the profinite pair. An exact valuation \(a\) uses modulus \(p^{a+1}\). For integer \(k\ge a+1\) terms divisible by \(p\) vanish modulo that power, and the unit terms agree with \(F_r(x,k)\). Everything else uses finite-field power tests. Thus each certificate eliminates every integer pair in its congruence classes above a finite lower cutoff on \(k\).

A finite subcollection of these cylinders covers measure greater than \(1-\epsilon\), by countability of the clopen basis and continuity of probability from below. Take common periods for both coordinates and the maximum of their cutoffs. Finite unions of periodic classes have their limiting density uniformly in all translated rectangles as both side lengths tend to infinity, by counting complete blocks. The finite range below the source cutoff costs \(O(1/K)\). Sending \(\epsilon\) to zero proves the assertion. The passage from individual targets to a finite certificate uses the positive valuation explicitly; an infinite union of density-zero statements alone would not justify the result.

## A reproducible finite example

Already at primes through 43, the single test \(v_p(S_k(x))=1\) excludes exactly \(870817301648350357177610646589080103423683348079/1119613845054044495712581992014974483798254513500\) of all periodic pairs with \(k\ge2\), approximately **77.77836%**. The [standard-library verifier](verify.py) computes complete prefix sums mod \(p^2\). A block of length \(p^2\) has increment divisible by \(p\), so \(p\) such blocks suffice. For a fixed exponent residue \(k\bmod p(p-1)\), this gives an exact passing fraction over \(x\bmod p^3\). To average the product over \(k\), fix \(k\bmod L\) for \(L=\operatorname{lcm}_{p\le43}(p-1)\). When \(p\nmid L\) the remaining exponent digit modulo \(p\) is independent for each such prime, and when \(p\mid L\) no digit is missing. This preserves the correlations among exponent periods. The finite example uses only valuations equal to one. The proof of the full density limit also allows finite valuations greater than one and uses separate power-residue tests.

~~~sh
python3 research/schaffer/33_joint_density/verify.py
~~~
