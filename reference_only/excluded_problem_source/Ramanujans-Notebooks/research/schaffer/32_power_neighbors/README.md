# Power sums one above or two below powers of two and three

For \(6\le k\le207\), \(b\in\{2,3\}\), \(a\ge1\) and \(x=b^a+1\) or \(b^a-2\) with \(x>1\), the sum \(S_k(x)=\sum_{j=1}^x j^k\) is never an integer perfect power. Both the endpoint and the target exponent are unrestricted in this statement. It extends [Pro's odd-index square calculation](../pro/PRO8_PRIME_POWER_SEED_RETURN.md) to odd targets and even indices in this range. The argument below gives an effective bound for each higher index on these four families as well. A [cyclotomic factorization condition](cyclotomic.md) applies to arbitrary shifted endpoints, including powers of primes other than 2 and 3.

## Bound the prime target

It suffices to test a prime exponent \(q\). Write \(d=k+1,\ t=b^a,\ G(u)=S_k(1+u)\) as a polynomial with \(G(0)=1\). For \(x=t+1\) put \(u=t\). For \(x=t-2\) put \(u=-t\), so Faulhaber reflection gives \(G(u)=(-1)^{k+1}S_k(x)\). Denote the appropriate sign by \(\eta\). Let \(L\) be the common denominator of the coefficients of \(G\). The binomial Bernoulli formula gives \(g_i=\binom di B_{d-i}(2)/d\) for \(i>0\). The von Staudt–Clausen denominator theorem gives \(L\mid d\prod_{p\le d}p\), in particular \(b^{v_b(L)}\le bd\). Also \(L\mid d!\), by the integer-valued binomial basis of the degree-\(d\) polynomial \(G\). Stage 17 supplies a [simple root of the source polynomial](../17_all_indices/roots.md) for every \(k\ge6\).

All \(t\le(6d)^2\) are tested directly for every target. For larger \(t\), \(G(u)\equiv1\bmod b^{a-v_b(L)}\), a modulus larger than \(t/(3d)-1\) and divisible by 4 for \(b=2\) or 3 for \(b=3\). If \(q=2,\eta=-1\), this is incompatible with \(G(u)=-y^2\). Otherwise take the signed root \(Y=\eta y\) for odd \(q\), or either square root for \(q=2\). If \(q>4d\), raising to \(q\) is a bijection on the units modulo this power of 2 or 3. Hence \(Y\equiv1\) modulo the same power. Since \(y>1\), \(y\ge t/(3d)-1\ge\sqrt t\). But \(S_k(x)<(2t)^d<t^{2d}\). This excludes \(q>4d\).

## Bound the displacement exponent for each remaining prime

Let \(f=G^{1/q}=\sum c_i u^i\), \(c_0=1\). For \(h=\lceil d/q\rceil\), clear the denominators of \(c_0,\ldots,c_h\) by their least common multiple \(B\) and put \(A=B\sum_{i=0}^h c_i u^i\). From \(qG f'=G'f\) one obtains
\[
 qn c_n=\sum_{i=1}^{\min(n,d)}((q+1)i-qn)g_i c_{n-i}.
\]
Some \(c_j\ne0\) for \(h<j\le h+d\). Otherwise the recurrence forces a polynomial root \(f\), impossible at a simple root of \(G\). Choose the first such \(j\).

Set \(s=v_b(L)+2v_b(q)\). The series \(f(b^s v)\) is \(b\)-adically integral. When \(b\nmid q\), binomial coefficients with upper argument \(1/q\in\mathbb Z_b\) are integral by continuity from ordinary integers. When \(b=q\), the denominator of the \(i\)-th binomial coefficient costs at most \(2i\) powers of \(b\), supplied by \(G(b^s v)-1\). At \(a>s\) the series converges to a \(q\)-th root of \(G(u)\) in \(\mathbb Q_b\). For odd prime \(q\) it is the unique such root because neither \(\mathbb Q_2\) nor \(\mathbb Q_3\) has nontrivial odd-order roots of unity. For \(q=2\) choose the matching sign of \(Y\).

Write \(C=B\lceil\|G\|_1\rceil+\|A\|_1\), where the norms sum absolute coefficients. If \(BY-A(u)\ne0\), its absolute value is at most \(Ct^h\). Its \(b\)-adic valuation is at least \(v_b(B)+(h+1)(a-s)\) by the integral tail of the series. Thus
\[
 a\le\lceil\log_b C\rceil +(h+1)s-v_b(B).
\]
If \(BY=A(u)\), the nonzero integral polynomial \(L(A^q-B^qG)\) vanishes at \(u=\pm b^a\). Its first nonzero coefficient, in degree \(j\), is \(-LqB^q c_j\). Dividing out \(u^j\) shows \(a\le v_b(LqB^q c_j)\). The verifier uses the maximum of these two bounds and \(s\).

These reductions also give \(a\ll d^2\log(2d)\) with an effective absolute constant for \(b=2,3\). Indeed \(J=q^2d!\) makes \(f(Jv)\) integral, so \(B\mid J^h\). The binomial representation of \(G\) gives \(\log(2+\|G\|_1)=O(d\log(2d))\). Cauchy's coefficient bound on \(|v|=(4J(1+\|G\|_1))^{-1}\) gives \(\log(2+C)=O((h+1)d\log(2d))\). In the exact-equality case use the ordinary Cauchy root bound for \(L(A^q-B^qG)\), whose logarithmic height is \(O(q(h+1)d\log(2d))\). Here \(qh\le d+q\) and \(q\le4d\).

## Exact verification

Run the Python standard-library checker with PARI/GP on the executable path.

~~~sh
python3 research/schaffer/32_power_neighbors/verify.py --last 207 --gp gp
~~~

For each remaining exponent \(a\) it tests the power-sum prefix at \((b^a+\mathrm{offset})\bmod\ell\) for up to 16 primes \(\ell>k+2\) with \(q\mid\ell-1\). A full length-\(\ell\) block sums to zero since \(\ell-1\nmid k\). Zero and nonzero \(q\)-th powers survive. PARI's exact integer power test checks each survivor in characteristic zero, and tests every exponent on the small endpoints. The Faulhaber polynomials used by PARI are generated independently of the rational recurrence in Python.

~~~text
Last source index 207 prime-target branches 65598
Bounded exponent candidates 60425286 largest exponent bound 78453
Exact small-mark tests 11048 exact modular survivors 261
Nontrivial perfect powers 0
~~~

For mixed products of powers of 2 and 3 see the [unbounded 3-smooth endpoint theorem](smooth.md). A [separate argument in the source exponent](source_indices.md) counts the indices admitting small endpoints across all target exponents.
