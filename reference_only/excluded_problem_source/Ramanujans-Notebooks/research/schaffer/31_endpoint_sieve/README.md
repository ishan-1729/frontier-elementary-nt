# The seventeenth power sum through \(10^{27}\)

For every \(2\le x\le10^{27}\), \(S_{17}(x)\) is not a square. Together with the [classification for higher targets](../18_seventeenth/README.md), this excludes every perfect power at index 17 throughout that interval, with no bound on the base. For endpoints above the interval, the individual residue restrictions below continue to hold.

## Conditions from all primes over a rational prime

Set \(T=x(x+1)/2\) and use the [septic field and polynomial](../24_seventeenth_descent/README.md) \(K=\mathbb Q(w)\), \(P(T)=45 S_{17}(x)/T^2\). If \(S_{17}(x)\) is square, \(\gamma=(20T-w)/(20-w)\) is either a global square or defines an everywhere-unramified quadratic extension split at the real place. The [reciprocity argument](../25_seventeenth_norm/README.md) proves this without a class-number assumption. Thus any prime ideal with odd-order ideal class splits in the quadratic extension, since its Frobenius is the image of its class under a homomorphism to a group of order two. In both cases \(\gamma\) is a local square at that prime.

Using the initial 284 elements from the [previous ideal certificate](../30_two_prime_bases/septic.md) and 816 [additional integral-basis elements](additional_elements.json.xz), the verifier checks that **all 1,092 prime ideals** above rational primes \(p\le3000\) have odd-order class, including primes of higher residue degree. It factors every principal ideal and applies the odd-index relation and propagation checks. It then enumerates every prime decomposition up to the bound independently.

For \(7\le p\le3000\), factor \(P\) modulo \(p\) and let \(h\) run over the distinct monic irreducible factors. Each factor gives a prime ideal because the presentation index for \(w\) has no primes above 5 as divisors, and 20 is invertible modulo \(p\). The residue of \(w/20\) is a root \(r\) of \(h\). Since \(P(1)=45\), \(1-r\) is a unit. The residue of \(\gamma\) is \((T-r)/(1-r)\). Its finite-field norm is \(h(T)/h(1)\). A nonzero element of a finite field of odd order is a square precisely when its norm to the prime field is a square, by Euler's criterion. Allow zero as well. The necessary tests are therefore
\[
                   h(T)h(1)\ \in\ \{u^2:u\in\mathbb F_p\}
                    \quad\text{for every } h .
\]
They also apply at repeated factors, as a necessary condition using just the residue. At 137 only four endpoint residues survive. This procedure retains the separate factors even when their product is already a square. At 3 and 5 test \(P(T)=45z^2\) modulo 27 and 25. At 2 use the [necessary condition](../25_seventeenth_norm/README.md) \(T\equiv1\bmod4\). These give \(x\equiv1\bmod3\), \(x\equiv1,3\bmod5\), and \(x\equiv1,6\bmod8\).

## Exact interval enumeration

The [verifier](verify.py) chooses 17 pairwise coprime moduli, including 8, with product
\[
 M=287221301190179857219688574075240>10^{27}.
\]
For each allowed residue at each modulus it forms its Chinese-remainder contribution modulo \(M\). Split the moduli into two groups. Their sorted lists of sums, with a minus sign for the second group, have 3,456,000 and 3,538,944 entries. A positive integer \(x\le10^{27}\) represented by these residues is either \(a-b\) or \(a+M-b\) for one entry of each list. Binary search enumerates all \(b\) in the corresponding closed interval for each \(a\). The upper and lower endpoints are included.

Exactly 42,575,958 candidates pass the first 17 moduli. Testing 1609 leaves 1,112,596. Eleven further primes leave only \(x=1\). The verifier checks the final value by an integer square root. The [two shifted endpoint factors](balanced.md) give an additional unconditional restriction above this interval. Each must split into two coprime divisors both exceeding four million. The possible nonprincipal-ideal case at larger endpoints is still unresolved.

~~~sh
python3 research/schaffer/31_endpoint_sieve/verify.py --gp gp
~~~

The argument uses the local splitting information from [Pro's normalized class-group calculation](../pro/PRO6_PRINCIPAL_PRIMES_AND_GLOBAL_KERNEL.md). The [count with a growing source index](growing_indices.md) treats all target exponents in a different range.
