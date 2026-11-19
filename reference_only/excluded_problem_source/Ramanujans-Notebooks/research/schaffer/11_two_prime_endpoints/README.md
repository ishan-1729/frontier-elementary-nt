# A global exclusion for endpoints supported on two small primes

**Theorem.** Let \(x>1,\ k\ge1\), and let \(q\) be any odd prime. Remove the largest \(q\)-th power dividing \(x(x+1)\). If the remaining factor has at most two distinct prime divisors, each at most 13, then \(S_k(x)\) is not a \(q\)-th power. Both \(k\) and \(q\) are unbounded. For the marked endpoint cover this tests an intrinsic property of the integer endpoint, so rational rescaling of the roots leaves it invariant.

For \(S_{15}\), the [valuation cover](../09_s15_thue/README.md) makes this applicable to **every** \(d=1,5,7\) branch at every odd prime exponent. Every nontrivial source at an odd prime exponent must therefore lie on \(d=35\). The [bounded prime-exponent result](../10_prime_exponents/README.md) also treats that branch throughout its stated range.

## Reduction to explicit endpoints

The hypothesis is precisely \(x(x+1)=w z^q\) with \(w\) a positive integer supported on a set of two primes at most 13. If fewer are present, extend the support to such a set. Theorem 1.2 of Bennett, Győry, Mignotte and Pintér restricts \(x\) to the 21 explicitly listed values in the verifier. Their Theorem 1.1 further restricts \(x\) to \(1,2,3,4,7,8\) if \(q\ge7\): distribute the prime powers of \(w\) between the coprime endpoints \(A u^q,B v^q\) and order the coefficients. The theorem forces \(u=v=1\) and the smaller endpoint is one of these six integers. Equal coefficients would both be 1 by coprimality, and two distinct positive \(q\)-th powers differ by more than 1. [Bennett, Győry, Mignotte and Pintér, Theorems 1.1 and 1.2](https://personal.math.ubc.ca/~bennett/BGMP.pdf)

For \(x=2\) a putative \(y^q=1+2^k\) gives \(2^k=(y-1)(1+y+\cdots+y^{q-1})\). The second factor is odd and greater than 1. For \(x=3,4,7,8\) the exact 2-adic valuation of \(S_k(x)\) is always 1, 2 or 4. To check this for every \(k\), take residues modulo 32. For \(k\ge5\) even bases contribute zero and every odd base has order dividing 8. The twelve checks \(1\le k\le12\) include the initial values and a complete period for each endpoint. None of the valuations is divisible by an odd prime. The endpoint \(x=1\) gives \(S_k(1)=1\).

## Remaining cubes and fifth powers

Only \(q=3,5\) remain at the other fifteen endpoints. The [verifier](verify.py) eliminates each case over a complete period of \(k\), allowing every positive exponent in the represented residue classes.

For a test prime \(p\) with \(q\mid p-1\mid2520\), set \(r=k\bmod(p-1)\) and write \(x=hp+s\), \(0\le s<p\). The cyclic multiplicative group of the field gives
\[
S_k(x)\equiv -h\,\mathbf 1_{r=0}+\sum_{j=1}^{s}j^r \pmod p .
\]
In each complete block the multiple of \(p\) contributes zero even when \(r=0\), since the original exponent is positive. Compare this residue with all \(q\)-th power residues including zero. For each endpoint and each \(q\), the verifier starts with all 2,520 classes of \(k\) and intersects the survivors for successive test primes up to 421. Every intersection is empty by the end. The finite certificate depends on no upper bound on \(k\) or \(x\) beyond the endpoint list supplied by the published theorem.

The general conjecture still requires other endpoint supports and the square exponent. The square and fourth-power exceptional families are outside the theorem's odd-prime hypothesis.
