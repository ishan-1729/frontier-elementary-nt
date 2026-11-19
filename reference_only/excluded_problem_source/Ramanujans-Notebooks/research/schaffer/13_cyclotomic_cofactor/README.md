# All odd prime target exponents for the thirteenth power sum

**Theorem.** For every odd prime \(q\), \(S_{13}(x)=y^q\) with \(x,y>0\) implies \(x=y=1\). Consequently the same holds for every target exponent with an odd prime divisor. The [relative-unit proof](../16_relative_units/README.md) supplies the remaining square case.

The [logarithm bound and first sieve](../12_endpoint_logarithms/README.md) left \(q=5,23\) among odd primes. The obstruction to using the binomial theorem at these primes is that both divide \(691-1\). Factoring the binomial equation locates the entire exceptional prime contribution in one of two factors. One allocation gives a known diagonal Nagell–Ljunggren equation. The other gives a new finite residue test in which the original Faulhaber cofactor matters.

## Endpoint reduction

For each of 5 and 23, repeat the [stage 10 consecutive-endpoint and cofactor sieve](../10_prime_exponents/README.md) on all 81 assignments. For 5 the primes 11, 31, 41, 61 leave successively 16, 8, 3, 2 assignments. For 23 the primes 47, 139, 277 leave 4, 3, 2. In each case the no-marked-prime assignment gives only \(x=1\) by Bennett. In the other assignment all four primes belong to the even endpoint. Orient as in stage 10, with \(n=u^q,\ n+1=Bv^q\) where \(n=x\) or \(-1-x\). Set \(X=-u,\ Z=-v,\ r=(q-1)/2\). Then
\[
X^q-1=BZ^q,\quad B=210\cdot72555^r,\qquad
T=\frac{X^q(X^q-1)}2,\quad R(T)=691 w^q.
\]
Here \(T=x(x+1)/2\) and \(R(T)=960T^5-2800T^4+4592T^3-4720T^2+2764T-691\). Both signed endpoint roots are nonzero. Because \(B>2\), \(X=0,\pm1\) cannot give such a solution.

## Allocating the primes after cyclotomic factorization

Write \(F=X-1\), \(G=\Phi_q(X)=1+X+\cdots+X^{q-1}\). Then \(FG=BZ^q\), \(\gcd(F,G)\mid q\), and \(v_q(G)=\epsilon\), where \(\epsilon=1\) if \(X\equiv1\bmod q\) and 0 otherwise. Indeed \(G\equiv(X-1)^{q-1}\bmod q\), and in the first case expanding \(((1+F)^q-1)/F\) gives \(G\equiv q\bmod q^2\). For \(q=5\), \(5\mid B\) forces \(X\equiv1\bmod5\), so \(\epsilon=1\).

If \(p\ne q\) divides \(G\), then the order of \(X\bmod p\) is \(q\). Among the primes of \(B\), only 691 can have that property, since the others are \(2,3,5,7\) and \(q\nmid(p-1)\) for each of them. Define \(\delta=1\) if \(691\mid G\), and 0 otherwise. The factorization now gives
\[
G=g V^q,\quad g=q^\epsilon 691^{\delta r},
\qquad F=f U^q,\quad
f=\prod_{p\mid Bq}p^{(v_p(B)-v_p(g))\bmod q},
\]
where each exponent representative is in \([0,q-1]\). The roots are integers, with signs permitted because \(q\) is odd. This construction applies equally to positive and negative \(X\).

If \(\delta=0\), one obtains \(\Phi_q(X)=q^\epsilon V^q\). The established diagonal Nagell–Ljunggren exclusion covers \(q=5,23\) for integral \(X,V\) with \(|X|,|V|>1\). For 5 the small-prime result is due to Bugeaud, Hanrot and Mignotte. For 23 use Mihăilescu's exclusion in the computed class-number range, which includes every prime from 17 through 23. His introduction states the signed equation and the small-prime input, and his conclusion gives the computed range. Here \(X\) is odd and \(|X|\ge3\), so \(G>q\) even when \(X<0\), from \(G=(|X|^q+1)/(|X|+1)\) in that case. Thus \(V>1\). [Mihăilescu, *Class Number Conditions for the Diagonal Case of the Equation of Nagell and Ljunggren*, §§1, 6](https://www.researchgate.net/publication/226863681_Class_Number_Conditions_for_the_Diagonal_Case_of_the_Equation_of_Nagell_and_Ljunggren) · [Bugeaud, Hanrot and Mignotte, Proc. London Math. Soc. 84 (2002)](https://doi.org/10.1112/plms/84.1.59)

## Restoring the Faulhaber cofactor

For \(\delta=1\), let \(H\) be the set of all \(q\)-th powers modulo the indicated test prime \(p\), including zero. Enumerate residues satisfying \(X\in1+fH\) and \(\Phi_q(X)\in gH\). For each survivor impose \(R(X^q(X^q-1)/2)\in691H\).

| \(q\) | \(\epsilon\) | \(p\) | \(X\bmod p\) after the first two tests | After the cofactor test |
|---:|---:|---:|---|---|
| 5 | 1 | 11 | none | none |
| 23 | 0 | 277 | 48, 190, 237 | none |
| 23 | 1 | 139 | none | none |

The [standard-library verifier](verify.py) enumerates the residue sets directly. Together with the published diagonal result and the preceding stages this exhausts the two missing exponents at arbitrary height. The square case is resolved separately in stage 16.
