# All perfect powers in the fifteenth power sum

**Theorem.** If \(x,y\) are positive integers and \(m\ge2\), then
\[
1^{15}+2^{15}+\cdots+x^{15}=y^m
\quad\Longleftrightarrow\quad x=y=1.
\]
For \(S_{13}\), the same conclusion holds when \(m\) has a prime divisor outside \(\{2,5,23\}\). There is no bound on the endpoint or the target exponent in either assertion. We make no claim about historical priority.

The missing step after [the finite prime-exponent sieve](../10_prime_exponents/README.md) was an upper bound on a prime exponent \(q\) for **every** marked assignment. Squaring the two consecutive endpoint equations makes their coefficient ratio independent of \(q\). An explicit theorem on two rational logarithms then gives a small enough bound to exhaust.

## A bound from consecutive endpoints

Let \(s,t>0\) be consecutive integers in either order, \(x=\min(s,t)\), and suppose
\[
(s/t)^e=\beta^q/\alpha,\qquad e=1\ \text{or}\ 2,
\]
with positive rationals \(\alpha,\beta\). Write \(h(a/b)=\log\max(a,b)\) in lowest terms. Suppose \(\alpha\ne1,\ \beta\ne1\) and take \(H\ge\max(1,h(\alpha))\). Put \(b=h(\beta)>0\). The nonzero real linear form
\[
\Lambda=q\log\beta-\log\alpha
\]
satisfies \(|\Lambda|=e\log(1+1/x)<e/x\). Also \(qb=h(\beta^q)\le H+e\log(x+1)\) by the product inequality for rational heights.

If the two rationals are multiplicatively dependent, their prime valuation vectors span a rational line. Its integer lattice has a primitive generator, giving a rational \(g>1\) of which both are integral powers. Therefore \(h(g)\le H\) and \(\Lambda\) is a nonzero integer multiple of \(\log g\). In lowest terms \(g=a/c>1\), so \(\log g\ge(a-c)/a\ge \exp(-H)\). Consequently \(x<e\exp H\) and \(qb<(e+1)H+e\log(e+1)\).

If \(\alpha,\beta\) are on opposite sides of 1, then \(|\Lambda|\ge|\log\alpha|\ge\exp(-H)\), giving the same small bound. For independent rationals on the same side of 1, invert both if needed, which preserves heights and \(|\Lambda|\). They now both exceed 1. Apply Laurent's Corollary 1 with field parameter \(D=1\), table parameter 10, coefficient pair \(1,q\), and height parameters with logarithms \(H,\max(1,b)\). The resulting \(b'\) is at most \(q+1\). Using \(32.3<33\) gives
\[
\log|\Lambda|\ge-33H\max(1,b)\,M(q)^2,\qquad M(q)=\max(10,\log(q+1)+1).
\]
Combining this with \(\log|\Lambda|<\log(2e)-\log(x+1)\) gives
\[
q b < 33 e H\max(1,b)\,M(q)^2+H+e\log(2e).
\tag{1}
\]
All logarithms are natural. [M. Laurent, *Linear forms in two logarithms and interpolation determinants II*, Corollary 1 and Table 1, Acta Arith. 133 (2008)](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/133/4/82485/linear-forms-in-two-logarithms-and-interpolation-determinants-ii)

This calculation applies directly to any marked cover \(s^e=A U^q,\ t^e=B V^q\) by taking \(\alpha=B/A,\ \beta=U/V\). The fixed-coefficient form of the [marked factor cover](../03_marked_graphs/Schaffer_marked_graphs.md) thus supplies bounds label by label. If \(\beta=1\), the ratio \((s/t)^e=1/\alpha\) determines the coprime positive endpoints outright. If \(\alpha=1\) and \(\gcd(q,e)=1\), both endpoints must be positive \(q\)-th powers, which cannot be consecutive. The bound in (1) depends on the label height and does not supply a single absolute bound for unbounded \(k\).

## Specialization to the two odd indices

Take \(q\ge5\) prime. The [valuation calculation](../10_prime_exponents/README.md) distributes the marked odd primes of \(x(x+1)\) between the odd endpoint \(s\) and the even endpoint \(t\). Write \(L,R\) for the products of those primes at each endpoint, respectively. Use the positive endpoints throughout this height argument.

- For \(S_{15}\) every marked prime is in \(\{5,7\}\), and \(s=L^{(q-1)/2}u^q,\ t=2R^{(q-1)/2}v^q\). Thus \(e=2,\ \beta=Lu^2/(Rv^2),\ \alpha=4L/R\).
- For \(S_{13}\), write \(L=C_L D_L,\ R=C_R D_R\), where the \(C\)'s contain the assigned primes from \(\{3,5,7\}\) and the \(D\)'s contain the assigned prime from \(\{691\}\). An empty product is 1. The odd endpoint is \(C_L^{(q+1)/2}D_L^{(q-1)/2}u^q\), and the even endpoint is twice the analogous expression with subscript \(R\) and root \(v\). Squaring gives \(e=2,\ \beta=Lu^2/(Rv^2),\ \alpha=4C_RD_L/(C_LD_R)\).

If there is no marked odd prime, the endpoints have coefficients 1 and 2. The known positive solution \((1,1)\) of \(\lvert u^q-2v^q\rvert=1\) is the only one by [Bennett's uniqueness theorem](https://personal.math.ubc.ca/~bennett/B-Crelle2.pdf). It gives \(x=1\).

Otherwise \(Lu^2\) and \(Rv^2\) are coprime, as their prime divisors belong to distinct endpoints. Their maximum is at least 3. Thus \(b\ge\log3>1\). The valuation of \(\alpha\) at 2 is 2. For \(S_{15}\) its height is at most \(\log140<5\). For \(S_{13}\) it is at most \(\log(4\cdot3\cdot5\cdot7\cdot691)<13\).

Putting \(e=2\) in (1) gives \(q<66H M(q)^2+H+3\), since \(2\log4<3\). The dependent case gives the smaller bound \(q<3H+2\log3\). For \(S_{15}\), take \(H=5\) and evaluate at \(Q=100{,}000\), where \(M(Q)<13\) and \(330\cdot13^2+8<Q\). For \(S_{13}\) take \(H=13,\ Q=200{,}000\), where \(M(Q)<14\) and \(858\cdot14^2+16<Q\). The ratio \(M(q)^2/q\) decreases for \(q\ge100{,}000\). On the logarithmic branch this follows from \(2q<(q+1)(\log(q+1)+1)\), and on the constant branch it is immediate. Hence the inequalities exclude the entire tail \(q\ge Q\).

## Exhaustion and recovery

For each remaining prime, [the finite-field certificate and its two global returns](../10_prime_exponents/README.md) apply with the following executed ranges. The script enumerates assignments and all possible \(q\)-th powers at each test prime, including zero. It rejects a mixed assignment only when its two endpoints and cofactor have no simultaneous residue. It verifies the brackets for the unit-coordinate cases in the surviving binomial equation. For the other coordinates the proof uses [Bartolomé and Mihăilescu, Theorem 1.2(A)](https://arxiv.org/pdf/1412.5798). Both numerical ranges are within that theorem's bound. Its totient condition is automatic for \(S_{15}\) and removes precisely 5 and 23 among primes \(q\ge5\) for \(S_{13}\).

| Sum | Prime exponents tested | Number | Largest \(h\) used for a test prime \(hq+1\) |
|---|---|---:|---:|
| \(S_{15}\) | \(5\le q<100{,}000\) | 9,590 | 94 |
| \(S_{13}\) | \(5\le q<200{,}000,\ q\notin\{5,23\}\) | 17,980 | 144 |

For \(S_{15}\), the [square proof](../07_notebook_continuation/Notebook_Methods_and_Complete_Cases.md) traps the square root strictly between consecutive points on the \(1/36\) grid at \(x\ge16\) and checks the smaller endpoints. The [cube proof](../02_uniform_k15/Schaffer_uniform_attempt.md) resolves its residual cubic Thue equation using Delaunay–Nagell and an exact continued-fraction certificate. For \(S_{13}\) the [separate cubic calculation](../01_k13_cubic/Schaffer13_cubic_complete_proof.md) handles \(q=3\). Every \(m\)-th power is a \(q\)-th power for each prime \(q\mid m\). Finally \(S_k(1)=1\) for both indices.

Run the [verifier](verify.py) with Python 3.10 or later. It uses the standard library, checks the numerical logarithm estimates with rational lower bounds on exponentials, reruns both extended sieves, and reruns the archived \(S_{15}\) square and cube certificates. The all-index conjecture remains open in this work.
