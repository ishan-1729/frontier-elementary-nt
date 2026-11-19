# Square values on Pell orbits, and fifteen fourth-power indices

Let \(T=x(x+1)/2\) and \(h=(k+1)/2\). For **each of the 102 odd \(k\) from 5 through 207**, if \(S_k(x)\) is a square and the squarefree part of \(T\) equals the squarefree part of \(h\), the only possible positive endpoint is \(x=1\). It occurs exactly when \(h\) is square. At \(k=19\), for example, the class \(T=10u^2\) contains no positive square source. Its fundamental Pell solution is \(9+\sqrt{80}\), and 32 candidate Pell indices suffice.

For these **15 source indices**, \(T=u^2\) also has only \(x=1\) giving a square:
\[
        k=7,15,17,31,35,49,63,71,97,99,127,143,161,195,199 .
\]
The [fourth-power triangular-class reduction](../40_even_targets/README.md) forces \(T=u^2\) for a fourth power throughout this range. Therefore **every positive fourth power at the fifteen displayed indices has \(x=y=1\), at all heights**. This includes the eleven indices \(31,49,63,71,97,99,127,143,161,195,199\). At [63](k63/README.md), [99](k99/README.md), [143](k143/README.md), [195](k195/README.md), [255](k255/README.md) and [323](k323/README.md), integral polynomial multipliers also settle squares at arbitrary endpoints. The [unrestricted height certificates](unrestricted.md) narrow that problem when \(h\) is twice a square.

## Trace estimate

For an arbitrary odd \(k\ge5\), let \(d=h-2\), and write
\[
      S_k(x)=T^2R(T),\qquad \deg R=d,\qquad
      \operatorname{lc}R=2^{h-1}/h .
\]
On \(T=su^2\), with \(s\) positive squarefree, take \(n=2x+1\) and \(\epsilon=n+u\sqrt{8s}>1\). Here \(n^2-8s u^2=1\), the conjugate of \(\epsilon\) is \(\epsilon^{-1}\), and
\[
       T=(\epsilon-\epsilon^{-1})^2/32.
\]
When \(s=2\), \((n-4u)(n+4u)=1\) excludes positive \(u\) immediately. In every other case \(K=\mathbb Q(\sqrt{2s})\) is real quadratic. Put
\[
   \alpha^2=\frac{2^{9-4h}}h,\qquad
   g(X)=\alpha^{-2} R\bigl((X-X^{-1})^2/32\bigr).
\]
Then \(g\) is a monic Laurent polynomial of highest exponent \(2d\) with only even exponents. Take the positive \(\alpha\). The trace argument requires \(\alpha\in K\), equivalently \(\operatorname{sf}(2h)\in\{1,\operatorname{sf}(2s)\}\). This holds for \(s=\operatorname{sf}(h)\) whenever the Pell equation is nontrivial. For \(s=1\) it holds exactly when \(h\) is a square or twice a square.

Let \(\sigma=\overline\alpha/\alpha\in\{1,-1\}\). In the positive formal square root of \(g\) at infinity retain the terms \(c_jX^{d-2j}\) for \(0\le j\le\lfloor d/2\rfloor\), where \(c_0=1\). The \(c_j\) are rational. Write \(Q(X)\) for the sum with strictly positive exponents and \(c\) for the constant term, taking \(c=0\) when \(d\) is odd. Set
\[
      W(X)=Q(X)+\sigma Q(X^{-1})+c,\qquad E(X)=g(X)-W(X)^2.
\]
The highest exponent of \(E\) is at most \(d-1\). It is nonzero. For the 102 tested indices the verifier checks its full coefficient vector. In general \(R\) has a simple nonreal root for odd \(k\ge7\) by the [project root theorem](../17_all_indices/roots.md), using \((2x+1)^2=8T+1\). The map \(X\mapsto(X-X^{-1})^2/32\) has finite critical values \(0,-1/8\), so such a root lifts to simple roots, impossible in a Laurent square. For \(k=5\), the simple root \(1/4\) of \(R(T)=(4T-1)/3\) gives the same conclusion.

Let \(\|P\|_1\) denote the sum of the absolute Laurent coefficients. Write \(B\) for the denominator of \(R\), \(b\) for a common denominator of the \(c_j\), and \(\alpha=(v/w)\sqrt{\operatorname{sf}(2h)}\) in reduced rational scale. Set \(D=\operatorname{lcm}(B,wb)\). If \(S_k(x)\) is square, \(z=\sqrt{R(T)}=y/T>0\) has denominator dividing \(B\) because the square of its reduced denominator divides \(B\). Also
\[
   z-\alpha W(\epsilon)
     = z-\operatorname{Tr}_{K/\mathbb Q}\bigl(\alpha Q(\epsilon)\bigr)-\alpha c .
\]
Multiplication by \(D\) takes the first two terms on the right to rational integers: \(\epsilon,\epsilon^{-1}\) are algebraic integers and the trace of an algebraic integer is integral. The final \(D\alpha c\) is an algebraic integer. This handles the possibly irrational constant term that would be missed by treating the whole truncation as a trace.

If the displayed difference is nonzero, its absolute value is at least \(1/(D^2(1+2|c|))\). Indeed a difference of magnitude at least 1 suffices. Otherwise its multiple by \(D\) has nonzero integral norm and its conjugate has magnitude at most \(D(1+2|\alpha c|)\le D(1+2|c|)\), since \(0<\alpha\le1\). On the other hand, for \(\epsilon>2\|W\|_1\) positivity of the monic leading term gives \(W(\epsilon)>\epsilon^d/2\), and hence
\[
   |z-\alpha W(\epsilon)|
      =\frac{\alpha^2|E(\epsilon)|}{z+\alpha W(\epsilon)}
      \le \frac{2\|E\|_1}{\epsilon}.
\]
If the difference vanishes, \(E(\epsilon)=0\), and Cauchy's root bound after clearing negative powers gives \(\epsilon\le1+\|E\|_1/|\operatorname{lc}E|\). Thus every square source in the indicated class satisfies the unconditional, effectively computable bound
\[
 \epsilon\le
 \max\left\{2\|W\|_1,\;
     2\|E\|_1 D^2(1+2|c|),\;
     1+\|E\|_1/|\operatorname{lc}E|\right\}.
\]
This is an application of Runge approximation on a covering curve. The fixed-denominator rational trace uses both embeddings of the *same* Pell field. When \(\alpha\notin K\), the inclusion needed for that trace is missing. [Levin, *Variations on a theme of Runge*](https://arxiv.org/abs/0805.1345)

Berndt's proof in Part V, Chapter 34, equation (7.15), takes a unit \(\eta\) and pairs it with its inverse under an involution to descend \(\eta+\eta^{-1}\) to a smaller field. Here that operation applies term by term to \(\alpha c_j\epsilon^{d-2j}\). The Faulhaber polynomial determines the coefficients, while the source Pell equation supplies \(\bar\epsilon=\epsilon^{-1}\). The fixed denominator and the treatment of the constant term above are the additional steps needed for an arithmetic exclusion. [Part V, printed p. 263, physical p. 271](../../../notebooks/notebook-5/pages/0271.webp)

## Finite certificate

The exact [cofactor constructor](../40_even_targets/polynomials.py) uses the centered Bernoulli polynomial, with an independent polynomial interpolation check against the integer power sums. If \(R=P/B\) and \(\ell=\operatorname{lc}P\), the verifier computes \(\ell g\) using
\[
 (X-X^{-1})^{2i}=\sum_{j=0}^{2i}(-1)^j\binom{2i}{j}X^{2i-2j}.
\]
It constructs \(E\) by squaring the integral coefficient vector for \(bW\) and checks its leading cancellations and nonvanishing exactly.

For nonsquare \(8s\), write \(\epsilon=\epsilon_1^j\), where \(\epsilon_1=a+v\sqrt{8s}\) is the least positive Pell unit in these integer coordinates. The verifier finds it by the continued fraction of \(\sqrt{8s}\), stopping at the first norm-one convergent. Every solution gives such a convergent since \(|a/v-\sqrt{8s}|=1/(v^2(a/v+\sqrt{8s}))<1/(2v^2)\); denominators of the convergents increase after the initial ones. Powers exhaust the solutions: divide any positive unit by a suitable power of \(\epsilon_1\) into \([1,\epsilon_1)\); a unit there strictly above 1 has positive integral coordinates, contradicting minimality. [Pell's equation and its fundamental solution, Archive of Formal Proofs](https://isa-afp.org/entries/Pell.html)

The cutoff for \(\epsilon\) is rounded strictly upward to \(2^b\) in the table (the column `cutoff_bits` uses \(b\) only for this exponent). Since \(\epsilon_1>2a-1\), taking \(\lambda=\lfloor\log_2(2a-1)\rfloor\ge1\) gives the conservative enumeration \(1\le j\le\lfloor b/\lambda\rfloor\). For each small prime \(p\nmid B\), the Pell orbit is generated modulo \(p\) until the first return of its unit to 1. Every square source must have \(BP(sv_j^2)\) square or zero modulo \(p\), where \(v_j\) is the second Pell coordinate. Removing failures leaves only fifteen indices, each the seed \(j=x=1\). They are checked by integer square roots. Over the 108 branch rows the verifier considered 20,275 indices before these modular tests. The largest bound on \(j\) was 1,294. The exponent bits and candidate counts are given separately for every branch. [Machine-readable certificate](certificate.tsv) · [Standalone verification using the cofactor constructor](verify.py)

For comparison, the previously established [common fourth-power sieve](../40_even_targets/pell.md) gives \(j>2\cdot10^{418}\) at a nontrivial fourth power. For the fifteen listed indices this already contradicts our upper bound, even without testing the at most 1,294 indices for squares.

~~~sh
python3 research/schaffer/40_even_targets/polynomials.py
python3 research/schaffer/42_pell_trace/verify.py | diff - research/schaffer/42_pell_trace/certificate.tsv
python3 research/schaffer/42_pell_trace/unrestricted.py | diff - research/schaffer/42_pell_trace/unrestricted.tsv
~~~
