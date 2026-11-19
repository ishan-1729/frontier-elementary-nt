# Off-anchor valuations in odd power sums

Write \(S_k(x)=\sum_{j=1}^x j^k\). In the notation of this note, the conjecture says that for \(k\ge1,\ m\ge2,\ x,y>1\) with \((k,m)\) outside \(\{(1,2),(3,2),(3,4),(5,2)\}\), the only solution to \(S_k(x)=y^m\) is \((k,m,x,y)=(2,2,24,70)\). This formulation follows [Hajdu](https://math.unideb.hu/sites/default/files/inline-files/fixx23_revised2%20%281%29.pdf). The argument below covers infinitely many endpoint classes as \(k\) varies. No priority claim is made.

## A prime for each odd exponent

**Proposition.** For every odd \(k\ge5\) one can choose a prime \(p\ge 11\) with \(p\mid 2^k+1\) and \(p\nmid k\).

Let \(q\) be the largest prime divisor of \(k\). If \(q\ge5\), then \(2^q+1=(3-1)^q+1\equiv3q\bmod9\). Thus \((2^q+1)/3>1\) is prime to 3. Choose any prime divisor \(p\) of this quotient. The order of 2 modulo \(p\) divides \(2q\) and does not divide \(q\). Order 2 would force \(p=3\), so the order is \(2q\). Thus \(p>q\), and \(k/q\) is odd, giving both assertions. Otherwise \(k=3^a\) with \(a\ge2\). Take \(p=19\), since \(19\mid 2^9+1\) and \(k/9\) is odd.

## Exact lifting

**Theorem.** For odd \(k\ge5\), let \(p\ge5\) be any prime dividing \(2^k+1\) with \(p\nmid k\). Among the \(p\) lifts modulo \(p^2\) of \(x\equiv2\bmod p\), exactly \(p-1\) have \(v_p(S_k(x))=1\). For these endpoints \(S_k(x)\) is excluded from every integer power of exponent \(m\ge2\). More generally, for every \(h\ge1\) there is precisely one class \(r_h\bmod p^h\) above 2 for which \(S_k(r_h)\equiv0\bmod p^h\). Of its lifts, \(p-1\) have valuation exactly \(h\).

Let \(A_h(x)=S_k(x+p^h)-S_k(x)\) and \(c=k(1+2^{k-1})\equiv k/2\bmod p\), a unit. The multiplicative group modulo \(p\) gives \(\sum_{i=1}^p i^{k-1}\equiv0\bmod p\), because \(p-1\mid k-1\) would imply \(2^k\equiv2\equiv-1\bmod p\). Pairing \(i\) with \(p-i\), and using oddness of \(k\), gives

\[
2\sum_{i=1}^{p-1}i^k\equiv kp\sum_{i=1}^{p-1}i^{k-1}\equiv0\pmod{p^2}.
\]

Shift the complete block \(1,\ldots,p\) to \(3,\ldots,p+2\). The two boundary differences contribute \(pk(1+2^{k-1})\bmod p^2\). Shifting this new block by any multiple of \(p\) adds zero modulo \(p^2\), by the same vanishing sum of \((k-1)\)-st powers. This proves \(A_1(x)\equiv pc\bmod p^2\) for \(x\equiv2\bmod p\).

For \(h\ge1\), divide \(A_{h+1}(x)\) into the \(p\) translates of the block defining \(A_h(x)\). Expand \((x+i+t p^h)^k\), where \(1\le i\le p^h,\ 0\le t<p\). The linear term vanishes modulo \(p^{h+2}\): \(\sum_t t\) is divisible by \(p\) and \(\sum_i(x+i)^{k-1}\) is divisible by \(p\). The quadratic term vanishes because \(p\mid\sum_t t^2\) for \(p\ge5\) and \(2h+1\ge h+2\). Terms of degree at least three contain \(p^{3h}\). It follows inductively that

\[
A_h(x)\equiv p^h c\pmod{p^{h+1}}.
\]

In particular \(S_k\) modulo \(p^h\) is constant on each class modulo \(p^h\) in this disk. Start with \(r_1=2\), since \(S_k(2)=1+2^k\). If \(S_k(r_h)\equiv0\bmod p^h\), summing the last congruence over \(d\) blocks gives the next digit by

\[
S_k(r_h+d p^h)/p^h\equiv S_k(r_h)/p^h+d c\pmod p.
\]

Precisely one digit makes this zero. Every other digit yields the specified exact valuation. If \(S_k(x)=y^m\), every nonzero valuation is a multiple of \(m\), proving the exclusion. For the first digit the sole surviving choice is \(d\equiv-2(1+2^k)/(pk)\bmod p\), interpreting the integer division by \(p\) before inverting \(k\bmod p\).

This applies the [Ghost factor-cover](../03_marked_graphs/Schaffer_marked_graphs.md) requirement on the valuations of an actual source. Here \(p\) divides neither \(x\) nor \(x+1\). It tests the remaining power-sum factor directly. The proof uses residue-class dissection and cancellation in a full multiplicative group. The [notebook concordance](../notebook-methods/README.md) records where these mechanisms were examined.

## A two-prime family for every target exponent

For \(k=35\), take \(p=11,43\). The theorem applies at both primes. Choose exact orders 2 and 3 respectively. They are coprime, so every \(m\ge2\) is excluded even though neither selected order is 1. One resulting family is

\[
x\equiv 21770954053\pmod{54605089572},\qquad
k\equiv35\pmod{2020272870},\quad k,x>0.
\]

The first congruence combines \(x\equiv13\bmod 11^3,\ x\equiv29285\bmod43^4,\ x\equiv1\bmod12\). Direct summation at these small representatives gives valuations 2 and 3 for \(k=35\), stable throughout the endpoint classes by the block congruence. The final modulus on \(k\) is the least common multiple of \(\varphi(11^3)\) and \(\varphi(43^4)\). Euler's theorem holds termwise for units, and all multiples of the respective primes vanish modulo the stated prime powers since \(k\ge35\). This proves the assertion on the entire product of progressions. The small representatives and the modular computations can be reproduced with [verify.py](verify.py).

## Remaining gap

The construction selects a prime from \(S_k(2)\), and controls endpoints in the disk \(2\bmod p\). For any fixed \(k\) the finite set of selected primes leaves, for example, every endpoint congruent to 1 modulo their product outside these disks. Covering those endpoints requires a new argument. General local valuation and exponent sieves already appear in [Pintér, §§2 and 5](https://arxiv.org/html/2610.33567v1). The [even-exponent descent of Chen, Choi and Houle](https://arxiv.org/abs/2610.28768) treats even \(k\le42\). Neither cited range establishes uniform source coverage here.
