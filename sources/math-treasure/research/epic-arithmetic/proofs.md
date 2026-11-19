# Four arithmetic assertions: proofs and exact repairs

R150 / T026, 2026-10-26. Review scope: [review.md](review.md). Independently reviewed: R157 approved the four mathematical treatments without correction. The four original images were personally inspected in full. Their exact identities, source paths and hashes are recorded in source-audit.json. Historical attribution and mathematical validity are separate questions.

## 1. Multiplicative order (U-252bcc265bc6d12a)

Let p be prime and a a positive integer with p not dividing a. Multiplication by a permutes the nonzero residue classes modulo p: injectivity follows from prime cancellation, and a finite injective map is bijective. Multiplying the residues gives a^(p−1)(p−1)! ≡ (p−1)! modulo p. The factorial is invertible, so a^(p−1) ≡ 1. Thus the least positive exponent d with a^d ≡ 1 exists.

For any nonnegative integer n write n=qd+r, 0≤r<d. Then a^n≡a^r; if a^n≡1, minimality forces r=0. Conversely every multiple of d gives 1. Applying this to n=p−1 proves d|(p−1). In particular the positive exponents in the image are precisely the positive multiples of d. The case p=2 is included: a is odd and d=1. Exponent zero also gives 1; negative exponents require interpreting powers in the group of invertible residues, and are not needed for the printed assertion.

No cyclic-group structure theorem or primitive-root theorem is imported. Prime cancellation follows directly from Euclid's lemma, or Bézout's identity for a residue coprime to p.

## 2. Powers plus one (U-686102bc8027347d)

The printed odd-order assertion is false at p=2: a=1, d=1 and 2 divides a^n+1 for every nonnegative integer n. More generally all odd a have d=1 and this same exception.

The corrected statement is: for an **odd prime** p not dividing the positive integer a, and its order d, the congruence a^n≡−1, for nonnegative integers n, has solutions exactly when d is even. In that event its solutions are exactly n≡d/2 modulo d.

Proof: if a^n≡−1, then d|2n and d does not divide n, because 1≠−1 modulo odd p. If d were odd, d|2n would imply d|n, a contradiction. Conversely, for even d, u=a^(d/2) has u²≡1 but u≠1 by minimality. Since (u−1)(u+1)≡0 and p is prime, u≡−1. Finally a^n≡a^(d/2) is equivalent, by cancellation and division of exponents by d, to equality of their residues modulo d. This also covers n=0, which is not a solution for odd p.

Primality is material: modulo 8, a=3 has order 2 but a^(2/2)=3 is not −1. The image does not state a composite-modulus theorem.

## 3. Fractional squares (U-549ee9d1fd4a577a)

Precisely, if the integer N is divisible by 3 but not by 9, then no rational x,y satisfy x²+y²=N. This includes N=21, and allows zero or negative x,y; fractional squares mean squares of rational numbers. Negative N is already impossible by nonnegativity.

Choose a common positive denominator and reduce the triple: x=X/Z, y=Y/Z with integers X,Y,Z, Z>0 and gcd(X,Y,Z)=1. A putative representation gives X²+Y²=NZ². Squares modulo 3 are 0 and 1, so the equation forces 3|X and 3|Y. Primitive normalization then forces 3 not dividing Z. But the left side is divisible by 9, whereas the right side is not, since 3|N, 9 does not divide N and 3 does not divide Z. This contradiction explicitly handles rational denominators; simply applying integer-square residues to rational x,y would not be valid.

A useful precisely bounded extension is that every nonzero rational sum x²+y² has even 3-adic valuation. Here v_3(r) means the exponent of 3 in the reduced numerator minus that in the denominator. For integers X,Y not both zero remove their common largest power 3^t. The remaining square sum is 1 or 2 modulo 3, since at least one entry is nonzero modulo 3. Thus v_3(X²+Y²)=2t, and division by Z² subtracts the even number 2v_3(Z). No sufficiency for being a sum of squares follows from this one necessary condition. For example 9=0²+3² shows why divisibility by 3 alone is insufficient for the quoted obstruction. Integer divisibility language in the source must not be applied without definition to arbitrary rational N.

## 4. Positive rational cube porism (U-02f5245865a45f1e)

**Theorem.** For every a,b in Q with a>b>0 there is an explicit terminating rational algorithm returning c,d in Q with c,d>0 and c³+d³=a³−b³. Neither integral outputs nor primitive outputs are asserted. The owner filename's tentative FLT application is not proved by this assertion.

Set N=a³−b³>0. At each stage, for a positive rational pair A>B with A³−B³=N, put

    C = A(A³−2B³)/(A³+B³),
    D = B(2A³−B³)/(A³+B³).                         (1)

The denominator is strictly positive; all operations are rational. Write U=A³,V=B³. The exact polynomial identity

    U(U−2V)³ + V(2U−V)³ = (U−V)(U+V)³             (2)

follows by expansion: each side is U⁴+2U³V−2UV³−V⁴. Consequently C³+D³=N. Also D>0 since U>V>0. If U>2V, then C>0 and the algorithm returns (C,D).

Equality U=2V is impossible for positive rationals: if A/B=r/s in lowest positive integer terms, r³=2s³ has 2-adic valuations 3v_2(r)=1+3v_2(s), a contradiction modulo 3. Thus the only other case is V<U<2V, where C<0. Set the next pair to

    A_new=D, B_new=−C.

Both entries are positive and A_new³−B_new³=N, so A_new>B_new. This is the same problem with its invariant N unchanged.

**Strict contraction.** In this continuing case let t=A/B>1. Then

    B_new/B = t(2−t³)/(t³+1) < 1/2.

Indeed the ratio is positive, and the difference after multiplying the desired inequality by 2(t³+1)>0 is

    t³+1−2t(2−t³)
      = 2t⁴+t³−4t+1
      = (t−1)(2t³+3t²+3t−1)>0.

The last factor exceeds 7 when t>1. After j continuing steps, B_j≤b/2^j (strict for j>0). Choose any integer K≥0 with b³/8^K<N; such a K exists by the Archimedean property, or by repeated exact division by 8. If K steps have occurred, then B_K³<N and A_K³=N+B_K³>2B_K³, so the next evaluation of (1) returns positive outputs. There are at most K continuing steps and K+1 formula evaluations. This is a proof of termination for rational arithmetic, not a numerical iteration heuristic or an elliptic-group assumption.

This proof also identifies where rationality matters: the boundary A³=2B³ cannot occur for rational A,B. Real existence alone would neither establish rational outputs nor justify this boundary exclusion.

**Construction mechanism.** On X³+Y³=N, the known point (A,−B) has tangent line A²X+B²Y=N. Parameterize it as (X,Y)=(A+B²u,−B−A²u). Substitution gives

    X³+Y³−N = −3ABN u² − N(A³+B³)u³.

The nonzero intersection parameter is u=−3AB/(A³+B³), giving exactly (C,D) in (1). Expansion (2) independently proves preservation. This elementary tangent construction is classical in mechanism; no elliptic-group theorem is needed and no novelty claim is made.

**Example.** For (a,b)=(2,1), formula (1) returns (4/3,5/3); their cube sum is 7. For (a,b)=(6/5,1), the first C is negative, and the algorithm performs the proved contraction before returning two strictly positive rational outputs. The retained verification record includes this trace. The cases a=b or b=0 are outside the printed positive strict-order domain, and are not silently handled by claiming (1) proves them.
