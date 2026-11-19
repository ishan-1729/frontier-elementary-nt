# P02 (parts 3–4): Local certificates, global obstructions, and the research landscape
written: 2026-10-27T05:47:36Z
question: Your Runge proof for 1^15 + ... + n^15 = y^2 is correct: I re-checked every identity and inequality exactly.
status: complete (the requested continuation is complete; proposed uniform exclusion theorems are explicitly marked TARGET or UNPROVED)

This file completes Parts 3 and 4 of P02. Parts 1 and 2 were delivered in chat. It preserves P02's focus on Schäffer's equation as a test problem. P03 develops the broader theory requested subsequently; a useful theory should be judged by several such evaluations, not organized around this conjecture alone.

Throughout, S_k(n)=sum_(j=1)^n j^k, with positive integers n,k. The power exponent is m>=2. All proofs labelled PROVED below use elementary algebra, explicit congruences, convergent p-adic series, or finite certificates. A literature statement is identified separately when its proof is not reproduced.

# PART 3. LOCAL FORMULAS, THE k-DIRECTION, AND TWO-PARAMETER CERTIFICATES

## 3.1. Exact 2-adic valuations

PROVED. Let a=v_2(n(n+1)). Then

    v_2(S_k(n)) = a−1             if k=1 or k is even;
    v_2(S_k(n)) = 2a−2            if k>=3 is odd.

Proof. Let M be the even member of {n,n+1}, so v_2(M)=a, and put

    A_r(M)=sum_(i=1)^(M−1) i^r.

For a=1, there are M/2 odd summands, an odd number, so A_r(M) is odd for every r>=1.

For even r>=2 and a>=2, pair i with M−i:

    A_r(M)=(M/2)^r+sum_(i=1)^(M/2−1) [i^r+(M−i)^r].

The paired sum is congruent to 2*i^r modulo 2^a. Also r(a−1)>=a, so the central term vanishes modulo 2^a. Therefore

    A_r(M) ≡ 2*A_r(M/2) mod 2^a.

Induction on a gives v_2(A_r(M))=a−1.

For odd r>=3 and a>=2, expansion of the same paired sum gives

    A_r(M)=(M/2)^r+r*M*A_(r−1)(M/2)+M^2*E

for an integer E. The middle term has valuation

    a+(a−2)=2a−2.

The central term has valuation r(a−1)>2a−2, and M^2*E has valuation at least 2a. Hence v_2(A_r(M))=2a−2.

For r=1, the identity A_1(M)=M(M−1)/2 gives valuation a−1 directly.

Finally, S_k(n) is either A_k(M) or A_k(M)+M^k. In every case M^k has larger valuation than A_k(M), including a=1. This proves the formula.

Consequences:

- For k=1 or even k, the valuation obstruction to a square occurs when a is even. For example, n≡3 or 4 mod8 is excluded for all these k.
- For odd k>=3, the 2-adic valuation is always even. Its parity alone never excludes a square.
- If n≡0 or 3 mod4, then a>=2. Any equality S_k(n)=y^m forces m|a−1 for k=1 or even k, and m|2a−2 for odd k>=3.
- If n≡1 or 2 mod4, then a=1, and these divisibility conditions say nothing about m.

The valuation formulas are established in the literature; see [Hajdu's paper](https://math.unideb.hu/sites/default/files/inline-files/fixx23_revised2%20%281%29.pdf). Its [corrigendum](https://math.unideb.hu/sites/default/files/inline-files/skxfixx_corr2%20%281%29.pdf) repairs the proof of its 3-adic statement and explicitly leaves the stated results unchanged. The arguments here are included in full so that no such result is a black box.

## 3.2. The normalized odd part gives a stronger 2-adic obstruction

PROVED. Suppose n≡0 or 3 mod4, write the even member of {n,n+1} as

    M=2^a*d,  a>=2,  d odd,

and, for even r>=4, define

    c_r=1 if r≡0 mod4;
    c_r=5 if r≡2 mod4.

Then

    S_k(n)/2^(a−1) ≡ d*c_k mod8             for even k>=4,

and

    S_k(n)/2^(2a−2) ≡ k*(2k−1) mod8        for odd k>=5.

Proof of the even formula. First consider M=4d. Expansion into d blocks of length 4 gives, modulo 16,

    A_r(4d) ≡ d*A_r(4).

Indeed, the terms of order at least two in the block offset 4t are divisible by 16. The linear term contains 4*r*A_(r−1)(4), whose valuation is at least 2+1+2=5. At the base,

    A_r(4)=1+2^r+3^r.

Since r>=4, reduction modulo 16 gives 2 if r≡0 mod4 and 10 if r≡2 mod4. Thus

    A_r(4d)/2 ≡ d*c_r mod8.

Now let v_2(M)=a>=2. Expansion of the block from M to 2M−1 yields the exact identity

    A_r(2M)=2*A_r(M)
      +sum_(j=1)^(r−1) binom(r,j)*M^j*A_(r−j)(M)
      +M^(r+1).

The j=1 term has valuation at least

    1+a+(2a−2)=3a−1,

because r is even and r−1>=3 is odd. For j>=2 with j<r, the valuation is at least

    ja+(a−1)>=3a−1.

The endpoint term also has valuation at least 3a−1. Since 3a−1>=a+3 for a>=2,

    A_r(2M) ≡ 2*A_r(M) mod 2^(a+3).

Dividing by 2^a proves that the normalized residue is unchanged when M is doubled. Starting from 4d proves

    A_r(2^a*d)/2^(a−1) ≡ d*c_r mod8.

The optional added endpoint M^r does not change this congruence, because ra>=a+2.

For the odd formula, use the reflection i -> M−i in the full sum. For odd k,

    2*A_k(M)=k*M*A_(k−1)(M)
      +sum_(j=2)^k (−1)^(k−j)*binom(k,j)*M^j*A_(k−j)(M),

where A_0(M)=M−1.

For k>=5, the j=2 term has valuation at least

    2a+(2a−2)=4a−2.

For 3<=j<=k−1, the valuation is at least

    ja+(a−1)>=4a−1.

The last term has valuation ka. Every one of these bounds is at least 2a+2 when a>=2 and k>=5. Dividing the identity by 2^(2a−1) therefore gives

    A_k(M)/2^(2a−2)
      ≡ k*d*[A_(k−1)(M)/2^(a−1)]
      ≡ k*d^2*c_(k−1)
      ≡ k*c_(k−1) mod8.

Since c_(k−1)≡2k−1 mod8, this is the claimed congruence. The endpoint term M^k again has enough divisibility not to change it.

For odd k≡1,3,5,7 mod8, the normalized odd part is respectively

    1, 7, 5, 3 mod8.

An odd square is 1 modulo 8. Therefore:

PROVED COROLLARY. If k>=5 is odd, n≡0 or 3 mod4, and S_k(n) is a square, then k≡1 mod8.

In particular, for every even u>=4,

    S_(u^2−1)(n) is not a square when n≡0 or 3 mod4,

because u^2−1≡3 or 7 mod8.

This excludes half the endpoints uniformly for the family of Part 2, without Runge. The normalized congruence also appears as Lemma 2.3 of [Pintér's 2026 preprint](https://arxiv.org/html/2610.33567v1).

## 3.3. The remaining endpoint classes are actually soluble over Z_2

PROVED. Let k>=3 be odd and n≡1 or 2 mod4. Then

    S_k(n) ≡ 1 mod8,

and S_k(n) is a square in Z_2.

Proof. Even summands vanish modulo 8, and j^k≡j mod8 for odd j. If h=ceil(n/2), the sum of the first h odd integers is h^2. In these endpoint classes h is odd, so S_k(n)≡h^2≡1 mod8.

For completeness, every integer U≡1 mod8 has a square root in Z_2. Starting from x_3=1, suppose an odd x_j satisfies x_j^2≡U mod2^j, j>=3. The two choices

    x_j and x_j+2^(j−1)

have squares differing by 2^j modulo 2^(j+1). One choice solves the next congruence. The successive choices form a Cauchy sequence because their differences are divisible by 2^(j−1), and its limit is a square root of U.

Thus the obstruction is stronger than the absence of a valuation reduction: on these endpoint classes there is no 2-adic square obstruction at any depth. An argument that removes them must use other arithmetic information.

## 3.4. An odd-prime endpoint formula, including the leading unit

PROVED. Let p be an odd prime, let k>=3, and suppose p−1 divides k−1. Let

    a=v_p(n(n+1))>=1,
    t=v_p(k),
    k0=k/p^t.

Let M be the member of {n,n+1} divisible by p, and write M=p^a*d, with p not dividing d. Then, with s=t+2a−1,

    v_p(S_k(n))=s,

and

    S_k(n)/p^s ≡ −k0*d^2/2 modp.

Proof. We first need the complete-block congruence, valid for every integer r>=1:

    A_r(p^a*d)
      ≡ d*p^(a−1)*sum_(b=1)^(p−1) b^r mod p^a.

At a=1 this is the decomposition into blocks modulo p. For the induction, expand (i+t*p^a)^r over 0<=i<p^a and 0<=t<p. The constant part gives p*A_r(p^a). The linear part contains sum_(t=0)^(p−1) t, which is divisible by p because p is odd. Every higher part is divisible by p^(2a), hence by p^(a+1). Thus

    A_r(p^(a+1)) ≡ p*A_r(p^a) mod p^(a+1).

The factor d follows by another block decomposition modulo p^a.

In particular, v_p(A_r(M))>=a−1. If p−1 divides r, the nonzero residues have r-th power 1 modulo p, so

    A_r(M)/p^(a−1) ≡ −d modp.

The hypothesis on k makes k odd. Apply the reflection identity of the preceding section:

    2*A_k(M)=k*M*A_(k−1)(M)+higher terms.

The first term has valuation t+2a−1 and normalized unit −k0*d^2.

For 2<=j<=k−1,

    v_p(binom(k,j))>=t−v_p(j),

because binom(k,j)=(k/j)*binom(k−1,j−1). The j-th higher term therefore has valuation at least

    t+(j+1)a−1−v_p(j).

For an odd prime and j>=2, v_p(j)<=j−2. Hence

    (j−1)a−v_p(j)>=1,

and this lower bound is at least t+2a.

The last term has valuation ka. Since t<=k−2 and a>=1,

    ka>=t+2a.

Thus the first term has strictly smallest valuation. Division by 2, a p-adic unit, proves both the valuation and leading-unit assertions. The optional addition of M^k has strictly larger valuation and changes neither conclusion.

For squares, the necessary conditions at such a prime dividing n(n+1) are

    t is odd,
    −k0/2 is a square modulo p.

This uses only a finite square-residue test, not quadratic reciprocity.

For p=3, the exponent hypothesis holds for every odd k. Therefore, if 3 divides n(n+1), a square requires

    v_3(k) odd,
    k/3^v_3(k) ≡ 1 mod3.

For example, k=15 has v_3(k)=1 but 15/3≡2 mod3, so a square S_15(n) requires n≡1 mod3.

For the family k=u^2−1, if 3 divides u then v_3(k)=0. Combining this with the 2-adic obstruction proves:

    If 6 divides u, u>=6, and S_(u^2−1)(n) is a square,
    then n≡1 or 10 mod12.

For arbitrary power exponents, if p does not divide k and a=1, the valuation is 1, excluding every m>=2.

An additional simple formula is useful for even k: if p is odd, p−1 divides k, and a=v_p(n(n+1))>=1, the block congruence gives v_p(S_k(n))=a−1. The endpoint term has larger valuation. Thus a=2 excludes all perfect powers in that situation as well.

## 3.5. Bernoulli and Kummer obstructions extend the list, but do not complete it

KNOWN ELEMENTARY-CONGRUENCE RESULT, not reproved here. Let p>=7 be prime, r be even with 4<=r<=p−3, and e=v_p(B_r)>=1. If

    k≡r+1 mod p^e(p−1),
    a=v_p(n(n+1))>=1,
    2a>e,

then

    v_p(S_k(n))=2a+e.

This is Theorem 1.1 of [Pintér's 2026 preprint](https://arxiv.org/html/2610.33567v1). Its proof uses explicit Kummer congruences and the first nonzero term of the Bernoulli polynomial at an endpoint.

If e is odd, this excludes squares on the stated two-parameter progressions. These progressions can meet n≡1 or 2 mod4, so odd primes do address parts of the 2-adically soluble region.

The limitation is coverage: the theorem does not say that every remaining pair (n,k) belongs to one of these progressions. A collection of strong necessary conditions is not a complete exclusion theorem without a proof that all remaining cases are covered.

## 3.6. No fixed finite set of primes can settle the unbounded endpoint problem

PROVED. Let P be any finite set of primes. There are arbitrarily large positive integers n such that, simultaneously for every odd k>=3, S_k(n) is a square in Z_p for every p in P.

Proof. Take

    n≡1 mod [4*product_(p in P, p odd) p].

At p=2, Section 3.3 applies. At an odd p in P, complete residue blocks have zero sum of odd powers modulo p, by pairing j and −j. Since n≡1 modp,

    S_k(n)≡1 modp.

The square root 1 modulo p lifts to Z_p, since the derivative of Y^2−S_k(n) at 1 is 2, a unit. There are arbitrarily large n in the displayed progression.

This proves a limitation even when unlimited depth is allowed at each prime in P. Increasing the depth at a fixed finite list of primes cannot remove these endpoints.

The statement does not rule out primes chosen in dependence on the remaining variables, a global height bound, or a descent. Nor does it rule out a finite-prime certificate after an explicit bound has made the search region finite.

There is a general perfect-power version.

PROVED. For every fixed m>=2 and finite set P, there are arbitrarily large n for which S_k(n) is an m-th power in Z_p for all p in P and all k>=1.

For q=p^c, one has the universal endpoint period

    S_k(n+q^2)≡S_k(n) modq.

Indeed, the additional interval contains q complete blocks of residues modulo q. Put c_p=2*v_p(m)+1 and take

    n≡1 mod product_(p in P) p^(2c_p).

Then S_k(n)≡1 modp^c_p. At Y=1, the polynomial Y^m−S_k(n) has value of valuation at least c_p>2*v_p(m), so the elementary Hensel criterion gives an m-th root in Z_p.

## 3.7. The full k-direction is not p-adically continuous when p<=n

PROVED. Fix n>=p and a positive integer k0. On positive integers in a fixed class modulo p−1, k -> S_k(n) has no continuous extension to Z_p.

Take

    k_h=k0+(p−1)*p^h.

Then k_h tends to k0 p-adically and tends to infinity in the ordinary sense. For p not dividing j,

    j^k_h -> j^k0 in Q_p.

For p dividing j,

    j^k_h -> 0.

Consequently,

    lim_h S_(k_h)(n)
      =S_(k0)(n)−p^k0*S_(k0)(floor(n/p)),

which is different from S_(k0)(n). The subtracted integer is positive and nonzero, hence is nonzero in Q_p.

For p=2, the same sequence eventually preserves parity and the convergence for odd bases follows from their even powers tending to 1; the conclusion remains valid.

The analytic object is the p-depleted sum

    U_k(n)=sum_(1<=j<=n, p not dividing j) j^k.

For odd p, U_k(n) is analytic on classes k=k0+(p−1)t, t in Z_p, because each base j^(p−1) lies in 1+pZ_p. For p=2, use k=k0+2t, so j^2 lies in 1+8Z_2.

The exact decomposition is

    S_k(n)=sum_(a=0)^(floor(log_p n)) p^(ak)*U_k(floor(n/p^a)).

The U terms interpolate; the factors p^(ak) generally do not. Kummer or Morita interpolation does not remove this discontinuity of the undepleted finite sum.

There is still a useful precision-dependent statement: modulo p^N, if k>=N then all terms with p dividing j vanish, so

    S_k(n)≡U_k(n) mod p^N.

This is sufficient for finite periodic certificates. It is not a single analytic interpolation theorem for the full sum at all positive integers.

## 3.8. Even genuine interpolation does not isolate perfect-power exponents

Suppose p>n, so all bases are units and the exponent interpolation is genuinely analytic. The equation is still

    F(t,y)=S_(k0+(p−1)t)(n)−y^m=0,

one equation in the two p-adic variables t and y.

PROVED LOCAL STATEMENT. Near a solution (t0,y0) with y0!=0, the equation admits an analytic branch y=y(t) on a sufficiently small p-adic neighborhood of t0.

One can use the p-adic implicit-function theorem, since

    partial F/partial y=−m*y0^(m−1)!=0 in Q_p.

Equivalently, shrink the neighborhood until S(t)/S(t0) is sufficiently close to 1 and write

    y(t)=y0*exp[(1/m)*log(S(t)/S(t0))].

Choosing the neighborhood so that the logarithm divided by m lies in the usual convergence domain makes this an explicit convergent expression.

Therefore the existence of a local m-th root is an open condition near a nonzero m-th power. Strassmann in k does not isolate the global perfect-power exponents merely by introducing F(t,y).

A successful ghost argument must first impose additional global arithmetic on y: a finite multiplicative parametrization, extra independent markers from a covering, a descent, or a height bound. Merely counting k as a single ghost variable overlooks the freely varying y-coordinate.

The same issue persists in a jointly analytic (n,k) model. At a nonzero y, the equation can locally be a two-dimensional graph y=y(n,k).

## 3.9. A precise two-parameter finite certificate

PROVED. Let q=p^a be a prime power.

1. For every n>=0 and k>=1,

       S_k(n+q^2)≡S_k(n) modq.

2. If k,k'>=a and k≡k' modlambda(q), then

       S_k(n)≡S_(k')(n) modq.

Here lambda(q) is the exponent of the finite unit group modulo q.

The first assertion was proved in Section 3.6. For the second, a base divisible by p contributes zero to both sums; a unit contributes the same value by the definition of lambda(q).

Thus, for a fixed m and finitely many prime powers q_i, the local perfect-power conditions are a finite table in

    n modulo lcm_i(q_i^2),
    k modulo lcm_i(lambda(q_i)),

after separately handling k<max_i v_(p_i)(q_i).

Each empty cell proves uniform nonexistence on the corresponding arithmetic progression of pairs. Each surviving cell, however, represents infinitely many integer pairs.

A finite table is not a finite residual set. To make the residual set finite, one needs an additional theorem giving one of:

- an explicit bound on both unbounded parameters;
- a descent that decreases a specified nonnegative integer height;
- termination of an adaptive refinement procedure;
- a covering whose global arithmetic supplies further equalities.

None of these termination assertions follows merely from the first- or second-order local isolation engine.

## 3.10. Why shallow exponent sieves retain the cubic-square branch

PROVED. Fix n>=1. Let M>=2 be a modulus such that v_p(M)<=3 for every prime p<=n. Then every k>=3 with

    k≡3 modlambda(M)

satisfies

    S_k(n)≡S_3(n)=[n(n+1)/2]^2 modM.

Proof. For each prime power p^a dividing M, a unit base has the same k-th and third powers by the exponent congruence. If p divides a base j<=n, then p<=n and a<=3, so both powers vanish modulo p^a. Sum and combine by the Chinese remainder theorem.

Thus direct square-residue sieves using only prime moduli, or only such shallow prime powers, leave infinitely many k for every fixed n. Deeper information at a prime dividing a summand is needed to break this inherited cubic-square congruence.

There is no contradiction with the normalized formula of Section 3.2: testing its odd part modulo 8 is a test on S itself modulo 2^(2a+1), which can be much deeper than 3.

## 3.11. What the local theory already gives the even-u evaluation

PROVED, by combining Parts 2 and 3. If u>=4 is even and S_(u^2−1)(n) is a square, then all the following necessary conditions hold:

    n+1/2 < (4u^2)^(4u^2);
    n≡1 or 2 mod4.

For every odd prime p satisfying p−1 | u^2−2 and p | n(n+1), put t=v_p(u^2−1) and k0=(u^2−1)/p^t. Then also

    t is odd;
    −k0/2 is a square modulo p.

In particular, if 6 divides u, then n≡1 or 10 mod12.

These statements are uniform, elementary restrictions. They do not prove that n<=1. The remaining obstacle is global and persists in a region that is genuinely soluble at 2 and, for any fixed finite list of primes, contains arbitrarily large locally soluble endpoints before a height bound is imposed.

# PART 4. THE LANDSCAPE AND THE MOST USEFUL ELEMENTARY THEOREM TARGETS

## 4.1. First divide by the source of arithmetic restrictions

BELIEF. The useful division is not simply “small k versus large k.” It is by the mechanisms available after the obvious factors have been removed.

### Region A: endpoint valuations give a finite list of power exponents

For n≡0 or 3 mod4, the exact 2-adic valuation constrains m. On an endpoint progression where a=v_2(n(n+1)) is fixed, it gives a fixed finite set of possible m. One can then build finite two-parameter congruence certificates.

The particularly favorable case is n≡3 or 4 mod8: a=2. For even k or k=1, valuation 1 excludes every perfect power. For odd k>=3, valuation 2 forces m=2. This is an effective way of reducing an all-m problem to one square-residue table.

The prospect here is a complete theorem on an infinite progression, not merely a bound on finitely many sampled endpoints.

### Region B: arithmetic Runge holds

Part 1 proved the exact direct criterion, apart from the reducible k=3 cases:

    k+1=u^q, with q prime, q|u, and q|m.

For squares this is k=u^2−1 with u even. The methods of Part 2 produce an explicit bound for every member and an exact finite certificate once the remaining range is excluded.

The degree reduction through T=n(n+1), the best actual denominator, and sign bounds for the exact residual polynomial should be retained. Coarse uniform coefficient bounds are useful for existence of an algorithm; exact identities are usually much better for a particular equation.

### Region C: coprime factorization produces several simultaneous powers

For even k, the factors n, n+1, and 2n+1 make this region especially promising. The resulting simultaneous conditions should be retained together, as explained below.

For odd k and m>2, the squared endpoint factors still give useful higher-power conditions. For odd k and m=2 they do not: the exponent restriction on each endpoint becomes vacuous. This distinction is a major reason squares are a particularly difficult evaluation.

### Region D: a low-genus quotient has a certified small orbit

The evident reflection symmetries give elliptic quotients in the low-degree cases k=7,9,11. When an explicit elementary descent determines the quotient's relevant points, this can be decisive. Your k=9 proof is an example of this role.

The claim should be restricted to these evident quotient constructions. It would be unjustified to infer, only from degrees, that no larger k ever has an accidental elliptic map or another useful special factorization.

### Region E: n is small compared with k

Last-term dominance and saddle-point estimates provide accurate real approximations. Part 2 shows why relative estimates alone do not exclude squares. The needed additional statement is an arithmetic separation, a special denominator cancellation, or a global multiplicative parametrization.

### Region F: locally soluble points on a curve with one infinity place

This includes every square equation with even k on the original model, and prime power exponents l with l not dividing k+1. For odd k, it also includes many square equations outside the even-square family.

Direct Runge does not apply. If the obvious endpoint conditions are also locally soluble and the quotient has no small elementary orbit description, this is a genuinely difficult region. A whole-conjecture method would need a new covering, a global descent, or a new effective arithmetic separation theorem there.

## 4.2. The coprime splitting for even k, with all bad primes explicit

PROVED FINITE-MULTIPLIER STATEMENT.

Fix an even k>=2. Write

    C*S_k(n)=n(n+1)(2n+1)*A(n),

where C is a positive integer and A is an integer polynomial. Assume A(0), A(−1), and A(−1/2) are nonzero. These nonvanishing conditions can be checked from the explicit polynomial; they express simplicity of the three extracted roots.

Let d=degree A, and let S be the primes dividing the nonzero integer

    2*C*A(0)*A(−1)*[2^d*A(−1/2)].

For every solution S_k(n)=y^m, there are positive integers a,b,c and positive m-th-power-free integers alpha,beta,gamma supported on S such that

    n=alpha*a^m,
    n+1=beta*b^m,
    2n+1=gamma*c^m.

Proof. The three integer factors n,n+1,2n+1 are pairwise coprime. If a prime outside S divides n, it does not divide C or A(n), because A(n)≡A(0) modulo that prime. Therefore its valuation in S_k(n) is exactly its valuation in n and is divisible by m.

The same argument uses A(−1) for n+1. For 2n+1, the prime is odd, and reducing A(n) at n=−1/2 uses the nonzero integer 2^d*A(−1/2). Again the endpoint valuation is divisible by m.

Extracting m-th powers gives the assertion. Each multiplier has exponents in {0,...,m−1} at primes of S, so the set is finite and explicit for fixed k,m.

Every original solution therefore satisfies

    beta*b^m−alpha*a^m=1,
    gamma*c^m−2*alpha*a^m=1,
    2*beta*b^m−gamma*c^m=1.

The third displayed equation is a linear combination of the first two. There are two independent relations, not three independent markers. Nevertheless, all three paired equations can expose different local, sign, and divisibility information.

The finite-multiplier result uses only polynomial evaluation and prime valuations. It does not require unique factorization in a number field, class groups, or a unit theorem.

The hard part is uniformity. The support S depends on k through coefficients and special values involving Bernoulli numbers. Nothing above gives one fixed finite S for every k. If m is also unbounded, the multiplier list varies with m as well.

The recent [Chen–Choi–Houle preprint](https://arxiv.org/html/2610.28768v1) develops a coupled-equation descent for even k and then uses modular and logarithmic tools. The elementary splitting is available to this programme; replacing those later tools uniformly is an additional research problem.

## 4.3. The corresponding odd-k splitting explains a structural difficulty of squares

PROVED. Fix odd k>=3 and an identity

    C*S_k(n)=n^2(n+1)^2*A(n)

with C positive, A integral, and A(0)*A(−1)!=0. Outside the primes dividing C*A(0)*A(−1), every solution S_k(n)=y^m satisfies

    m | 2*v_p(n),
    m | 2*v_p(n+1).

Putting M=m/gcd(m,2), this yields

    n=alpha*a^M,
    n+1=beta*b^M,

with alpha,beta in a finite explicit set of M-th-power-free multipliers.

The proof is the same endpoint valuation argument.

For m=2, M=1, so this produces no perfect-power restriction on n or n+1.

For m=4, M=2, so one obtains coupled square or Pell-type conditions.

For even m>4, M>=3, and for odd m, M=m. The endpoint factors then carry substantially more arithmetic information.

This is a structural explanation for the usefulness of higher even exponents. It also reinforces the correction in Part 1: Pintér's 2007 theorem for odd k<170 concerns even m>4, and its stated range does not include the square cases.

## 4.4. What the obvious hyperelliptic quotients actually say

For odd k>=5, write

    S_k(n)=T^2*R_k(T),
    T=n(n+1),
    degree R_k=(k−3)/2.

After dividing a square equation by T^2 and putting X=n+1/2, the curve becomes

    w^2=R_k(X^2−1/4).

The two evident involutions give quotient models

    w^2=R_k(T)

and

    V^2=(T+1/4)*R_k(T),  V=X*w.

When the polynomials on the right are squarefree, their genera are respectively

    floor((k−5)/4),
    floor((k−3)/4).

The low-degree elliptic instances occur at k=7,9,11: for k=7 the second quotient is elliptic; for k=9 both are; for k=11 the first is. If repeated factors occur, they must be removed before using these genus formulas.

These constructions are elementary identities. A proof that the relevant quotient points are exhausted still needs its own descent or coverage certificate. The existence of a map to an elliptic curve does not itself provide a usable one-index parametrization.

## 4.5. The hardest region for the present evaluation

BELIEF, supported by the proved limitations above.

The most resistant square region has all of the following features:

- odd k, with n≡1 or 2 mod4, so the value is already a square in Z_2;
- no useful odd-prime endpoint obstruction at the selected primes;
- endpoint factors whose exponents are already even, so coprimality does not force new squares;
- no low-rank quotient with an elementary coverage proof;
- either no arithmetic Runge condition, or a Runge bound leaving an unbounded family of parameter-dependent finite intervals.

For higher powers, analogous difficulty remains when the covering multipliers grow with k and the exponent m has no independent bound.

A new mechanism in this region would have to add genuinely global information. Possible forms are:

1. A covering that imposes extra independent arithmetic conditions and still admits elementary coverage.
2. A height-decreasing transformation valid across a parameter family.
3. An integer-separation theorem for the relevant approximants, including their denominators and conjugates.
4. An adaptive-prime argument with a proved termination theorem and an explicit size bound.

Adding more local Taylor coefficients without changing the dimension of the smooth local solution set cannot supply that information by itself.

## 4.6. Ranking the attack lines

The following ranking is a BELIEF about expected mathematical return multiplied by the chance of an elementary complete result within weeks. It is not a claim about the intrinsic importance of the theories involved.

1. **Uniform endpoint-progressions with finite congruence certificates.** High feasibility, substantial power on complete infinite families of parameter triples. The proof engine is fully elementary and the required output is independently verifiable.

2. **Runge with exact denominators and residual identities, followed by a rigorous local cover.** High feasibility for individual favorable equations and selected covers; moderate feasibility for structurally uniform families. Part 2 already proves a uniform global bound, but a full family requires a further coverage theorem.

3. **Coupled perfect-power descent retaining all endpoint factors.** High potential power and a natural source of several arithmetic restrictions. A finite multiplier enumeration is elementary. A complete uniform solution of the resulting systems is the difficult missing step.

4. **Special low-genus quotients and explicit endomorphism or isogeny descent.** Very effective where present, but structurally sparse in this particular family. This remains an excellent testbed for improvements to the general ghost theory.

5. **A k-direction Strassmann argument without further parametrization.** Low feasibility as stated, for a proved reason: the full function may not interpolate, and even when it does, y varies along an analytic branch. It becomes promising only after a new global parametrization creates enough independent conditions.

6. **Last-term or saddle-point estimates by themselves.** Valuable for rigorous bounds and for constructing approximants, but presently insufficient for a complete arithmetic exclusion. Their useful combination with norm or lattice information belongs in lines 2 and 3.

The first two are the most immediate complete-proof projects. The third is the stronger candidate for a conceptual advance that transfers beyond power sums.

## 4.7. A precise top theorem: uniform finite certificates on an endpoint progression

PROVED CERTIFICATE THEOREM.

Let a>=2 be fixed. Let b,N be positive integers such that every positive n≡b modN satisfies

    v_2(n(n+1))=a.

For example, one can take N divisible by 2^(a+1) and b≡2^a or 2^a−1 modulo 2^(a+1).

Let

    E={m>=2 : m divides 2(a−1)}.

Choose finitely many prime powers q_i=p_i^e_i, and put

    K=max(6,e_1,...,e_s),
    U=lcm(N,q_1^2,...,q_s^2),
    V=lcm(2,lambda(q_1),...,lambda(q_s)).

For each integer residue r with 0<=r<U and r≡b modN, each residue c modulo V, and each m in E, choose any k_c>=K with k_c≡c modV. Suppose that for every such triple (r,c,m), at least one q_i satisfies

    S_(k_c)(r) is not an m-th-power residue modulo q_i.

Then there are no solutions

    S_k(n)=y^m,
    n≡b modN,
    n>=1, k>=K, m>=2.

If the finitely many exponent equations 1<=k<K on that progression are separately classified by elementary arguments, the same data give a complete theorem for all k,m on the progression.

Proof. Suppose a solution exists. The valuation theorem gives m|a−1 for k=1 or even k, or m|2a−2 for odd k>=3. In either case m belongs to E.

Let r be the residue of n modulo U and c the residue of k modulo V. For each i, endpoint periodicity gives

    S_k(n)≡S_k(r) mod q_i.

Since k,k_c>=K>=e_i and k≡k_c modulo lambda(q_i), exponent periodicity gives

    S_k(r)≡S_(k_c)(r) mod q_i.

But S_k(n)=y^m is an m-th-power residue modulo every q_i, contradicting the hypothesis for the triple (r,c,m).

The value of the table is independent of the choice of k_c by the same exponent periodicity. Every check is finite integer arithmetic. This proves the theorem.

The theorem is useful even when the full U-by-V table is too large to print. A certificate can be a finite tree of congruence refinements. Each leaf states an incompatible power-residue test; coverage is verified by exact residue partitions or Chinese-remainder calculations.

HARDEST NEW LEMMA FOR A CHOSEN FAMILY. Produce a covering of all relevant residue cells by such incompatible tests, and prove that any remaining low-k equations have only the listed solutions. The valuation and periodicity engine does not guarantee that this cover exists.

For a family of endpoint progressions indexed by another unbounded parameter, the new assertion must be uniform coverage. Repeating finite searches for larger samples of the parameter is evidence, not that lemma.

## 4.8. The second theorem project and its precise remaining lemma

The proved second theorem is the uniform Runge bound from Part 2:

    If u>=4 is even and S_(u^2−1)(n) is a square,
    then n+1/2<(4u^2)^(4u^2).

Together with Section 3.11 it gives an entirely elementary finite reduction for every specified u.

TARGET. Find a uniform family of exact residual inequalities, congruence certificates, or descent transformations that excludes all n>=2 satisfying these bounds and local restrictions, for every even u>=4.

The missing lemma is a global coverage or arithmetic-separation statement with u still unbounded. It is not the nonsquareness of the power-sum polynomial, which Part 1 already proved. It is not existence of a Runge cutoff, which Part 2 already proved. It is not a deeper 2-adic square test on n≡1 or 2 mod4, which Section 3.3 proves cannot work alone.

An attainable intermediate project is a complete elementary certificate for the next few chosen u, followed by an attempt to identify a single symbolic inequality or congruence pattern that explains the certificates. The literature status of each chosen square case should be checked separately before making a novelty claim.

## 4.9. Concrete computations for Sol

For the local-certificate direction:

1. Work first on endpoint progressions with a=v_2(n(n+1))=2, where m=2 is the only surviving power exponent for odd k.
2. Generate exact tables for selected prime powers, including depths greater than 3 at some primes dividing summands. Do not restrict to prime moduli: Section 3.10 explains the unavoidable surviving cubic branch.
3. Combine the tables as a coverage tree, retaining exact residue classes in both n and k. Output every leaf's modulus, exponent class, endpoint class, and power-residue contradiction.
4. Verify the tree with an independent small checker using only modular exponentiation, residue enumeration, and gcd/Chinese-remainder arithmetic.
5. Report surviving cells explicitly; do not report an empty sampled rectangle as an empty infinite residue cell.

For the Runge direction:

1. For u=4,6,8,10,12,14 and selected larger even u, generate the exact centered polynomial part q_D, its actual denominator at half-integers, the T-normalized polynomial part, and the exact residual polynomial.
2. Record Q_D(0), the first omitted centered coefficient, and the best rigorously certified cutoff from a sign or coefficient argument.
3. Compute the local admissible n classes after the normalized 2-adic test and the odd-prime valuation-and-unit tests. Use exact unit tests, not only valuation parity.
4. Search for symbolic patterns in the residual polynomial or its finite differences. Any proposed uniform sign pattern should be checked over an interval of D and then proved as an identity or coefficient inequality.
5. Keep the actual denominator in every approximation comparison. A real error less than 1 is not enough when the comparison values lie on a finer rational lattice.

For the coupled-descent direction:

1. Compute the bad-prime support from the explicit special values or resultants.
2. Retain all three endpoint factors for even k and both independent coupled equations.
3. Separate multiplier branches by local valuations and signs before introducing a unit basis.
4. On a promising branch, state the remaining global coverage lemma explicitly: unit generation, subgroup index, or height descent. These are distinct proof obligations.

## 4.10. Final assessment of P02

PROVED: the direct Runge region is exact; the even-u square family has a uniform elementary bound; normalized 2-adic congruences remove n≡0 or 3 mod4; endpoint formulas provide further odd-prime valuation and unit restrictions; finite periodic certificates can prove complete arithmetic-progression results.

PROVED LIMITATIONS: the other two endpoint classes are actual 2-adic squares for odd k; any fixed finite prime set leaves infinitely many endpoints locally soluble; full exponent interpolation fails when p<=n; and a locally varying y prevents a one-variable Strassmann argument from isolating perfect-power exponents by itself.

UNPROVED: a complete elementary proof of the even-u square family, or of the whole conjecture, is not supplied. The strongest transferable lesson is to separate coverage, local isolation, and global termination. Each can be explicit and elementary, but none follows automatically from the other two.

## Sources and provenance

- [L. Hajdu, On a conjecture of Schäffer concerning the equation 1^k+...+x^k=y^n](https://math.unideb.hu/sites/default/files/inline-files/fixx23_revised2%20%281%29.pdf).
- [Hajdu's corrigendum](https://math.unideb.hu/sites/default/files/inline-files/skxfixx_corr2%20%281%29.pdf).
- [Á. Pintér, Local obstructions and certified computations for perfect powers in power sums, arXiv:2610.33567v1](https://arxiv.org/html/2610.33567v1).
- [I. Chen, S. Choi and J. Houle, Schäffer's conjecture for even k<=42, arXiv:2610.28768v1](https://arxiv.org/html/2610.28768v1). As noted in Part 1, the v1 main theorem statement must retain the known exceptional solution.

The displayed valuation and leading-unit formulas were also checked by exact arithmetic for n<=100, k<=40, and p in {3,5,7,11,13}. Those computations are a consistency check; the proofs above supply the unbounded statements.
