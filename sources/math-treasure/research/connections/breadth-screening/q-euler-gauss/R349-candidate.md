# R349 candidate: the full recursive q-Euler–Gauss divisibility criterion

2026-10-27. **Complete proposed proof, awaiting independent review. No solved or priority claim.**

## Exact primary target

Narayan–Kiran, *On Some Generalisations of Gauss Sequences*,
[arXiv2511.19503v5](https://arxiv.org/pdf/2511.19503v5), September22 2026,
Definition4.1 (actualPDFp14) and Conjecture4.10 (actualPDFp19), consider a sequence a_n(q)∈Z[q], n≥1,
with, for every n≥1,

    product_(d|n, μ(d)=1) a_(n/d)(q^d)
      ≡ product_(d|n, μ(d)=−1) a_(n/d)(q^d) mod[n]_q,    (EG)
    [n]_q=1+q+...+q^(n−1).

They define A_1(q)=1 and, for n>1,

    A_n(q)=product_(d|n,d>1,μ(d)=1) a_(n/d)(q^d)
              *product_(d|n,μ(d)=−1) A_(n/d)(q^d).      (A)

Their precise conjecture is

    [n]_q divides A_n(q)*sum_(d|n) μ(d)a_(n/d)(q^d)
      in Z[q], for every n≥1.                          (C)

The source proves at most two distinct prime factors and reports finite
computational checks beyond them. The candidate below proves the full
statement, including repeated prime factors, zero polynomials and zeros
at roots of unity. It does not characterize all q-analogs of integer
Euler–Gauss sequences, give a new cyclic-sieving action, or assert the
unmultiplied q-Gauss congruence.

The archive mechanisms are [reviewed finite-character cancellation](../../../sum-products/README.md)
and the [reviewed root-of-unity/geometric return](../../cyclovariety-torus/README.md).
This is a new connection deduction from those mechanisms, not a claim that
the owner formulated Conjecture4.10. The source definitions and recursive
multiplier are Narayan–Kiran's. Finite Möbius inversion, prime-toggle pairing
and evaluation at roots of unity are classical; no novel general method or
first-priority claim is made. Original Math files were not reread here.

## 1. A finite divisor-lattice lemma with zero values allowed

Let F be any field. Fix n>1 and a proper positive divisor t of n. Give
values x_k∈F for every positive divisor k of n. For k|n write

    E_k=product_(d|k,d>1,μ(d)=1) x_(k/d),
    O_k=product_(d|k,μ(d)=−1) x_(k/d),
    G_k=sum_(d|k) μ(d)x_(k/d).

Define B_1=1 and, for k>1,

    B_k=E_k*product_(d|k,μ(d)=−1) B_(k/d).              (B)

Empty products are1. Suppose

    x_k E_k=O_k for every k|n with k not dividing t.    (H)

**Lemma.** Under these hypotheses, B_n G_n=0. More precisely, if B_n≠0,
then x_k=x_gcd(k,t) for every k|n. No nonzero assumption on the x_k is imposed.

**Proof by strong induction on n.** If B_n=0, the claimed product is zero,
so assume B_n≠0. Every prime p dividing k occurs among the divisors with
μ(p)=−1 in(B). Thus B_k has B_(k/p) as a factor. Iterating along a prime
factorization of n/k shows that B_k is a factor of B_n for every k|n.
Consequently every B_k is nonzero. In particular E_n≠0.

For any proper k|n with k not dividing t, put t_k=gcd(k,t)<k. The restricted
array indexed by divisors of k satisfies(H) with t replaced by t_k:
for j|k, j|t_k is equivalent to j|t. The induction hypothesis applies,
giving B_k G_k=0. Since B_k≠0 in the field, G_k=0.

Finite Möbius inversion now gives, for every proper k|n,

    x_k=sum_(e|k) G_e
       =sum_(e|k,e|t) G_e
       =sum_(e|gcd(k,t)) G_e
       =x_gcd(k,t).                                   (1)

Here G_e=0 was proved for every e|k not dividing t. For completeness,
the inversion identity follows by reversing the two finite divisor sums
and using sum_(d|r)μ(d)=0 for r>1,1 for r1. It is valid in every field,
including positive characteristic.

Choose a prime p with v_p(t)<v_p(n), possible because t is a proper divisor.
On squarefree divisors d of n, toggle whether p divides d. This is an
involution reversing the sign of μ(d) and preserving gcd(n/d,t): the two
possible p-exponents of n/d are v_p(n) and v_p(n)−1, both at least v_p(t),
and all other prime-exponents agree. Thus the multisets

    {gcd(n/d,t): μ(d)=1}
      and {gcd(n/d,t): μ(d)=−1}                        (2)

are equal, including multiplicities. Put y_k=x_gcd(k,t). Equality(2)
gives both

    sum_(d|n) μ(d)y_(n/d)=0,
    product_(d|n,μ(d)=1)y_(n/d)
      =product_(d|n,μ(d)=−1)y_(n/d).                   (3)

For every proper k|n, (1) identifies x_k=y_k; also y_n=x_t. The second
identity in(3) therefore says x_t E_n=O_n. Hypothesis(H) at k=n says
x_n E_n=O_n. Since E_n≠0 in a field, x_n=x_t. Combined with(1), this
proves the stronger conclusion for all k|n. The first identity in(3)
then gives G_n=0, completing the induction and the lemma. □

The only cancellations occur after B_n≠0 has forced the relevant factors
nonzero. This is why arbitrary zeros, including x_1=0, do not invalidate
the argument. We do not claim this field proof for arbitrary rings with
zero divisors.

## 2. Apply the lemma at every required root of unity

Fix n>1 and any ζ≠1 with ζ^n=1. Let m>1 be its multiplicative order,
and set t=n/m, a proper divisor of n. For every k|n define

    x_k=a_k(ζ^(n/k)).                                 (4)

If k does not divide t, then w=ζ^(n/k) is not1. Also w^k=ζ^n=1.
Hence [k]_w=0. Evaluate the assumed q-Euler–Gauss congruence(EG) at
index k and q=w. It becomes exactly

    x_k E_k=O_k,

because a_(k/d)(w^d)=a_(k/d)(ζ^(n/(k/d)))=x_(k/d).
The lemma's hypotheses hold over F=C. Notice that indices k dividing t
have w1, where [k]_1=k need not vanish; **no congruence at q1 is used
to assert an equality there**.

The recursive definition(B) agrees term for term with(A) under the same
substitutions. Induction on k|n therefore gives

    B_k=A_k(ζ^(n/k)).                                 (5)

In particular B_n=A_n(ζ), while

    G_n=sum_(d|n) μ(d)a_(n/d)(ζ^d).

The lemma proves that A_n(ζ)G_n=0 for this ζ. This holds for **every**
nontrivial nth root, not only primitive nth roots.

## 3. Return to integer-polynomial divisibility

The polynomial

    R_n(q)=A_n(q)*sum_(d|n) μ(d)a_(n/d)(q^d)

lies in Z[q] and vanishes at every nontrivial nth root of unity. These
are the n−1 distinct roots of the monic polynomial[n]_q. Divide R_n by
[n]_q using monic polynomial division in Z[q]. Its remainder has integer
coefficients and degree less than n−1, and it vanishes at n−1 distinct
complex numbers. The remainder is therefore zero. This proves(C) in
Z[q]. At n1, [1]_q=1 and the divisibility is automatic. □

The proposed proof establishes Conjecture4.10 with its full stated sequence
quantifier. It does not replace the recursive multiplier by1 or add a
nonvanishing hypothesis. It uses no unproved infinite arithmetic, analytic
or topological assumption. Independent scrutiny of this complete argument
is still required before promotion to a checked result.

## 4. Exact controls and honest status

Run

    python experiments/connections/breadth-screening/q-euler-gauss/verify.py

The exact verifier checks prime-toggle gcd multiset identities and recursive
monomial factor propagation for n2–150 and every proper divisor t. It
exhaustively enumerates scalar divisor arrays for n in{2,4,6,8,12,18,24,30}
over F₂,F₃, and n6,12 over F₅, for every proper t|n, retaining every zero
pattern satisfying(H). It verifies B_n G_n=0 and the stronger collapsed
array conclusion whenever B_n≠0. These are bounded controls for a universal
field lemma, not a substitute for its proof.

Polynomial controls include three known divisor-sum q-Gauss families for
n1–30, and the source's strictly q-Euler–Gauss Example4.11 for n2–60.
The latter has additive failures, including n12, and checks the recursive
zero multiplier without illicit cancellation. Computations use exact
integer/modular/SymPy polynomial arithmetic, not numerical roots of unity.

An initial known-q-Gauss test fixture incorrectly indexed its free divisor
coefficient as g_d instead of g_(n/d); it failed its own q-Gauss assertion.
The fixture was corrected. The universal field/zero-pattern checks had
passed before that failure. This was a verifier-control error, not an
unreported change in the candidate theorem or proof.

The retained primary v5 explicitly labels Conjecture4.10 conjectural; its
September22 date and version are verified in the arXiv record. Exact
Definition4.1, Lemma4.3, Lemma4.5, Theorem4.6 and the complete at-most-two-
prime proof were read, as were Example4.11 and the target discussion.
The source's other long sequence/CSP classifications were not fully audited.
Bounded October27 title/author/arXiv-number/Conjecture4.10/divisibility searches
found no verified later proof. This does not establish exhaustive novelty,
first priority, current external acceptance or that every equivalent
argument was absent from prior literature.

Decision: **freeze this full-target proof candidate for independent review**,
then rotate the breadth scout. The review should check recursive prime-step
factor propagation, the restriction to t_k=gcd(k,t), the nontrivial-root/
proper-divisor equivalence, the zero-factor case and the final monic return.
Do not label the source conjecture solved until that scrutiny succeeds.
