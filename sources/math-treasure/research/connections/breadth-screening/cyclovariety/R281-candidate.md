# R281: cyclotomic torus intersection — candidate for independent scrutiny

T026/R281; 2026-10-26; breadth scout. **Stable proof candidate, not independently
reviewed, not labelled solved or novel.** The screening question has a concrete
advance path, so sustained verification is handed to root and the depth peer.

## Exact primary target and definitions

Cornelissen–Hokken–Ringeling, *The asymptotic Mahler measure of Gaussian
periods*, [arXiv:2507.09303v3](https://arxiv.org/pdf/2507.09303v3), 19 June 2026,
61 PDF pages; published as [PLMS 133 (2026), e70198](https://doi.org/10.1112/plms.70198),
25 August 2026. Definition A on **PDF page 3** and the computing description
on **PDF page 15** define

    L=Z[ζ_k], d=φ(k), basis (1,ζ_k,...,ζ_k^(d−1)),
    A_k={coefficient vectors of ζ_k^j : 0≤j<k}⊂L≅Z^d,
    F_k(x)=Σ_(a∈A_k) x^a,  P_k(x0,x)=x0+F_k(x),
    C_k={P_k=0}⊂(C*)^(d+1).

**A_k includes ALL k-th roots, not only primitive roots.** This is explicit
in Definition A and essential: multiplication by ζ_k permutes this set.
T^d denotes the unit complex torus, equivalently Hom(L,S¹).

**Section 10(iii), PDF page 56**, asks whether C_k(C)∩T^(d+1) is nonempty
for every k. Section 5.3, **PDF pages 23–24**, proves listed subfamilies and
states a belief in the general case. The source itself is current primary
evidence of an open question as of August 2026. Searches on 2026-10-26 by
cyclovariety, torus intersection, non-expansiveness, the authors and paper
title found no verified later resolution; that is a bounded search, not an
exhaustive current-status or novelty determination.

The proposed argument below proves the stronger assertion F_k(x)=0 for some
x∈T^d whenever k>1. A path then gives the desired |F_k|=1 point; **the F_k=0
point by itself has x0=0 and is NOT a point of C_k inside the multiplicative
torus**. That return step is kept explicit.

## General symmetry/sign lemma

Let L be a free abelian group of finite rank, U∈Aut(L), and assume U−I is
invertible over Q. Let A⊂L be finite, preserved by U, with no element in 2L.
Let

    F(χ)=Σ_(a∈A) χ(a), χ∈T=Hom(L,S¹),

and assume A is nonempty. Then F has a zero in T.

The coefficients are real (indeed all one). No real-valuedness of F on the
whole torus is assumed. In particular the proof does not apply a real
intermediate value theorem directly to a complex-valued path.

### Proof, with the action and logarithm lift expanded

Put V=Hom(L,R), Λ=Hom(L,Z), so T=V/Λ via
χ_t(a)=exp(2πi t(a)). Define U* on V and Λ by

    (U*t)(a)=t(Ua).

This is the toral-action convention used throughout. Since U permutes A,
f(t):=F(χ_t) satisfies f(U*t)=f(t). Also f(−t)=conj(f(t)), f(0)=|A|>0,
and f(t+η)=f(t) for η∈Λ.

Suppose for contradiction that f never vanishes. It is smooth, being a finite
sum of exponentials. A continuous logarithm on V is given explicitly by

    g(t)=log|A| + ∫_0^1 [d/ds f(st)]/f(st) ds.           (1)

For each t, the integral exponentiates to f(t)/f(0): both functions solve
the scalar ODE y′(s)=(d f(st)/ds)/f(st)·y(s), with initial value 1.
Thus exp(g(t))=f(t). Continuity follows from the continuous integrand and
nonvanishing on each compact radial family. No globally chosen argument on
the torus has yet been assumed.

For η∈Λ, exp(g(t+η)−g(t))=1. The difference is continuous in t, takes values
in the discrete set 2πiZ, and hence is constant on V. Write it as
2πiν(η). Comparing two successive translations proves ν:Λ→Z is additive.
This is the winding homomorphism, obtained here directly from (1).

Similarly g(U*t)−g(t)∈2πiZ is constant; it is zero at t=0, so

    g(U*t)=g(t),       ν(U*η)=ν(η).                     (2)

The first equality and a translation by η prove the second. Since U*−I is
invertible over Q, its image on Λ has finite index. Equation (2) says ν
vanishes on that image. A homomorphism from Λ to torsion-free Z vanishing
on a finite-index subgroup is zero: for any η, some positive m has
mη∈(U*−I)Λ, so mν(η)=0. Therefore ν=0 and g is Λ-periodic, defining a
continuous logarithm G on T.

Conjugation gives exp(g(−t))=exp(conj(g(t))). Their difference lies in
2πiZ, is continuous, and is zero at t=0 because g(0)=log|A| is real.
Consequently g(−t)=conj(g(t)). If χ is a two-torsion sign character, choose
t with 2t∈Λ. Then −t differs from t by a lattice vector, so periodicity
and conjugation imply g(t)=conj(g(t)). Hence G(χ) is real and

    F(χ)=exp(G(χ))>0 for every χ∈Hom(L,{±1}).           (3)

But the exact average of F over this finite sign torus is zero. For each
a∉2L, at least one coordinate is odd, so pairing characters that differ in
that coordinate cancels χ(a). Summing over A gives

    2^(−rank L) Σ_(χ∈Hom(L,{±1})) F(χ)=0.              (4)

This contradicts strict positivity in (3). The nonvanishing assumption is
false, proving the lemma. ∎

## Cyclotomic application and exact return

For k>1 let U:L→L be multiplication by ζ_k. It is an integral automorphism:
its inverse is multiplication by ζ_k^(k−1). The endomorphism U−I is
invertible over Q since multiplication by the nonzero field element ζ_k−1
is invertible on Q(ζ_k). The set A_k of **all** roots is a single U-orbit.
Every a∈A_k is primitive in L: a=U^j(1), and integral automorphisms preserve
primitivity of the basis vector 1. In particular a∉2L. All hypotheses of the
lemma therefore hold, giving x*∈T^d with F_k(x*)=0.

Choose any continuous path γ:[0,1]→T^d from (1,...,1) to x*. Such paths
exist by path connectedness of the torus (choose real arguments for its
endpoint coordinates and interpolate them). The **real** continuous function
s↦|F_k(γ(s))| takes values k>1 at 0 and 0 at 1. The real intermediate value
theorem gives s with |F_k(γ(s))|=1. Put x=γ(s) and x0=−F_k(x). Then every
coordinate, including x0, has modulus one and is nonzero, and P_k(x0,x)=0.
This is the full return to the precise geometric question.

For k=1, U−I=0 so the lemma is inapplicable. Here F₁(x₁)=x₁ and C₁ has
the explicit unit-torus point (x0,x1)=(−1,1). Thus k=1 is handled separately.

Only the geometric intersection question is asserted as a proof candidate.
The paper's equivalence with non-expansiveness of an algebraic action requires
its own primary dependency audit before importing that dynamical conclusion.
The paper's other Mahler-measure, reflexivity and period questions are untouched.

## Archive bridge, controls, attribution and handoff

The archive's [reviewed subset/radix character products](../../../sum-products/README.md)
and [root-of-unity cancellation proof](../../../monster-limit/theta-proof.md)
motivate keeping finite character averages exact. The decisive additional
ingredient here is the symmetry-forced trivial winding class and conjugation
normalization; this is a new deduction in the present screening, with priority
unestablished. It is not a theorem attributed to Ishan's original notes.
The source's lattice, cyclovariety and posed question belong to the cited authors.

The [standard-library exact probe](../../../../experiments/connections/breadth-screening/cyclovariety/verify.py)
checks k=2,...,256: integral cyclotomic reduction, orbit closure and distinctness,
the cyclic basis, Φ_k(1)≠0, primitivity/nonzero mod2, and zero sum of orbit
vectors. It averages entire sign cubes when φ(k)≤12 and retains explicit
unit-valued sign witnesses for k=105,165,195, beyond the source's stated
totient condition; these individual witnesses directly give torus points.
The universal proof is the written argument, not the finite range. Exact
witness and matrix conventions are in the saved verification record.

Independent review priorities: verify the source's ALL-roots definition,
action U* and lattice translation convention, (1)'s logarithm lift and
continuity, invariance forcing ν=0, conjugation at sign characters, primitive
orbit argument and the |F|=1 return. Recheck current status and novelty by
equivalent finite symmetry/Laurent torus formulations before any solution
claim. Root and `/root/open_connections_depth` have this stable candidate;
the scout now hands sustained scrutiny to them and rotates to other material.
