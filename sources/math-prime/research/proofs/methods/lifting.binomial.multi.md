# lifting.binomial.multi: two-index finite-polynomial zero count
Parent: [lifting.binomial](lifting.binomial.md)     Children: none yet

**Status: WRITTEN (2026-10-15).** The full proof, as numbered statements, is
[proof.md](lifting.binomial.multi/proof.md). Two independent in-session referee passes found no counterexample to any
numbered statement, and their corrections are applied; there has been no audit outside the session.

The lemma is proved under a hypothesis that is checked on finitely many coefficients of each orbit. Forcing that
hypothesis a priori is OPEN, and no determinant identity can do it.

## Statement
- **Setting.** Fix a prime p and q = p^e, with e ≥ 1, and e ≥ 2 when p = 2. For integers b_ij put
  f(m,n) = Σ_{i≤m, j≤n} C(m,i) C(n,j) q^{i+j} b_ij. Equivalently f(m,n) = ℓT₁^mT₂^n v for commuting integer matrices
  T_l = I+qK_l, with b_ij = ℓK₁^iK₂^j v.
- **Stable reduction.** Let val(i,j) = ν_p(q^{i+j}b_ij) − ν_p(i!) − ν_p(j!), let s* be its minimum and E the set of
  minimizers. The stable reduction of f is

  R̄ = Σ_{(i,j)∈E} ū_ij X^(i) Y^(j) ∈ F_p[X,Y],   with X^(i) = X(X−1)⋯(X−i+1) and ū_ij the unit part of q^{i+j}b_ij/(i!j!) mod p.

  It is the reduction of every truncation G_N = (N!)² Σ_{i,j<N} q^{i+j}b_ij C(X,i)C(Y,j) once N > max over E of
  max(i,j), after dividing by its exact p-content 2ν(N!)+s* (Lemmas 3–4). Its total degree is at most
  max(0, ⌊(s*(p−1)−1)/(e(p−1)−1)⌋).
- **Theorem 8.** If the stable reductions R̄₁, R̄₂ of f₁, f₂ have no common factor, of total degrees d₁, d₂, then
  #{(m,n) ∈ N² : f₁ = f₂ = 0} ≤ Σ_P mult_P(R̄₁,R̄₂) ≤ dim F_p[X,Y]/(R̄₁,R̄₂) ≤ d₁d₂. The sum is over the F_p-rational
  common zeros P of R̄₁, R̄₂. A gcd form counts the zeros off the common factor.
- **Theorem 10.** The same bound holds on Z². A shift of the seed translates the stable reductions:
  R̄_{T^{(m₀,n₀)}v} = R̄_v(X+m₀, Y+n₀). This needs T_l invertible over Z_(p), which T_l = I+qK_l always is.
- **Theorem 11 (saturation).** Combinations (A₁F₁+A₂F₂)/p^τ of the truncations may be added to the ideal, under a
  stated precision condition. The bound then uses the larger ideal.
- **Theorem 18 (frame criterion; 4×4 matrices, two markers independent mod p).** Let
  Q_v(c) = det[v, K₁v, K₂v, (c₁K₁+c₂K₂)²v]. Suppose ν_p(Q_v(c)) ≤ s_F for every c with p ∤ gcd(c₁,c₂), and
  e ≥ s_F + 1 + [p ≤ 3]. Then the orbit of v has at most two common zeros.
  - In a field, Q_β = N(β)·Q₁, so the bound holds for every orbit whose seed has norm prime to p.
  - The proof is one explicit saturation step.
  - The term [p ≤ 3] is necessary at p = 3: there is an explicit system with three zeros at e = s_F + 1.
  - At p = 2 one always has s_F ≥ 1, so q ≥ 8 is needed.

## Proof idea
It follows [lifting.binomial](lifting.binomial.md). No completion and no power series is used.
1. **Exact content.** Finite differences at (0,0) recover q^{i+j}b_ij. For an integer polynomial they are
   divisible by i!j!, by Stirling numbers. This gives the exact content of G_N, which in one index is the paper's
   remark that R(2)−2R(1)+R(0) is even.
2. **Degree bound.** Legendre's formula ν_p(i!) ≤ (i−1)/(p−1) makes every term outside E divisible by an extra p.
   So R̄ has bounded degree and is stable in N.
3. **Colength.** For a finite set S of common zeros, choose N after S. The evaluation image of Z_(p)[X,Y] on S is free
   of rank |S|, so F_p[X,Y]/(R̄₁,R̄₂) maps onto an algebra of dimension |S|. In one index this replaces "the monic
   product Π(X−n) divides R". The resultant step planned in GENERALIZATIONS.md §3.2 is not needed.
4. **Affine Bézout.** Two coprime polynomials of degrees d₁, d₂ over F_p have quotient of dimension at most d₁d₂. The
   proof is a dimension count in degree ≤ t.
5. **Signed indices.** The Vandermonde identity for binomials, with the high-order tail controlled, moves the seed to
   the corner of any finite set of zeros, as the one-index proof does.

## The determinant identity, and where the plan of §3.2 fails
- **The identity.** The plan in [GENERALIZATIONS.md](../GENERALIZATIONS.md) §3.2 asked for a determinant identity
  that would force the coprimality of R̄₁ and R̄₂, as Δ₀ ≠ 0 did in one index. The analogue is
  det[w, K₁w, K₂w, K₁K₂w] = det(T₁)^m det(T₂)^n det[v, K₁v, K₂v, K₁K₂v]. In a quartic field this equals
  N(β)·det[1, η₁, η₂, η₁η₂] when ε_l = 1+qη_l.
- **What it forces (Proposition 13).** It forces bounded content of each marker function. It bounds the depth of any
  constant-coefficient dependence between them. And at every common zero, one of three 2×2 minors of the marker jets
  (the Jacobian J and two minors involving the mixed column) has valuation ≤ ν(det).
- **What it does not force.** It forces neither coprimality nor finiteness.
  - In A = Z[x,y]/(x²,y²) the identity holds with valuation 0 at every orbit point. Yet the markers f₁ = qm and
    f₂ = q²mn have infinitely many common zeros, and R̄₁ = X, R̄₂ = XY share a factor (Example 14).
  - In Z[2^{1/4}] at p = 7, 16 of 47 seeds x−y·2^{1/4} have stable reductions with a common linear factor,
    although ν₇ of the determinant is 0 there ([necessity_examples.out.txt](lifting.binomial.multi/necessity_examples.out.txt)).
    ν₇(J) is 1 or 2 at their zeros.
  - Lido–Mauri §4 show the same for a quintic Thue equation.
- **Why.** The reductions depend on which (i,j) minimize val. The identity controls only the mixed column K₁K₂w, not
  K₁²w or K₂²w. When the linear parts at a zero are proportional mod p (p | J), coprimality fails.
- **What survives.**
  - Saturation (Theorem 11), the finite-precision form of Lido–Mauri's saturated ideal. One step is not always
    enough: 59 of 167 non-coprime synthetic cases needed more.
  - The frame family Q_v(c). It controls the second-order term in the one bad direction, and gives Theorem 18.
  - The Jacobian criterion (Proposition 17): a zero where J is not too divisible by p is the only zero in its
    orbit. This is Skolem's 1934 criterion in finite form.

## Tests
Every script was re-run by an independent verifier, and an independent re-implementation
([independent_check.py](lifting.binomial.multi/independent_check.py) →
[.out.txt](lifting.binomial.multi/independent_check.out.txt)) reproduces every per-instance result.
- **Synthetic** ([synthetic_tests.py](lifting.binomial.multi/synthetic_tests.py) →
  [.out.txt](lifting.binomial.multi/synthetic_tests.out.txt)): 328 instances, p ∈ {2,3,5,7}, d = 4…8.
  - 161 have coprime reductions, with 154 common zeros in all and 0 violations. The bound is attained in 71
    instances.
  - Stability in N and translation invariance hold in every instance tested.
  - One saturation step repairs 108 of the 167 non-coprime cases.
- **Necessity** ([necessity_examples.py](lifting.binomial.multi/necessity_examples.py) →
  [.out.txt](lifting.binomial.multi/necessity_examples.out.txt)). With T₂ = T₁ there are 21 common zeros on a line
  against d₁d₂ = 3, so the hypothesis is needed. The natural non-coprime seeds in Z[2^{1/4}] at p = 7 are each
  repaired by one saturation step, to colength 2.
- **Unit example** ([quartic_units_example.py](lifting.binomial.multi/quartic_units_example.py) →
  [.out.txt](lifting.binomial.multi/quartic_units_example.out.txt)): the orbit of ±1 under ⟨1+2^{1/4}, 1+√2⟩ in
  Z[2^{1/4}], i.e. x⁴−2y⁴ = ±1.
  - With cosets of h₁Z×h₂Z the per-prime sums are 8, 240 and 6 at p = 3, 5, 7, against 3 known solutions. That
    decomposition is redundant, because of −1 and relations between the units mod p.
  - With cosets of G₇ = {(a,b) : ε₁^aε₂^b ≡ ±1 mod 7}, the prime 7 alone gives a sharp sum of 3, using only the
    multiplicity-1 case (independent_check, part 5).
  - This is not a proof for the equation: neither the unit index nor the coverage by orbits was checked.
- **The proof's own checks** ([lemma_checks.py](lifting.binomial.multi/lemma_checks.py) →
  [.out.txt](lifting.binomial.multi/lemma_checks.out.txt)) cover the content and reduction formulas (1784
  truncations), the one-index specialization, Examples 14–15, and planted tests of Theorem 18 with a negative control
  at the threshold.
- **Referees' checks.** [referee_checks.py](lifting.binomial.multi/referee_checks.py) →
  [.out.txt](lifting.binomial.multi/referee_checks.out.txt) and
  [adversarial_checks.py](lifting.binomial.multi/adversarial_checks.py) →
  [.out.txt](lifting.binomial.multi/adversarial_checks.out.txt). They include 738 planted certificates for
  Theorem 18, a search for three dependent orbit points under its hypotheses (none found), and the explicit
  three-zero system at p = 3, e = s_F+1.
- **The S4 order Z[ω], ω⁴ = ω+1** (signature (2,1), units ω and ω−1).
  - The frame form is anisotropic mod p at p = 5, 11, 19, 23, 31. After refinement the hypothesis also holds at
    p = 13 and 17.
  - At p = 5, every coset of ⟨ω³¹², (ω−1)¹⁰⁴⟩ in the unit group contains at most two solutions of
    x⁴−xy³−y⁴ = ±1.

## Requirements
- Two markers in two indices. This is the case of a quartic Thue equation whose field has one complex place,
  signature (2,1).
- A prime p with T_l ≡ I mod q. In a field, take powers of the units.
- One of three finite checks on the seed: coprime stable reductions (Theorem 8), a saturation certificate
  (Theorem 11), or the frame hypothesis (Theorem 18).

## Limits
- **The hypothesis is checked per orbit, not forced.** Example 14, and the natural Z[2^{1/4}] seeds, show that no
  argument using only the determinant identity can force it.
- **Theorem 18 is in practice an S4 tool.** A quartic field of signature (2,1) is D4 or S4. A D4 field contains a
  real quadratic subfield, whose unit makes the frame form vanish in its direction when it is a generator, as 1+√2
  does in Z[2^{1/4}]. Theorem 11 and Proposition 17 still apply zero by zero.
- **It counts zeros only.** It does not locate them, and it does not cover the solutions by orbits. Coverage is the
  descent's job: [ghost.unit.multi_unit](ghost.unit.multi_unit.md).
- **Two markers in two indices only.** More markers or more indices (quintic Thue equations with 3 markers in 2
  indices; 3 or more indices) are not treated.
- **The four rank-2 targets are out of reach.** None has a Thue field with a complex place
  ([rank2_signatures_2026-10.md](../targets/rank2_signatures_2026-10.md)). The comparison curves 1611.a1, 8604.a1
  and y²=x³+17 do. On all three the lemma was used (2026-10-15): Proposition 17 with a congruence sieve sufficed for
  all 32 equations (twelve each for 1611.a1 and 8604.a1, eight for y²=x³+17), and Theorem 18 and saturation were not
  needed. One instance on y²=x³+17 has λ = 1 (q = 8).

## Relation to Lido–Mauri
G. M. Lido and L. Mauri, *A multivariate Strassmann theorem*, arXiv:2605.03683 (2026).
- **Their results.** For power series f₁,…,f_r in n variables over a complete valuation ring, Theorem 1.6 bounds the
  common zeros by the sum of the local dimensions of k[x]/(f̄₁,…,f̄_r) at maximal ideals with residue field k.
  Theorem 2.5 characterizes finiteness of the zeros over the completed algebraic closure through the reduction of
  the saturated ideal.
- **What Theorem 8 is.** It is their Theorem 1.6 with n = r = 2 and R = Z_p, restricted to integer zeros. For two
  polynomials in two variables, coprimality is the same hypothesis as finite dimensionality. Theorem 11 is the
  sufficient direction of their Theorem 2.5, at finite precision. Example 15 is their Example 3.2.3 in this setting.
- **What they have that we lack.** General n and r; zeros in R^n; necessity; algorithms for saturation from
  approximate data; a solved quintic Thue equation.
  - The series printed in their §4 contain errors.
  - Their stated reduction and their per-class bound of 2 are correct by recomputation
    ([verify_literature_checks.out.txt](../targets/rank2_signatures_2026-10/verify_literature_checks.out.txt)).
- **What this card adds.**
  - Only finite integer polynomials are used, with the cutoff chosen after the zeros, as the project's standard asks.
  - The content and the stable reduction are explicit, and so is the translation lemma for signed orbits.
  - The frame criterion (Theorem 18) gives an a priori one-step certificate. Nothing like it appears in their text;
    the wider Skolem literature was not searched for second-order criteria.

## Theorems using it
- [triangular_tetrahedral](../theorems/triangular_tetrahedral/README.md), Theorem A.
- [triangular_square_pyramidal](../theorems/triangular_square_pyramidal/README.md), Theorem A.
- [mordell_k_minus17](../theorems/mordell_k_minus17/README.md), Theorem A.

All three use Proposition 17 only, on the classes of kernel lattices of the units, with a congruence sieve.

## Provenance
- **Novelty.** None is claimed. The counting principle is classical:
  - Strassmann (one variable);
  - Skolem, from 1933, with congruence and several-variable forms; Proposition 17 is his 1934 criterion, as stated
    by Lombardo, arXiv:2208.03830, Theorem 1.7;
  - Ljunggren (1942), who used a two-index Skolem argument;
  - Edixhoven–Lido (2023), Theorem 4.12;
  - Lido–Mauri (2026).

  What is new here is the form: finite polynomials only.
- **Sources.**
  - The literature: [rank2_signatures_2026-10/literature.md](../targets/rank2_signatures_2026-10/literature.md).
  - The one-index original: [main.tex](../../shared/prior/elementary_diophantine_paper/main.tex), section "A finite
    lifting argument for marked returns" and subsection "A sharper count of all marked returns".
