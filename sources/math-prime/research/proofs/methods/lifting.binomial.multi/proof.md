# lifting.binomial.multi: the two-index finite-polynomial return lemma

Author: Claude (Opus 5.5), workflow task B1, 2026-10-15.

Refereed the same day by two independent in-session passes. Their checks are in
[referee_checks.py](referee_checks.py) → [referee_checks.out.txt](referee_checks.out.txt) and
[adversarial_checks.py](adversarial_checks.py) → [adversarial_checks.out.txt](adversarial_checks.out.txt). No
counterexample to any numbered statement was found, and their corrections are applied below. There has been no
audit outside this session.

Checks: [lemma_checks.py](lemma_checks.py) with its output [lemma_checks.out.txt](lemma_checks.out.txt)
(all ten sections PASS; 20–45 s, depending on the machine). Task B2's scripts in this folder (`binomial2.py`, `synthetic_tests.py`,
`necessity_examples.py`, `quartic_units_example.py`) test the same statements; their definitions agree
with §1 (their stabilization index `N0` is the index N_E of Lemma 4).

No novelty is claimed. The analytic statement is classical: Strassmann (one variable), Skolem's method,
and the multivariate Strassmann theorem of Lido and Mauri (2026). What is written here is a form of it
that uses only finite integer polynomials, with the cutoff chosen after the finite set of zeros, as in
the paper's one-index lemma.

## 0. Summary

1. **Proved.** Let f₁, f₂ be two-index binomial sequences (§1). Suppose their stable reductions R̄₁, R̄₂
   are coprime in F_p[X,Y], of total degrees d₁, d₂. Then the common zeros in N² number at most
   Σ_P mult_P ≤ dim F_p[X,Y]/(R̄₁,R̄₂) ≤ d₁d₂ (Theorem 8). The same holds on Z² when T₁, T₂ are
   invertible over Z (Theorem 10). Two further forms need no global coprimality: a local form and a
   gcd form (Theorem 8(a),(b)). A saturated form adds certificates (aR₁+bR₂)/p^t, under an explicit
   precision condition (Theorem 11). The content is exact, k_N = 2ν(N!) + s*_N, which equals 2ν(N!) + s* once
   N ≥ N_E (Lemma 3). The stable reduction is the explicit polynomial Σ_{(i,j)∈E} ū_ij X^(i) Y^(j) (Lemma 4).
2. **The determinant identity.** det[w, K₁w, K₂w, K₁K₂w] = det(T₁)^m det(T₂)^n det[v, K₁v, K₂v, K₁K₂v];
   in a quartic field it equals N(w)·D₀ = ±N(β)·D₀ with D₀ = det[1, η₁, η₂, η₁η₂]. It forces the
   three properties (1)–(3) of the plan (Proposition 13). It does **not** force coprimality, nor even
   finiteness: Example 14 satisfies the identity with s = 0 at every orbit point and has infinitely
   many common zeros. Example 15 shows that coprimality of the raw reductions also fails where the
   count is finite.
3. **What survives.** Use the family of frames det[w, K₁w, K₂w, (c₁K₁+c₂K₂)²w]. Their value is a
   binary quadratic form Q(c). Suppose ν_p(Q(c)) ≤ s_F for every primitive c, and
   e ≥ s_F + 1 + [p ≤ 3]. Then the Z²-orbit of v (T_l ∈ GL₄(Z)) has at most **2** common zeros, and in
   general at most 2 lie in the quadrant above any common zero (Theorem 18). In a field this covers every
   orbit whose seed β has p ∤ N(β), since Q_β = N(β)·Q₁. The proof is one explicit saturation step,
   which contains a finite Jacobian criterion (Proposition 17). With the refinement T → T^{p^a} the
   bound becomes 2p^{2a} (Corollary 19).
4. **Examples.** Take the S4 quartic order Z[ω], ω⁴ = ω + 1, of signature (2,1), with units ω and ω − 1.
   - The frame form is anisotropic mod p (s_F = 0) at p = 5, 11, 19, 23, 31.
   - At p = 13 and 17 the hypothesis holds with s_F = 1 after a refinement (second referee; §8.5).
   - For Z[2^{1/4}] with units 1+ω, 1+ω², the hypothesis fails at every prime, because 1+ω² lies in the
     subfield Q(√2).

## 1. Setting and notation

1.1. p is a prime and q = p^e, with e ≥ 1 if p is odd and e ≥ 2 if p = 2. Write ν = ν_p, with
ν(0) = +∞, and Z_(p) for the rationals with denominator prime to p. Put ι := ν(2) (1 if p = 2,
else 0) and ι' := 1 if p ∈ {2,3}, else 0. The content of a polynomial with rational coefficients is
the least ν of its coefficients in the monomial basis.

1.2. b_ij (i,j ≥ 0) are integers. Every statement below holds verbatim for b_ij ∈ Z_(p), because the
proofs use only ν(b_ij) ≥ 0; then G_N and R_N lie in Z_(p)[X,Y] instead of Z[X,Y]. This is needed for
combinations such as f₂ − ρf₁ with ρ ∈ Z_(p).

   f(m,n) = Σ_{i≤m, j≤n} C(m,i) C(n,j) q^{i+j} b_ij      (m,n ≥ 0).

Matrix form: T_l = I + qK_l are commuting integer matrices, l is a row and v a column, and
b_ij = l K₁^i K₂^j v. The finite binomial theorem (I commutes with K_l) gives l T₁^m T₂^n v = f(m,n).

1.3. X^(i) := X(X−1)⋯(X−i+1), with X^(0) = 1; C(X,i) := X^(i)/i!. The truncations are

   F_N(X,Y) = Σ_{i,j<N} q^{i+j} b_ij C(X,i) C(Y,j),   G_N = (N!)² F_N.

Since C(m,i) = 0 for integers 0 ≤ m < i, F_N(m,n) = f(m,n) for 0 ≤ m,n < N. Also
G_N = Σ (N!/i!)(N!/j!) q^{i+j} b_ij X^(i) Y^(j) ∈ Z[X,Y]. Write G_N = p^{k_N} R_N with R_N not
identically 0 mod p, and R̄_N := R_N mod p.

1.4. val(i,j) := ν(q^{i+j} b_ij) − ν(i!) − ν(j!) for b_ij ≠ 0, and +∞ for b_ij = 0. Then
s* := min val(i,j), and s*_N := min_{i,j<N} val(i,j). E := {(i,j) : val(i,j) = s*} is the set of
minimizers, D* := max_{(i,j)∈E} (i+j), and N_E := 1 + max_{(i,j)∈E} max(i,j).

1.5. θ := e(p−1) − 1, an integer ≥ 1 (p odd: p−1 ≥ 2; p = 2: e ≥ 2). For n ≥ 1 put
Λ(n) := ⌈(nθ+1)/(p−1)⌉. Let D(s) := max(0, ⌊(s(p−1)−1)/θ⌋) and N₀ := D(s*) + 1. Explicitly
Λ(1) = e and Λ(2) = 2e − ι. Λ(3) is 3e − 2 for p = 2, 3e − 1 for p = 3 and 3e for p ≥ 5, so
Λ(3) ≥ 3e − ι − ι'.

1.6. A *stable reduction* is R̄_N for all large N; these agree up to a nonzero scalar (Lemma 4). The
normalized stable reduction R̄ is defined in Lemma 4.

## 2. Content and stable reduction

**Lemma 1 (growth).** (a) For i ≥ 1, ν(i!) ≤ (i−1)/(p−1). (b) For (i,j) ≠ (0,0) with b_ij ≠ 0,
val(i,j) ≥ ((i+j)θ + 1)/(p−1), hence val(i,j) ≥ Λ(i+j). Λ is nondecreasing and unbounded. (c) Every
val(i,j) is a nonnegative integer. If some b_ij ≠ 0, then s* exists, E is finite, and
E ⊆ {i+j ≤ D(s*)}. Hence D* ≤ D(s*) and N_E ≤ D* + 1 ≤ N₀. (d) If s̃ := min_{max(i,j)<N} val(i,j)
satisfies s̃ < Λ(N), then s* = s̃. So s* and E are determined by the b_ij with max(i,j) < N.

*Proof.*
1. Legendre: ν(i!) = Σ_{r≥1} ⌊i/p^r⌋ (count the multiples of p, of p², …). If i = Σ e_r p^r in base p,
   this sum equals (i − Σ e_r)/(p−1), and Σ e_r ≥ 1 for i ≥ 1. This gives (a).
2. For (i,j) ≠ (0,0), (a) gives ν(i!) + ν(j!) ≤ (i+j−1)/(p−1). This is (i−1+j−1)/(p−1) when both are
   ≥ 1, and (i+j−1)/(p−1) when one is 0. Since ν(b_ij) ≥ 0,
   val(i,j) ≥ (i+j)e − (i+j−1)/(p−1) = ((i+j)θ + 1)/(p−1). val is an integer, so val(i,j) ≥ Λ(i+j).
   This gives (b), because θ ≥ 1.
3. val(0,0) = ν(b₀₀) ≥ 0, and the other values are ≥ Λ(1) = e ≥ 1. A nonempty set of nonnegative
   integers has a least element, so s* exists. If (i,j) ∈ E and (i,j) ≠ (0,0), then Λ(i+j) ≤ s*,
   that is (i+j)θ + 1 ≤ s*(p−1), so i+j ≤ D(s*). This gives (c).
4. Every (i,j) with max(i,j) ≥ N has i+j ≥ N, so val(i,j) ≥ Λ(N) > s̃. This gives (d). ∎

**Lemma 2 (Stirling divisibility).** Let g ∈ Z_(p)[X,Y]. Then g = Σ c_ij C(X,i) C(Y,j) uniquely over
Q, with c_ij = Δ_X^i Δ_Y^j g(0,0) ∈ i! j! Z_(p). Here Δ_X g(X,Y) = g(X+1,Y) − g(X,Y).

*Proof.*
1. The polynomials C(X,i)C(Y,j) have distinct leading monomials X^iY^j, so they form a Q-basis of
   Q[X,Y].
2. Pascal's rule gives Δ_X C(X,i) = C(X,i−1) for i ≥ 1, and Δ_X C(X,0) = 0. Hence Δ_X^a C(X,i) = C(X,i−a)
   for a ≤ i, and 0 for a > i. At X = 0 this is 1 if a = i and 0 otherwise. Applying Δ_X^iΔ_Y^j to
   the expansion and evaluating at (0,0) gives c_ij = Δ_X^iΔ_Y^j g(0,0).
3. From X·X^(i) = X^(i+1) + i·X^(i), induction on a gives X^a = Σ_{i≤a} S(a,i) X^(i). The S(a,i) are
   nonnegative integers (Stirling numbers of the second kind), with S(a+1,i) = S(a,i−1) + i·S(a,i).
   Hence X^aY^b = Σ_{i,j} S(a,i) S(b,j) i! j! C(X,i) C(Y,j).
4. For g = Σ g_ab X^aY^b with g_ab ∈ Z_(p), linearity gives c_ij = i! j! Σ g_ab S(a,i) S(b,j). ∎

In one index with i = 2 this is the paper's remark that R(2) − 2R(1) + R(0) is even (main.tex, proof of
Theorem three-return).

**Lemma 3 (exact content).** Let N ≥ 1 and suppose some b_ij ≠ 0 with i,j < N. Then content(F_N) =
s*_N and k_N = 2ν(N!) + s*_N. For N ≥ N_E, k_N = 2ν(N!) + s*.

*Proof.*
1. Each term t_ij := q^{i+j} b_ij C(X,i) C(Y,j) equals [q^{i+j} b_ij/(i! j!)]·X^(i)Y^(j). The polynomial
   X^(i)Y^(j) has integer coefficients and coefficient 1 at X^iY^j, so content(t_ij) = val(i,j)
   exactly. Hence content(F_N) ≥ s*_N.
2. Let (i₀,j₀) attain s*_N, and suppose content(F_N) ≥ s*_N + 1. Then p^{−s*_N−1}F_N ∈ Z_(p)[X,Y].
   Its binomial coefficient at (i₀,j₀) is p^{−s*_N−1} q^{i₀+j₀} b_{i₀j₀}. By Lemma 2 it lies in
   i₀! j₀! Z_(p), so val(i₀,j₀) ≥ s*_N + 1. This is a contradiction.
3. G_N = (N!)²F_N ≠ 0, by the basis property, and content(G_N) = 2ν(N!) + content(F_N) = k_N. For
   N ≥ N_E, E meets [0,N)², so s*_N = s*. ∎

Step 2 gives k_N − 2ν(N!) ≤ s*_N, and step 1 gives equality. The plan's inequality k_N − 2ν(N!) ≤ s*
holds only for N ≥ N_E. For example, with p = 5, e = 1, b₀₀ = 5³ and b₂₀ = 1, we have s* = 2 and N_E = 3,
but k₂ − 2ν(2!) = 3 (first referee, referee_checks §B).

**Lemma 4 (stable reduction).** Assume some b_ij ≠ 0. For (i,j) ∈ E put
u_ij := q^{i+j} b_ij /(p^{s*} i! j!), an element of Z_(p)^×. Define

   R̄ := Σ_{(i,j)∈E} ū_ij X^(i) Y^(j) ∈ F_p[X,Y],

the normalized stable reduction. Then:
(a) For N ≥ N_E, p^{−s*}F_N ∈ Z_(p)[X,Y] reduces to R̄. Moreover R_N = u_N·p^{−s*}F_N, where
    u_N := (N!/p^{ν(N!)})² is an integer prime to p, so R̄_N = ū_N R̄.
(b) R̄ ≠ 0. Its total degree is D* ≤ D(s*), and its top-degree part is Σ_{(i,j)∈E, i+j=D*} ū_ij X^iY^j.
(c) R̄ is determined by the b_ij with i+j ≤ D(s*), through their valuations and their unit parts mod p.

*Proof.*
1. p^{−s*}t_ij = [q^{i+j} b_ij/(p^{s*} i! j!)]·X^(i)Y^(j). The scalar has valuation val(i,j) − s* ≥ 0,
   with equality exactly for (i,j) ∈ E. So p^{−s*}t_ij ≡ 0 mod p for (i,j) ∉ E, and
   p^{−s*}t_ij ≡ ū_ij X^(i)Y^(j) for (i,j) ∈ E.
2. For N ≥ N_E all of E lies in [0,N)². Summing gives (a), since R_N = G_N/p^{k_N}
   = (N!)² F_N/p^{2ν(N!)+s*} by Lemma 3.
3. X^(i) is monic of degree i, also mod p. Distinct (i,j) have distinct top monomials X^iY^j, so the
   degree-D* part of R̄ cannot cancel. Lemma 1(c) gives D* ≤ D(s*) and part (c). ∎

For N < N_E the reduction may differ, because a minimizer is missing. B2's `binomial2.py` uses exactly
N_E. The a priori index N₀ = D(s*)+1 ≥ N_E needs only s*, p and e.

**Corollary 5 (degree bound; one index).**
(a) For N ≥ N_E, deg R̄_N = D* ≤ D(s*) = max(0, ⌊(s*(p−1)−1)/(e(p−1)−1)⌋). This is independent of N.
(b) In one index, f(n) = Σ C(n,j) q^j b_j and G_N = N! Σ_{j<N} q^j b_j C(X,j). The same proofs give
    k_N = ν(N!) + s*_N, R̄ = Σ_{j∈E} ū_j X^(j), and deg R̄ = max E ≤ D(s*). For p = 2, q = 4 (θ = 1),
    D(s*) = max(0, s*−1). Suppose b_j ≠ 0 and ν(b_j) ≤ s for some j ≤ 2. Then
    s* ≤ val(j) ≤ s + 3, where j = 2 gives 4 + ν(b₂) − ν(2!). So for N ≥ 3,
    k_N ≤ ν(N!) + s + 3. This is the paper's bound (main.tex eq. `return-content`). Also
    deg R̄ ≤ s + 2, which is the paper's count (Theorem `three-return`). Checked in lemma_checks §2,
    with no violation.

## 3. Two algebraic lemmas

**Lemma 6 (colength).** Let S ⊂ Z² be finite and J ⊆ Z_(p)[X,Y] an ideal all of whose elements vanish
on S. Let J̄ ⊆ F_p[X,Y] be its reduction and A := F_p[X,Y]/J̄. For P ∈ F_p² put
m_P := (X−P₁, Y−P₂), S_P := {s ∈ S : s ≡ P mod p}, and

   mult_P(J̄) := sup_k dim_{F_p} F_p[X,Y]/(J̄ + m_P^k) ∈ N ∪ {∞}.

Then:
(a) |S_P| ≤ dim F_p[X,Y]/(J̄ + m_P^k) for every k ≥ 2|S_P|, so |S_P| ≤ mult_P(J̄).
(b) |S| ≤ Σ_{P ∈ S mod p} mult_P(J̄).
(c) |S| ≤ dim A. Moreover, for any ideal I ⊆ F_p[X,Y] with dim F_p[X,Y]/I < ∞,
    Σ_{P∈F_p²} mult_P(I) ≤ dim F_p[X,Y]/I. In particular this holds for I = J̄ when dim A < ∞.

*Proof.*
1. Ev: Z_(p)[X,Y] → Z_(p)^S, g ↦ (g(s))_{s∈S}, is a ring homomorphism. Call its image Ev; it is a
   subring containing 1. Its kernel I_S contains J.
2. *Ev is free of rank |S|.* First, a Z_(p)-submodule of Z_(p)^n is free of rank ≤ n. The image of
   the last coordinate is 0 or p^t Z_(p), because every nonzero element of Z_(p) is p^k times a unit.
   Split off one preimage of p^t and induct on n. Second, choose κ ∈ Z with ℓ := X + κY injective on
   S; each pair s ≠ s' excludes at most one κ. Then g_s := Π_{s'≠s}(ℓ − ℓ(s')) ∈ Z[X,Y] has
   Ev(g_s) = c_s·1_s with c_s ≠ 0, so Ev contains |S| independent vectors.
3. F_p[X,Y]/Ī_S ≅ Z_(p)[X,Y]/(I_S + p) ≅ Ev/pEv ≅ F_p^{|S|}.
4. *Residue-class idempotents.* Fix P with S_P ≠ ∅. Each s' ∈ S∖S_P differs from P in some coordinate.
   Choose ℓ_{s'} ∈ {X − s'₁, Y − s'₂} accordingly, so that ℓ_{s'}(s) ≢ 0 mod p for every s ∈ S_P.
   Put h := Π_{s'∉S_P} ℓ_{s'}. It vanishes on S∖S_P and is a p-adic unit on S_P. Let
   c(T) := Π_{s∈S_P}(T − h(s)) ∈ Z_(p)[T]; then c(0) ∈ Z_(p)^×. The element c(Ev(h)) ∈ Ev has
   coordinate 0 on S_P and c(0) off S_P. Hence the indicator e_P := 1 − c(0)^{−1}c(Ev(h)) of S_P lies
   in Ev.
5. Ev = ⊕_P e_P Ev, and e_P Ev is the evaluation image on S_P, free of rank |S_P| by step 2. So
   B := Ev/pEv = ⊕_P V_P with V_P := e_PEv/p e_PEv and dim V_P = |S_P|. Each V_P is stable under
   multiplication by X and Y.
6. *Nilpotency on V_P.* Fix a lift P̃ ∈ Z² and put r := |S_P|. On e_P Ev, multiplication by X − P̃₁ is
   diagonal with entries s₁ − P̃₁ ∈ pZ_(p). The polynomial χ(T) := Π_{s∈S_P}(T − (s₁−P̃₁)) kills it,
   and χ(T) = T^r + pψ(T) with ψ ∈ Z_(p)[T]. So (X−P̃₁)^r e_PEv = −pψ(X−P̃₁) e_PEv ⊆ p e_PEv, that is
   (X − P₁)^r = 0 on V_P. The same holds for Y. Hence m_P^{2r−1} V_P = 0.
7. *Invertibility on the other classes.* For P' ≠ P with, say, P'₁ ≠ P₁, on V_{P'} we have
   X − P₁ = (X − P'₁) + (P'₁ − P₁): a nilpotent operator plus a nonzero scalar, hence invertible.
   So m_P^k V_{P'} = V_{P'} for all k.
8. By steps 5–7, B/m_P^k B ≅ V_P for k ≥ 2r. Since J̄ ⊆ Ī_S, A maps onto B, so
   A/m_P^kA = F_p[X,Y]/(J̄ + m_P^k) maps onto V_P. This gives (a); summing over P gives (b).
9. (c): |S| = dim B ≤ dim A. Now let I be any ideal with dim F_p[X,Y]/I < ∞, and let P ≠ P'. The
   ideal m_P + m_{P'} contains a nonzero constant, so raising 1 = a + b to the power 2k−1 gives
   1 ∈ m_P^k + m_{P'}^k. By the Chinese remainder theorem,
   F_p[X,Y]/∩_P(I+m_P^k) ≅ Π_P F_p[X,Y]/(I+m_P^k), and the left side is a quotient of F_p[X,Y]/I. So
   Σ_P dim F_p[X,Y]/(I + m_P^k) ≤ dim F_p[X,Y]/I for every k; let k → ∞. ∎

mult_P equals the dimension of the localization at m_P. This is standard and not used here. In one
index, step 3 is the paper's remark that the monic product Π(X − n) divides R.

**Lemma 7 (affine Bézout).** Let F be a field and g₁, g₂ ∈ F[X,Y] nonzero, of total degrees d₁, d₂.
If they have no common factor of positive degree, then dim F[X,Y]/(g₁,g₂) ≤ d₁d₂. If they have one,
the dimension is infinite.

*Proof.*
1. Let V_t be the space of polynomials of degree ≤ t; dim V_t = C(t+2,2). For t ≥ d₁+d₂ consider
   φ_t: V_{t−d₁} × V_{t−d₂} → V_t, (a,b) ↦ ag₁ + bg₂.
2. If ag₁ = −bg₂, then g₂ | ag₁, so g₂ | a. This uses that F[X,Y] is a UFD (Gauss's lemma; named
   fact) and that g₁, g₂ are coprime. Write a = cg₂; then b = −cg₁ and deg c ≤ t − d₁ − d₂, since
   degrees add in a domain. Conversely every such c gives a kernel element. So
   dim ker φ_t = C(t−d₁−d₂+2, 2).
3. dim im φ_t = C(t−d₁+2,2) + C(t−d₂+2,2) − C(t−d₁−d₂+2,2) = C(t+2,2) − d₁d₂. This is a direct
   expansion, checked in lemma_checks §4.
4. im φ_t ⊆ (g₁,g₂) ∩ V_t, so the image of V_t in F[X,Y]/(g₁,g₂) has dimension ≤ d₁d₂. These images
   increase with t and exhaust the quotient, so the quotient has dimension ≤ d₁d₂.
5. If g = gcd has positive degree in Y, then F[X] meets (g) only in 0, so F[X] ↪ F[X,Y]/(g), which is
   a quotient of F[X,Y]/(g₁,g₂). Otherwise g has positive degree in X and F[Y] ↪ F[X,Y]/(g). ∎

In particular, for two polynomials in two variables, coprimality is equivalent to finite colength.
That is the hypothesis of Lido–Mauri's Theorem 1.6 for n = r = 2 (§9).

## 4. The two-index return lemma

**Theorem 8.** Hypotheses:
- (T1) p, q as in 1.1, and two coefficient families b^(1), b^(2) (integers or Z_(p)), giving f₁, f₂.
- (T2) neither f_k is identically zero on N². Equivalently, some b^(k)_ij ≠ 0: by Lemma 2 applied to
  F_N, q^{i+j}b_ij = Δ_m^iΔ_n^j f(0,0).

Let R̄₁, R̄₂ be the stable reductions (Lemma 4), of total degrees d₁, d₂, and let
Ω := {(m,n) ∈ N² : f₁(m,n) = f₂(m,n) = 0}. Then:
(a) For every P ∈ F_p², #{z ∈ Ω : z ≡ P mod p} ≤ mult_P(R̄₁, R̄₂).
(b) Let g := gcd(R̄₁, R̄₂), δ := deg g and R'_k := R̄_k/g. Then
    #{z ∈ Ω : g(z mod p) ≠ 0} ≤ Σ_{P: g(P)≠0} mult_P(R'₁,R'₂) ≤ dim F_p[X,Y]/(R'₁,R'₂) ≤ (d₁−δ)(d₂−δ).
(c) If R̄₁, R̄₂ are coprime, then |Ω| ≤ Σ_{P∈F_p², R̄₁(P)=R̄₂(P)=0} mult_P(R̄₁,R̄₂)
    ≤ dim F_p[X,Y]/(R̄₁,R̄₂) ≤ d₁d₂.

*Proof.*
1. Let S ⊆ Ω be finite. Choose N ≥ max(N_E^(1), N_E^(2)) with S ⊂ [0,N)². The cutoff is chosen after
   S, as in the paper.
2. G_k := p^{−s*_k} F^(k)_N lies in Z_(p)[X,Y] (Lemma 3) and reduces to R̄_k (Lemma 4). For s ∈ S,
   G_k(s) = p^{−s*_k} f_k(s) = 0 (1.3).
3. Apply Lemma 6 to J := (G₁, G₂), so J̄ = (R̄₁, R̄₂). This gives (a), and the first inequality in (c)
   with P ranging over S mod p; a point P with |S_P| > 0 has mult_P > 0, hence R̄₁(P) = R̄₂(P) = 0.
4. For (b): suppose g(P) ≠ 0. Then g = g(P)(1 − n) with n ∈ m_P, and g is a unit modulo m_P^k, with
   inverse g(P)^{−1}(1 + n + ⋯ + n^{k−1}). Hence (R̄₁,R̄₂) + m_P^k = (R'₁,R'₂) + m_P^k, so
   mult_P(R̄₁,R̄₂) = mult_P(R'₁,R'₂). Sum (a) over these P. Then use the general clause of Lemma 6(c) with
   I = (R'₁,R'₂), whose quotient is finite by Lemma 7. The first clause does not apply, because (R'₁,R'₂)
   is not the reduction of an ideal that vanishes on S. Lemma 7 also gives the last inequality.
5. (c) is (b) with g = 1. Since every finite S ⊆ Ω satisfies the bounds, Ω is finite and satisfies
   them. ∎

Theorem 8 bounds the number of zeros, not their positions. Its hypothesis is checkable from finitely
many b_ij (Lemma 4(c) and Lemma 1(d)). Part (b) needs coprimality only at the F_p-rational points:
for example, a common factor X² + 1 with p ≡ 3 mod 4 costs nothing.

## 5. Signed indices

**Lemma 9 (translation).** (a) Let (m₀,n₀) ∈ N². Then f'(m,n) := f(m+m₀, n+n₀) is a sequence of the
same form, with

   b'_ij = Σ_{k≤m₀, l≤n₀} C(m₀,k) C(n₀,l) q^{k+l} b_{i+k, j+l}.

For every N ≥ 1, F'_N ≡ F_N(X+m₀, Y+n₀) mod p^{Λ(N)}Z_(p)[X,Y]. Consequently s*' = s* and
R̄' = R̄(X+m₀, Y+n₀).
(b) In the matrix form with T₁, T₂ ∈ GL_d(Z), let (m₀,n₀) ∈ Z² and v' := T₁^{m₀}T₂^{n₀}v ∈ Z^d. Then
s*(v') = s*(v) and R̄_{v'}(X,Y) = R̄_v(X+m₀, Y+n₀).

*Proof.*
1. Vandermonde: C(X+m₀, a) = Σ_{k=0}^{a} C(m₀,k) C(X, a−k) in Q[X]. Both sides have degree ≤ a and
   agree on N. Substituting into f(m+m₀, n+n₀) = Σ q^{a+b} b_ab C(m+m₀,a) C(n+n₀,b), a finite sum,
   gives the formula for b'_ij. These are integers, or elements of Z_(p).
2. F_N(X+m₀,Y+n₀) = Σ C(m₀,k) C(n₀,l) q^{i+j+k+l} b_{i+k,j+l} C(X,i) C(Y,j) over i+k < N, j+l < N.
   F'_N is the same sum over i,j < N with all k ≤ m₀, l ≤ n₀. The difference consists of the terms
   with i,j < N and max(i+k, j+l) ≥ N. The content of such a term is
   ≥ ν(b_{i+k,j+l}) + (i+j+k+l)e − ν(i!) − ν(j!) ≥ val(i+k, j+l), since i! | (i+k)!, and this is
   ≥ Λ(i+j+k+l) ≥ Λ(N).
3. Translation by integers is an automorphism of Z_(p)[X,Y] and preserves content. So
   content(F_N(X+m₀,Y+n₀)) = s*_N (Lemma 3). For N ≥ N₀, Λ(N) > s*, so content(F'_N) = s*_N = s*;
   hence s*'_N = s* for all large N, and s*' = s*. Reducing p^{−s*}F'_N ≡ p^{−s*}F_N(X+m₀,Y+n₀) mod p
   gives R̄' = R̄(X+m₀,Y+n₀), since translation commutes with reduction.
4. (b): write (m₀,n₀) = (m₊,n₊) − (m₋,n₋) with all four ≥ 0. Put u := T₁^{−m₋}T₂^{−n₋}v ∈ Z^d,
   which is integral because T_l ∈ GL_d(Z). Then v = T^{(m₋,n₋)}u and v' = T^{(m₊,n₊)}u, and f_u is
   not identically 0, because f_v is a shift of f_u. By (a), R̄_v = R̄_u(X+m₋,Y+n₋) and
   R̄_{v'} = R̄_u(X+m₊,Y+n₊) = R̄_v(X+m₀,Y+n₀). ∎

**Theorem 10 (signed indices).** Let T₁, T₂ ∈ GL_d(Z) commute, T_l = I + qK_l, and
f_k(m,n) := l_k T₁^m T₂^n v for (m,n) ∈ Z², with neither f_k identically 0. Then Theorem 8(a)–(c)
hold for Ω^± := {(m,n) ∈ Z² : f₁ = f₂ = 0}, where the residue of (m,n) is (m mod p, n mod p).

*Proof.* Let S ⊆ Ω^± be finite and (m₀,n₀) its componentwise minimum. Then S − (m₀,n₀) ⊂ N² is a set
of common zeros for the seed v' := T^{(m₀,n₀)}v.

Neither f_k vanishes identically at the shifted seed. Suppose every b_ij(Tv) = l_kK₁^iK₂^jT₁v vanishes,
for one forward step T = T₁ (the case of T₂ is the same). Since T₁ = I + qK₁, this says
b_ij(v) = −q·b_{i+1,j}(v) for all i,j. So every b_ij(v) is divisible by every power of q, hence 0. A
backward step reduces to a forward one from the new seed. By Lemma 9(b), R̄_{v',k} = τ(R̄_k), where τ is the
F_p-algebra automorphism X ↦ X+m₀, Y ↦ Y+n₀. τ preserves degrees, gcds and colengths, and maps
m_{P} onto m_{P−(m₀,n₀)}, so mult_{P−(m₀,n₀)}(R̄_{v'}) = mult_P(R̄). Apply Theorem 8 to v'. ∎

This is the paper's shift of the seed to r = min S (main.tex lines 1882–1891), now in two indices.

## 6. Saturation

**Theorem 11 (saturated form).** Keep the hypotheses of Theorem 8. A *certificate* is a tuple
(A₁, A₂, τ, N₁) with A_k ∈ Z_(p)[X,Y], τ ≥ 0 and N₁ ≥ 1, such that
- (C1) A₁F^(1)_{N₁} + A₂F^(2)_{N₁} ∈ p^τ Z_(p)[X,Y];
- (C2) Λ(N₁) ≥ τ + 1 (precision).

Its reduction is h̄ := p^{−τ}(A₁F^(1)_{N₁} + A₂F^(2)_{N₁}) mod p.
(a) For every N ≥ N₁, H_N := p^{−τ}(A₁F^(1)_N + A₂F^(2)_N) lies in Z_(p)[X,Y], reduces to h̄, and
    vanishes at every common zero in [0,N)².
(b) For finitely many certificates with reductions h̄₁,…,h̄_r, Theorem 8(a) and the bound
    |Ω| ≤ Σ_P mult_P ≤ dim hold for the ideal (R̄₁, R̄₂, h̄₁, …, h̄_r) in place of (R̄₁, R̄₂).
(c) (Signed.) Under the hypotheses of Theorem 10, a certificate at the seed v yields, for every
    (m₀,n₀) ∈ Z², a certificate at T^{(m₀,n₀)}v with reduction h̄(X+m₀, Y+n₀). Hence (b) holds on Z².
(d) The stable reductions are certificate reductions: take A_k = 1, A_{3−k} = 0, τ = s*_k and
    N₁ ≥ N₀.

*Proof.*
1. F^(k)_N − F^(k)_{N₁} consists of the terms with max(i,j) ≥ N₁. Its content is ≥ Λ(N₁) ≥ τ + 1, by
   Lemma 1 and (C2). The A_k have content ≥ 0, so A₁F_N + A₂F_N ≡ A₁F_{N₁} + A₂F_{N₁} mod p^{τ+1}. This
   gives the first two claims of (a). At a common zero (m,n) ∈ [0,N)², F^(k)_N(m,n) = f_k(m,n) = 0, so
   p^τ H_N(m,n) = 0.
2. (b): for a finite S ⊆ Ω choose N beyond every N₁, every N_E and the coordinates of S. Then
   J := (G₁, G₂, H_{1,N}, …, H_{r,N}) vanishes on S and J̄ = (R̄₁, R̄₂, h̄₁, …). Apply Lemma 6.
3. (c), forward shift (m₀,n₀ ≥ 0): by Lemma 9 step 2, F'^(k)_N ≡ F^(k)_N(X+m₀,Y+n₀) mod p^{Λ(N)}.
   With A'_k := A_k(X+m₀,Y+n₀), the combination A'₁F'^(1)_N + A'₂F'^(2)_N is congruent to
   (A₁F^(1)_N + A₂F^(2)_N)(X+m₀,Y+n₀) modulo p^{Λ(N)}. For N ≥ N₁ with Λ(N) ≥ τ+1 the right side lies
   in p^τZ_(p)[X,Y] and reduces after division to h̄(X+m₀,Y+n₀). Backward shift (v = T^{(m₀,n₀)}v'):
   substitute X ↦ X−m₀, Y ↦ Y−n₀ in the same congruence and use A_k(X−m₀,Y−n₀). A general shift is
   a backward shift followed by a forward one. Then argue as in Theorem 10.
4. (d) is Lemma 4. ∎

(C2) is the finite form of a fact about the completion: a combination known modulo p^{Λ(N₁)} can be
divided by p^τ only while τ < Λ(N₁). Compare Lido–Mauri's Algorithm 1, in which the precision drops at
each saturation step (§9).

## 7. The determinant identity and what it forces

**Proposition 12 (frame identities).** Let d = 4. Let M ∈ M₄(Z_(p)) commute with T₁ and T₂, for
example M ∈ Z_(p)[K₁,K₂], and put Φ(w; M) := [w, K₁w, K₂w, Mw].
(a) Φ(T₁^mT₂^nv; M) = T₁^mT₂^n Φ(v; M), so det Φ(T₁^mT₂^nv; M) = det(T₁)^m det(T₂)^n det Φ(v; M).
    Since T_l = I + qK_l, det T_l ≡ 1 mod q is a p-adic unit. So ν(det Φ(w;M)) = ν(det Φ(v;M)) at every
    orbit point: (m,n) ∈ N², or Z² when T_l ∈ GL₄(Z).
(b) Field case. Let L be a quartic field and O an order with Z-basis ω₁,…,ω₄. Let ε_l = 1 + qη_l be
    units with η_l ∈ O, and β ∈ O. Let T_l, K_l be the matrices of multiplication by ε_l, η_l, let v
    be the coordinate column of β, and let M be multiplication by μ ∈ Z_(p)[η₁,η₂]. Then
    det Φ(β; M) = N(β)·det[1, η₁, η₂, μ], with coordinate columns. This uses the definition
    N_{L/Q}(β) = det(multiplication by β), a named standard fact. For μ = η₁η₂ the value is N(β)D₀ with
    D₀ := det[1, η₁, η₂, η₁η₂].
    Column operations (ε_l − 1 = qη_l and ε₁ε₂ − ε₁ − ε₂ + 1 = q²η₁η₂) give
    det[1, ε₁, ε₂, ε₁ε₂] = q⁴D₀. This is the analogue of Δ = 64Δ₀ = q³Δ₀ (q = 4) in main.tex line 1865.
(c) The frame form is Q_w(c) := det Φ(w; (c₁K₁+c₂K₂)²) = c₁²Q₂₀ + 2c₁c₂Q₁₁ + c₂²Q₀₂, with
    Q_ij := det Φ(w; K₁^iK₂^j), by multilinearity in the last column. In the field case
    Q_β(c) = N(β)·(c₁²D₂₀ + 2c₁c₂D₀ + c₂²D₀₂), where D₂₀ := det[1,η₁,η₂,η₁²] and
    D₀₂ := det[1,η₁,η₂,η₂²].

*Proof.* (a) T_l commutes with K₁, K₂ and M. (b) [β, η₁β, η₂β, μβ] = (multiplication by β)·[1, η₁, η₂, μ].
(c) is multilinearity. ∎

**Proposition 13 (Laplace and its consequences).** Let l₁, l₂ ∈ Z⁴ be marker rows whose reductions mod
p are independent. Choose l₃, l₄ among the standard rows so that U := [l₁; l₂; l₃; l₄] is invertible mod
p. For an orbit point w, let W(w; M) be the 2×4 matrix (l_k applied to the columns of Φ(w; M)), and put
s_M := ν(det Φ(v; M)).
(0) Some 2×2 minor of W(w; M) has valuation ≤ s_M.

For the mixed frame M = K₁K₂ write s := s_{K₁K₂}, and r_k(w) := (b^(k)₀₀, b^(k)₁₀, b^(k)₀₁, b^(k)₁₁)(w)
for the k-th row. Here b^(k)_ij(w) = l_kK₁^iK₂^jw are the coefficients of the shifted sequence
(m',n') ↦ f_k(m+m', n+n'). Then:
(1) Bounded content: ν(r₁(w)) + ν(r₂(w)) ≤ s. Hence s*_k ≤ 2e + ν(r_k(v)), s*₁ + s*₂ ≤ 4e + s, and
    deg R̄_k ≤ D(2e + s), uniformly along the orbit.
(2) Bounded dependence: for c ∈ Z_(p)² with p ∤ c, ν(c₁r₁(w) + c₂r₂(w)) ≤ s − min(ν r₁(w), ν r₂(w)) ≤ s.
(3) At a common zero, write x_k := l_kK₁w, y_k := l_kK₂w, z_k := l_kK₁K₂w. Then
    min(ν J, ν M₁₃, ν M₂₃) ≤ s, where J = x₁y₂ − x₂y₁ (the Jacobian), M₁₃ = x₁z₂ − x₂z₁ and
    M₂₃ = y₁z₂ − y₂z₁.

*Proof.*
1. det(UΦ) = det U·det Φ, and ν(det U) = 0. Expand det(UΦ) by Laplace along its rows 1 and 2: it is a
   signed sum of products (2×2 minor of W) × (complementary 2×2 minor from rows 3, 4), and the second
   factors lie in Z_(p). So ν(det Φ(w;M)) ≥ min ν(minors of W). Proposition 12(a) then gives (0).
2. (1): each minor is a sum of products of one entry of r₁ and one entry of r₂. Also
   s*_k ≤ min_{(i,j)∈{0,1}²} val(i,j) = min(ν b₀₀, e+ν b₁₀, e+ν b₀₁, 2e+ν b₁₁) ≤ 2e + ν(r_k), and s*_k is
   the same at every orbit point (Lemma 9). D is nondecreasing.
3. (2): if p ∤ c₁, the minors of (c₁r₁ + c₂r₂, r₂) are c₁ times the minors of (r₁, r₂), so
   ν(c₁r₁ + c₂r₂) + ν(r₂) ≤ s. The case p ∤ c₂ is symmetric.
4. (3): the minors containing the w-column vanish at a common zero. ∎

**Example 14 (the identity forces neither coprimality nor finiteness).** Let A = Z[x,y]/(x², y²) with
basis 1, x, y, xy. Let K₁, K₂ be multiplication by x, y, and T_l = I + qK_l; these are unipotent and lie
in GL₄(Z). Take v = 1, l₁ = the x-coordinate and l₂ = the xy-coordinate. Then

   T₁^mT₂^n·1 = (1+qx)^m(1+qy)^n = 1 + mq·x + nq·y + mnq²·xy   for all (m,n) ∈ Z²,

so f₁ = qm and f₂ = q²mn.
- det Φ(w; K₁K₂) = det[1, x, y, xy] = 1 at every orbit point, so s = 0 everywhere. The rows of W at w
  are (mq, 1, 0, 0) and (mnq², nq, mq, 1), and M₁₃ = 1.
- The common zeros {(0,n) : n ∈ Z} are infinite in number.
- R̄₁ = X and R̄₂ = XY share the factor X. The frame form is Q(c) = 2c₁c₂, isotropic at (1,0) and (0,1).

Every argument that uses only Proposition 13 (with any s) also applies here. So no such argument can
prove coprimality of R̄₁, R̄₂ or finiteness of the zero set. Checked in lemma_checks §7: frame
determinant 1 on [−6,6]², minimum minor valuation 0, and 13 zeros, all on m = 0.

**Example 15 (coprimality fails where the count is finite; D₀ = 0 is harmless).** Let
A = Z[x,y]/(x² − y², xy) with basis 1, x, y, x²; note x³ = x·y² = (xy)y = 0. With K_l, T_l and v = 1 as
above,

   T₁^mT₂^n·1 = 1 + mq·x + nq·y + (C(m,2)+C(n,2))q²·x²   on Z².

- D₀ = det[1, x, y, xy] = 0, so the mixed frame is degenerate (s = ∞).
- (c₁x + c₂y)² = (c₁² + c₂²)x², so Q(c) = c₁² + c₂². For p ≡ 3 mod 4 this form is anisotropic mod p
  and s_F = 0.
- Markers: l₁ = the x-coordinate, and l₂ = the x-coordinate + p^k·(y-coordinate) + (x²-coordinate),
  which are independent mod p. Then f₁ = qm and f₂ = qm + p^kq·n + q²(C(m,2) + C(n,2)).
- For k ≥ 1 and p odd: R̄₁ = R̄₂ = X, so Theorem 8 gives nothing.
- Common zeros on Z²: m = 0 and nq(p^k + (n−1)q/2) = 0. These are (0,0), and also (0, 1 − 2p^{k−e})
  when k ≥ e: at most 2.
- One saturation step H = p^{−τ}(F₂ − F₁), τ = min(e+k, 2e), gives h̄ = Y (k < e),
  (X² − X + Y² + Y)/2 (k = e), or (X² − X + Y² − Y)/2 (k > e). Then dim F_p[X,Y]/(X, h̄) = 1, 2, 2, equal
  to the number of zeros. Checked for p = 3, e = 2, k = 1, 2, 3 in lemma_checks §8.

This model satisfies the hypothesis of Theorem 18 (p = 3, e = 2 ≥ s_F + 2). So coprimality of the raw
reductions is not implied by (F4)–(F5). Also D₀ ≠ 0 is not necessary for a finite count.

The first referee's model satisfies both frame hypotheses (referee_checks §I): A = Z[x,y]/(x²−xy, y²+xy, (x,y)³),
basis 1, x, y, xy, with p = 3, e = 2, v = 1, l₁ = the x-coordinate and l₂ = x + 3^k·y + xy coordinates. It has
D₀ = 1 and Q = c₁²+2c₁c₂−c₂², which is anisotropic mod 3. Yet R̄₁ = R̄₂ = X, and the zeros are {(0,0)} or
{(0,0), (0, 1+2·3^{k−2})}.

Example 15 is the same phenomenon as Lido–Mauri's Example 3.2.3: reductions (x², xy) share the factor x, and
one saturation step repairs it.

**7.5. Where a coprimality proof breaks.**
1. By Lemma 4, R̄_k = Σ_{E_k} ū_ij X^(i)Y^(j), so coprimality is a condition on the minimizer sets E_k
   and the units ū_ij. A frame identity controls only finitely many entries l_kK₁^iK₂^jw. It does not
   control which (i,j) are minimizers.
2. *Minimizers outside {0,1}².* At p = 2, val(2,0) = 2e − 1 + ν(b₂₀) is smaller than
   val(1,1) = 2e + ν(b₁₁) whenever ν(b₂₀) ≤ ν(b₁₁): the pure second-order terms beat the mixed one.
   The mixed identity says nothing about b₂₀ or b₀₂.
3. *Linear minimizers with proportional linear parts.* At a common zero with E_k ⊆ {(1,0),(0,1)},
   R̄_k = (x_kX + y_kY)/p^{λ_k} mod p with λ_k := min(νx_k, νy_k). These are coprime iff
   ν(J) = λ₁ + λ₂. Proposition 13(3) gives only min(νJ, νM₁₃, νM₂₃) ≤ s, which is satisfied with
   J ≡ 0 mod p (J = 0 in Example 14; J = p^k in Example 15, where s = ∞).
   - Task B2 found natural field instances (necessity_examples.out.txt). In Z[2^{1/4}] at p = 7, 16 of 47 seeds
     x − yω have stable reductions with a common factor.
   - For all 16, the mixed frame determinant det Φ(v; K₁K₂) is a 7-adic unit (s = 0), and one saturation step
     gives colength 2.
   - ν₇(J) = 1 at the common zero for 14 of them, including the three detailed in that output (x − yω = −4−ω,
     −3−ω, 3−ω). It is 2 for −2−3ω and 2−3ω (task B2's verifier, independent_check.out.txt).
4. So the planned step of GENERALIZATIONS.md §3.2, coprimality forced by a determinant identity,
   cannot be completed. The obstruction is real (Example 14), not a gap in the argument. The right
   target is the saturated ideal (Theorem 11). In Lido–Mauri §4 the three raw reductions of a quintic Thue equation are all
   multiples of t₂ + 1, and one saturation step repairs this. Theorem 18 shows that, for two markers in
   two indices and under its hypotheses (F4)–(F5), the frame family guarantees that one explicit
   saturation step suffices.

## 8. What survives

**Lemma 16 (refinement).** For a ≥ 0, T_l^{p^a} = I + qp^a K_l^(a) with

   K_l^(a) = Σ_{r≥1} [C(p^a,r) q^{r−1}/p^a] K_l^r ∈ M_d(Z),

and these commute. Moreover K_l^(a) ≡ K_l mod p^{e−ι}. Each class (m₀,n₀) + p^aZ² is the Z²-orbit of
the system (T^(a), seed T₁^{m₀}T₂^{n₀}v) when T_l ∈ GL_d(Z); for (m₀,n₀) ∈ [0,p^a)² the class
(m₀,n₀) + p^aN² is the corresponding N²-orbit. The system has the form of §1 with q replaced by qp^a.

*Proof.* For 1 ≤ r ≤ p^a, ν(C(p^a,r)) = a − ν(r): write C(p^a,r) = (p^a/r)·C(p^a−1, r−1), and note
ν(p^a − k) = ν(k) for 1 ≤ k < p^a. So the bracket has valuation (r−1)e − ν(r) ≥ 0, since
ν(r) ≤ r − 1. For r = 2 this valuation is e − ι. For r ≥ 3 it is ≥ e, since (r−2)e ≥ ν(r): ν₃(3) = 1,
and r − 2 ≥ log₂ r for r ≥ 4. The remaining claims follow from T^{m₀+p^am'} = T^{m₀}(T^{p^a})^{m'}. ∎

**Proposition 17 (Jacobian criterion).** Let w be a common zero (taken as the seed), and let
x_k, y_k, λ := min ν(x₁,x₂,y₁,y₂) and J := x₁y₂ − x₂y₁ be as in Proposition 13(3). Put δ := ν(J) − λ.
Suppose λ ≤ e − 1 − ι and δ ≤ e − 1 − ι. Then w is the only common zero in w + N², and the only one in
Z² if T_l ∈ GL_d(Z). Here d is arbitrary.

*Proof.*
1. Relabel the markers and, if necessary, exchange the indices so that ν(x₁) = λ. Put
   ρ := x₂/x₁ ∈ Z_(p) and f'₂ := f₂ − ρf₁. On {f₁ = 0}, f'₂ = 0 iff f₂ = 0, so (f₁, f₂) and (f₁, f'₂)
   have the same common zeros.
2. f'₂ has b'₀₀ = 0, b'₁₀ = x₂ − ρx₁ = 0 and b'₀₁ = y₂ − ρy₁ = J/x₁, of valuation δ. So
   val'(0,1) = e + δ < 2e − ι = Λ(2) ≤ val'(i,j) for i+j ≥ 2. Hence E' = {(0,1)} and R̄(f'₂) = ūY.
3. For f₁: b₀₀ = 0, val(1,0) = e + λ < Λ(2), and val(0,1) ≥ e + λ. So
   R̄(f₁) = L̄ := (x₁X + y₁Y)/p^λ mod p = āX + b̄Y with ā ≠ 0.
4. dim F_p[X,Y]/(L̄, Y) = 1. Theorem 8(c) (or Theorem 10) applied to (f₁, f'₂) gives at most one
   common zero, and w is one. ∎

Provenance. Proposition 17 is the finite-polynomial form of Skolem's 1934 multivariate criterion, as Lombardo
states it (arXiv:2208.03830, Theorem 1.7): a simple common zero with nonsingular Jacobian mod p is isolated in
its residue class. Only the presentation is new.

Remark (refinement). K^(a) ≡ K mod p^{e−ι} (Lemma 16), so the data x_k, y_k at level a agree with
level 0 modulo p^{e−ι}. If the hypotheses of Proposition 17 hold at level 0, they hold at every level.
If ν(J) is too large at level 0, one computes it at level a, where e + a replaces e.

**Theorem 18 (frame criterion).** Hypotheses:
- (F1) d = 4; K₁, K₂ ∈ M₄(Z) commute; T_l = I + qK_l with p, q, e as in 1.1.
- (F2) l₁, l₂ ∈ Z⁴ are rows whose reductions mod p are linearly independent.
- (F3) v ∈ Z⁴ and f_k(m,n) := l_kT₁^mT₂^nv.
- (F4) (Frame hypothesis.) With Q_v(c) := det[v, K₁v, K₂v, (c₁K₁+c₂K₂)²v], there is an integer
  s_F ≥ 0 such that ν(Q_v(c)) ≤ s_F for every c ∈ Z² with p ∤ gcd(c₁,c₂). By Proposition 12(c),
  Q_v is a binary quadratic form. (F4) is a finite check modulo p^{s_F+1}. For odd p it holds with
  s_F = k whenever Q_v/p^k is anisotropic mod p, i.e. its discriminant is a non-residue.
- (F5) e ≥ s_F + 1 + ι' (that is, e ≥ s_F + 1 for p ≥ 5, and e ≥ s_F + 2 for p ∈ {2,3}).

Conclusions:
(a) For every common zero (m_w,n_w) ∈ N², at most 2 common zeros lie in (m_w,n_w) + N².
(b) If T₁, T₂ ∈ GL₄(Z), f₁ and f₂ have at most 2 common zeros in Z².
(c) More precisely, at each common zero w either w is the only common zero (Case J), or the count is
    at most deg h̄_w ≤ 2 for the explicit polynomial h̄_w of step 6 (Case D).

*Proof.*
1. *Invariance.* At an orbit point w = T^{(m,n)}v, Q_w(c) = det(T₁)^m det(T₂)^n Q_v(c), with p-adic
   unit factors (Proposition 12(a)). For c ∈ Z_(p)² with p ∤ c, Q_w(c) mod p^{s_F+1} depends only on
   c mod p^{s_F+1}, which is represented by an integer vector prime to p. So ν(Q_w(c)) ≤ s_F for all
   such c.
2. *Setup at a zero.* Fix a common zero w and use it as the seed. Let x_k, y_k be as in
   Proposition 13(3) and λ := min ν(x₁,x₂,y₁,y₂). Take any primitive c and M_c := (c₁K₁+c₂K₂)².
   In W(w; M_c) the w-column is 0, and each remaining minor is a sum of products containing one of
   x_k, y_k. Proposition 13(0) gives λ ≤ s_c := ν(Q_w(c)) ≤ s_F. In particular λ < ∞.
3. *Normalization.* Relabel the markers and/or exchange the indices (with K₁, K₂) so that ν(x₁) = λ.
   (F2), (F4) and (F5) are unchanged: exchanging K₁, K₂ replaces Q_w(c₁,c₂) by −Q_w(c₂,c₁). Put
   ρ := x₂/x₁ and κ := −y₁/x₁, both in Z_(p); l'₂ := l₂ − ρl₁; f'₂ := f₂ − ρf₁, with the same common
   zeros as before; b'_ij := l'₂K₁^iK₂^jw = b^(2)_ij − ρb^(1)_ij; and A_ij := b'_ij for i+j = 2. Then
   b'₀₀ = b'₁₀ = 0, y'₂ := b'₀₁ = J/x₁ and ν(y'₂) = ν(J) − λ. Finally c := (κ, 1), s_c := ν(Q_w(c)) ≤ s_F
   and L := x₁X + y₁Y = x₁(X − κY).
4. *Case J: ν(J) ≤ s_c + λ.* Then ν(J) − λ ≤ s_c ≤ e − 1 − ι' ≤ e − 1 − ι, and λ ≤ e − 1 − ι.
   Proposition 17 shows that w is the only common zero.
5. *Case D: ν(J) > s_c + λ.* Use the rows l₁, l'₂ (the matrix [l₁; l'₂; l₃; l₄] has the same
   determinant as U) and the columns w, K₁w, K₂w, M_cw. The rows of W are (0, x₁, y₁, z₁) and
   (0, 0, y'₂, γ), where z₁ := l₁M_cw and γ := l'₂M_cw = κ²A₂₀ + 2κA₁₁ + A₀₂. Proposition 13(0) gives
   min(ν(x₁y'₂), ν(x₁γ), ν(y₁γ − y'₂z₁)) ≤ s_c. The first equals ν(J) > s_c. If the second is ≤ s_c,
   then σ := ν(γ) ≤ s_c − λ. Otherwise ν(y₁γ − y'₂z₁) ≤ s_c, while ν(y'₂z₁) ≥ ν(J) − λ > s_c; so
   ν(y₁γ) ≤ s_c and σ ≤ s_c − ν(y₁) ≤ s_c − λ. In both cases **σ + λ ≤ s_c ≤ e − 1 − ι'**.
6. *The saturation certificate.* At the seed w, and for every N ≥ 3:
   - F^(1)_N = qL + E₁, where E₁ collects the terms with i+j ≥ 2 (b^(1)₀₀ = 0), so content(E₁) ≥ Λ(2) = 2e − ι.
   - F'_N = qy'₂Y + q²(A₂₀C(X,2) + A₁₁XY + A₀₂C(Y,2)) + R', where R' collects the terms with i+j ≥ 3,
     so content(R') ≥ Λ(3).
   - U := (q/x₁)[A₂₀(X + κY − 1)/2 + A₁₁Y] ∈ Z_(p)[X,Y] has content(U) ≥ e − λ − ι ≥ 0.
   - Expanding both sides gives the identity
     q²(A₂₀C(X,2) + A₁₁XY + A₀₂C(Y,2)) = qU·L − (q²/2)(A₂₀κ + A₀₂)·Y + (q²γ/2)·Y².
   - Hence H_N := F'_N − U·F^(1)_N = ξY + (q²γ/2)Y² + Rem_N, with ξ := qy'₂ − (q²/2)(A₂₀κ + A₀₂) and
     Rem_N := R' − U·E₁.
   - content(Rem_N) ≥ min(Λ(3), (e − λ − ι) + (2e − ι)) ≥ 3e − λ − ι − ι'. By 1.5, Λ(3) is 3e for
     p ≥ 5, 3e − 1 for p = 3, and 3e − 2 for p = 2.
   - ν(q²γ/2) = 2e − ι + σ ≤ 2e − ι + (e − 1 − ι' − λ) < 3e − λ − ι − ι' ≤ content(Rem_N).
   - So τ := content(H_N) = min(ν(ξ), 2e − ι + σ). If ν(ξ) ≥ content(Rem_N), both terms in the Y
     coefficient have valuation above 2e − ι + σ, so τ = 2e − ι + σ. Thus τ does not depend on N, and
     h̄ := p^{−τ}H_N mod p = αY + βY² ∈ F_p[Y], with (α,β) ≠ (0,0) and independent of N.
   - H_N = A₁F^(1)_N + A₂F^(2)_N with A₁ = −(ρ + U) and A₂ = 1. Together with any N₁ ≥ 3 such that
     Λ(N₁) ≥ τ + 1, this is a certificate in the sense of Theorem 11.
7. *Stable reduction of f₁.* val^(1)(1,0) = e + λ < 2e − ι = Λ(2), since λ ≤ e − 1 − ι' ≤ e − 1 − ι.
   Also val^(1)(0,1) ≥ e + λ and b^(1)₀₀ = 0. So R̄₁ = L̄ := (L/p^λ mod p) = āX + b̄Y with ā ≠ 0.
8. *Count.* F_p[X,Y]/(L̄) ≅ F_p[Y] via X ↦ −(b̄/ā)Y, and h̄ involves Y only. So
   dim F_p[X,Y]/(R̄₁, R̄₂, h̄) ≤ dim F_p[X,Y]/(L̄, h̄) = deg h̄ ≤ 2. Theorem 11(b) gives (a), and
   Theorem 11(c) gives (b). ∎

(F5) is used in steps 4, 6 and 7. For p = 3, step 6 has content(Rem_N) ≥ min(Λ(3), 3e − λ) = 3e − max(1, λ).
So e ≥ s_F + 1 already suffices at every zero with λ ≥ 1. The extra unit is spent only when λ = 0, where
Λ(3) = 3e − 1 is attained by a C(X,3) term.
- **The proof needs it.** The negative control in lemma_checks §9 shows that one less is not enough for this
  proof. At p = 3 with e = s_F + 1, the certificate has degree 3 and fails in 70 of 109 two-zero trials, all
  with λ = 0. The other 39 pass: 25 with λ = 0 and 14 with λ ≥ 1.
- **The statement needs it.** The second referee found an explicit system with p = 3, e = 1 and s_F = 0
  (adversarial_checks.out.txt):
  - A = Q × Q[x]/(x³), η₁ = (−1; x²), η₂ = (1; x), v = 1, q = 3, and Q(c) = c₂² + (c₁−c₂)²;
  - markers (1,−1,0,1) and (0,0,1,0), independent mod 3;
  - f₁ = (−2)^m 4^n − 1 + 3m + 9C(n,2) and f₂ = 3n;
  - exactly three common zeros, (0,0), (1,0) and (3,0).
  So the term [p ≤ 3] in (F5) is necessary at p = 3.
- **p ≥ 5 and p = 2.** For p ≥ 5 the threshold e = s_F + 1 passed every trial, and no three-zero example was
  found at p = 2 or p ≥ 5. Whether (F5) is sharp there is open.
- **p = 2 in particular.** Q mod 2 = (Q₂₀c₁ + Q₀₂c₂)² always has a primitive zero, so s_F ≥ 1 and (F5)
  forces e ≥ 3 (q ≥ 8). Theorem 18 never applies with q = 4.

Remark (integrality of T). T_l = I + qK_l has det ≡ 1 mod q, so it is invertible over Z_(p), and the proofs
run over Z_(p) (1.2). So the hypothesis T_l ∈ GL_d(Z) in Lemma 9(b) and Theorems 10, 11(c) and 18(b) can
be dropped, with f extended to Z² as a Z_(p)-valued function (second referee). In particular, Theorem 18(a)
gives at most 2 common zeros in N² itself.

**Corollary 19 (refinement form).** Let T_l ∈ GL₄(Z) and a ≥ 0. Let s_F^(a) be the frame depth (F4)
computed with K^(a) (Lemma 16). If e + a ≥ s_F^(a) + 1 + ι', then each class (m₀,n₀) + p^aZ² contains
at most 2 common zeros, so the orbit contains at most 2p^{2a}.

*Proof.* Apply Theorem 18(b) to each class (Lemma 16). The frame form at the seed T^{(m₀,n₀)}v differs
from that at v by the unit factor det(T₁)^{m₀}det(T₂)^{n₀} (Proposition 12(a)), so (F4) holds with the
same s_F^(a). ∎

**8.5. Two quartic orders of signature (2,1)** (lemma_checks §§6, 9, 10; exact arithmetic, the forms
computed modulo p^16).
- Z[ω], ω⁴ = ω + 1: discriminant −283, two real roots, Galois group S4 (sympy `galois_group`). The
  units are ω and ω − 1, of norms −1 and −1. For T_l = multiplication by ε_l^{h_l} (h_l the order mod p)
  and seed 1, the frame form is anisotropic mod p with s_F = 0 and e = 1 at p = 5 and p = 11. The same
  holds for exponent lcm(p⁴−1, p³−1) at p = 5, 11, 19, 23, 31 (for a = 0 and a = 1).
  - It is not anisotropic mod p at p = 3, 7, 13, 17, 29.
  - At p = 13 and 17 the form is degenerate mod p but (F4) still holds, with s_F^(a) = 1, 1, 1 (p = 13) and
    1, 3, 2 (p = 17) at levels a = 0, 1, 2. So Corollary 19 applies at p = 13 with a = 1 and at p = 17 with
    a = 2 (first referee, referee_checks §K2).
  - At p = 3, 7, 29 the form has a p-adic zero, and (F4) fails (checked for a ≤ 2).

  At p = 5 the orders are h = (312, 104), and Theorem 18(b)
  applies with no refinement. So **every coset β·⟨ω^{312}, (ω−1)^{104}⟩ of a unit β contains at most two
  elements x − yω, i.e. at most two solutions of x⁴ − xy³ − y⁴ = ±1** (the markers are the ω²- and
  ω³-coordinates, and Q_β = N(β)·Q₁ with N(β) = ±1). This is not a count of all solutions: the number of
  such cosets (the unit index) is not computed here. The signed planted tests at p = 5 (a = 0, 1) and
  p = 11 (a = 0) agree.
- Z[2^{1/4}], ω⁴ = 2: discriminant −2048, Galois group D4, units 1+ω and 1+ω² (task B2's example).
  Here ε₂ = 1+ω² ∈ Z[√2], so every η₂ lies in Q(√2) and η₂² ∈ span_Q(1, η₂). Hence
  Q(0,1) = D₀₂ = det[1, η₁, η₂, η₂²] = 0 identically. (F4) fails at every prime and every refinement
  level. With other generators of a finite-index subgroup, a power of 1+ω² still occurs, and its
  direction makes ν(Q) grow with the refinement level (p-adic limit, not proved here). Theorem 18 does
  not reach this order. Proposition 17 and Theorem 11 still apply zero by zero.
- *Which fields Theorem 18 can reach.* A quartic field of signature (2,1) has Galois group D4 or S4. C4, V4 and
  A4 fields are totally real or totally complex, since complex conjugation fixes no root or all four.
  - A D4 field contains a quadratic subfield. It is real, since an imaginary one would make every embedding
    complex, so it contributes a unit of infinite order.
  - The frame form vanishes in that unit's direction when the unit is a generator, as for Z[2^{1/4}].
  - So Theorem 18 is in practice a tool for S4 fields. For D4 fields, failure with other generators is expected
    but not proved.

## 9. Comparison with Lido–Mauri

Source: G. M. Lido and L. Mauri, *A multivariate Strassmann theorem*, arXiv:2605.03683v1 (5 May 2026),
HTML at https://arxiv.org/html/2605.03683, read on 2026-10-15. The statements below are my paraphrases;
the one direct quotation is marked.

What they prove (paraphrase):
- **Theorem 1.6.** K is complete for a non-archimedean valuation, not necessarily discrete, with
  valuation ring R and residue field k. Let f₁,…,f_r ∈ R⟨x₁,…,x_n⟩ have reductions f̄_i. If
  Ā = k[x]/(f̄₁,…,f̄_r) is finite-dimensional, the common zeros in R^n are finite in number, and their
  number is at most the sum, over the maximal ideals m of Ā with residue field k, of dim Ā_m, which is
  ≤ dim Ā. They note that this generalizes Edixhoven–Lido, Geometric quadratic Chabauty (2023),
  Theorem 4.12.
- **Theorem 2.5.** The saturation is I^ec = {f : rf ∈ I for some nonzero r ∈ R} (their (1.8)). Five
  conditions are equivalent, among them: finiteness of the common zeros in the integers of the completed
  algebraic closure, and finite colength of the reduction of I^ec.
- **Section 3** (discrete valuation). The generators are known only modulo π^N. The saturation chain is
  I_l = (I : π^l). The next ideal I_{l+1} is computed from syzygies of the reductions (Lemma 3.1.5,
  Corollary 3.1.6, Algorithm 1), with the precision dropping as l grows. Proposition 3.2.1 is a
  sufficient dimension criterion for I_l = I^ec, and reduced Gröbner bases appear in
  Propositions 3.4.3–3.4.4 and Algorithm 2. They state (§3.2) that they know no general method to
  decide from approximate information whether the saturation chain has already stabilized.
- **Section 4.** The quintic Thue equation X⁵ − X⁴Y + X³Y² − X²Y³ − Y⁵ = 1 over Z[α], with unit rank 2,
  at p = 5.
  - In one residue class the three marker series F₂, F₃, F₄ all reduce to multiples of t₂ + 1, and they
    note that the quotient "is not finite, we cannot apply Theorem 1.6."
  - One saturation step gives a reduction generated by t₂ + 1 and t₁² − t₁, hence at most 2 zeros in that
    class.
  - They conclude that (−1,−1), (0,−1), (1,0) and (5,4) are all the solutions, with "at most 4" in total. That
    is the PDF; the HTML shows "44". The per-class bounds of the other residue classes are not displayed.
  - The series F₂ and F₃ printed in §4 contradict the paper's own units. The stated reduction and the bound 2
    are nevertheless correct by an independent recomputation
    ([verify_literature_checks.out.txt](../../targets/rank2_signatures_2026-10/verify_literature_checks.out.txt)).

What our finite lemma reproduces:
- Theorem 8 is the case n = r = 2, R = Z_p of their Theorem 1.6, restricted to integer zeros and to the
  series F_∞ = lim F_N (a Mahler-type series). For comparison only, and using the completion:
  p^{−s*}F_∞ ∈ Z_p⟨X,Y⟩ reduces to our R̄ by Lemma 4, since F_∞ ≡ F_N mod p^{s*+1} for N ≥ N₀.
- By Lemma 7, our coprimality hypothesis is equivalent to their finite dimensionality of Ā when
  n = r = 2.
- Our local bound (Lemma 6, Theorem 8(a),(c)) is their sum over k-rational maximal ideals, for
  integer zeros.
- Theorem 11 is the sufficient direction of their Theorem 2.5 at finite precision. Condition (C2) is
  the analogue of their precision accounting.

What they have and we lack:
- general n and r, and non-discrete valuations;
- zeros in R^n, not only integer points;
- the characterization (Theorem 2.5), i.e. necessity;
- the algorithms for saturation from approximate data and the dimension criterion;
- a fully solved Thue example.

What ours adds, as the project's elementary standard requires:
- Only finite integer polynomials, with the cutoff chosen after the finite set of zeros. No completion,
  Tate algebra, Weierstrass preparation or decomposition of complete rings is used. The named algebraic
  facts are the UFD property of F_p[X,Y], the PID structure of Z_(p), and the Chinese remainder theorem.
- The exact content and the explicit stable reduction (Lemmas 3–4), with explicit indices N_E ≤ N₀.
- The gcd form of Theorem 8(b), which does not assume dim Ā < ∞; their Theorem 1.6 as stated does.
- The translation lemma for signed orbits.
- The frame identities and Theorem 18: an a priori criterion, checkable on the seed, guaranteeing that
  one explicit saturation step gives colength ≤ 2. Nothing of this kind appears in the text I read.

## 10. Checks run (lemma_checks.out.txt)

1. Lemmas 3–4: 1784 truncations, 0 mismatches in the content, in k_N or in the stable formula. In 85
   tables N_E < N₀.
2. Corollary 5 (one index, q = 4): 0 violations of k ≤ ν(N!)+s+3 or of deg R̄ = D* ≤ D(s*) ≤ s+2.
3. Lemma 2: 0 failures.
4. Lemma 7: the Hilbert identity holds exactly. In 116 random coprime pairs over F_p, 0 have
   dim > d₁d₂ (88 have equality).
5. Lemma 6: in 114 residue classes, the generalized eigenspace has dimension |S_P|.
6. Lemma 9: 24 shifts, including negative ones, with genuine units of Z[ω], ω⁴ = ω + 1, at p = 5;
   0 mismatches.
7–8. Examples 14 and 15: as stated in §7.
9. Theorem 18:
   - two planted zeros at the threshold: 138 trials (p odd) and 63 (p = 2), all Case D with exactly
     2 zeros, deg h̄ = 2, and every proof step verified;
   - one planted zero: 65 Case J and 11 Case D, all consistent;
   - signed tests with genuine units: 18 trials, all consistent;
   - negative control: as stated after the proof.
10. The frame forms of §8.5.

## 11. Status

**PROVED**, as a written proof using finite integer polynomials only. In catalogue terms it is WRITTEN: refereed
twice in-session, with no outside audit. Hypotheses: p prime, q = p^e (e ≥ 1, and e ≥ 2 for p = 2),
b_ij ∈ Z or Z_(p), and neither sequence identically 0.
- Lemma 1 (growth/tail), Lemma 2 (Stirling), Lemma 3 (exact content, k_N = 2ν(N!) + s*_N), Lemma 4
  (explicit stable reduction, stable for N ≥ N_E, with N_E ≤ N₀ = D(s*)+1), Corollary 5 (degree
  ≤ D(s*); the one-index specialization reproduces the paper's bounds).
- Lemma 6 (colength, local and global) and Lemma 7 (affine Bézout).
- **Theorem 8**: if R̄₁, R̄₂ are coprime, the zeros in N² number at most
  Σ_P mult_P ≤ dim ≤ d₁d₂, with local and gcd forms.
- Lemma 9 and **Theorem 10** (Z² when T_l ∈ GL_d(Z)).
- **Theorem 11** (saturation certificates with precision Λ(N₁) ≥ τ+1; signed version).
- Propositions 12–13 (frame identities; consequences (1)–(3)).
- Example 14: the identity, even with s = 0 everywhere, forces neither coprimality nor finiteness.
- Example 15: coprimality fails under the hypothesis of Theorem 18, and D₀ = 0 is compatible with a
  finite count.
- Lemma 16, Proposition 17 (Jacobian criterion; Skolem 1934 in finite form), **Theorem 18** (frame criterion:
  at most 2 zeros in the orbit of v under (F4) and e ≥ s_F + 1 + [p ≤ 3]; the term [p ≤ 3] is necessary at
  p = 3), Corollary 19.
- The two computed quartic examples of §8.5, as computations.

**OPEN.**
1. Whether the coprimality hypothesis of Theorem 8 can be forced from the structure of a Thue field.
   It cannot be forced by the determinant identity (Example 14), and it is not the right target
   (Example 15).
2. Whether (F4) holds at some prime, after some refinement, for every quartic order of signature (2,1)
   and its unit group.
   - It fails structurally when a unit direction lies in a quadratic subfield: D4 fields, e.g. Z[2^{1/4}].
   - For S4 orders it is a p-adic anisotropy condition, and no a priori reason for it to hold is known.
   - For x⁴ − x − 1 the frame form is anisotropic mod p (s_F = 0) at 5 of the 10 primes tested. With
     refinement, the hypotheses of Theorem 18 or Corollary 19 hold at 7 of the 10 (p = 13 and 17 added, §8.5).
3. A finite-polynomial converse to Theorem 11 (Lido–Mauri Theorem 2.5): that a certificate exists
   whenever the p-adic zero set is finite.
4. Frames and the count for more markers or more indices: d ≥ 5 (for example quintic Thue equations
   with 3 markers in 2 indices, as in Lido–Mauri §4) and r ≥ 3 indices. Theorem 18 is for d = 4 only.
5. Outside this lemma: coverage of the solutions by finitely many orbits (the descent of
   ghost.unit.multi_unit), and the unit index needed to turn §8.5 into a count of solutions.
6. Whether (F5) is sharp at p = 2 and at p ≥ 5. No orbit with three zeros at e = s_F + 1 was found.
7. Iterated saturation. One step is not always enough: 59 of the 167 non-coprime synthetic cases of
   task B2 keep infinite colength after one step. No bound on the number of steps is known (compare
   Lido–Mauri's open Question 3.2.4).

## Sources
- research/shared/prior/elementary_diophantine_paper/main.tex: section `sec:returns` (lines 1339–1448:
  Theorem `three-return`, Corollary `matrix-return`), subsection "A sharper count of all marked returns"
  (lines 1849–1910: Δ = 64Δ₀, det[β, Kβ, K²β] = mΔ₀, signed shift).
- research/proofs/GENERALIZATIONS.md §3.2; research/proofs/methods/lifting.binomial.md.
- Lido, Mauri, arXiv:2605.03683v1, https://arxiv.org/html/2605.03683 (fetched 2026-10-15). Their
  reference list includes Strassmann (1928), Skolem (1934), Stroeker–Tzanakis (1988), Lombardo (2022)
  and Edixhoven–Lido (2023).
- Computations: lemma_checks.py / lemma_checks.out.txt (this folder).
