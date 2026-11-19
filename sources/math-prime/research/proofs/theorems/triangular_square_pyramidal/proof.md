# Triangular = square-pyramidal (8604.a1) by strategy R1: two-unit ghosts and finite-polynomial exclusion

Written 2026-10-15 by Claude (Opus 5.5) as the second run of strategy R1 of
[GENERALIZATIONS.md](../../GENERALIZATIONS.md) §3.2, after [triangular_tetrahedral](../triangular_tetrahedral/proof.md)
(1611.a1). The machinery is the same; the data are new. Scripts: [scripts/](scripts/). Each script has its output
beside it as `*.out.txt`. Every check is exact, and floating point only seeds searches. Hypothesis (H) of §1.2, which
this proof first took from the rank, is proved in [hypothesis_H.md](hypothesis_H.md) (added 2026-10-15, reviewed in
round 3 of [review.md](review.md)).

## 0. Statements

**0.1 The curve.** E: y² = x³ − 9x + 81. This is LMFDB 8604.a1, whose model is already monic, so X = x and Y = y
throughout. Its discriminant is −2787696 = −2⁴·3⁶·239. LMFDB gives rank 2, an mwrank value. The proof does not use
it (see 0.4).
- Write T(n) = n(n+1)/2 and P(m) = m(m+1)(2m+1)/6 = 1² + 2² + ⋯ + m².
- T(n) = P(m) holds exactly when (x, y) = (6m+3, 18n+9) is a point of E, by the identity
  (18n+9)² − ((6m+3)³ − 9(6m+3) + 81) = 648·(T(n) − P(m)).

**0.2 Class forms and field forms.** The 2-descent of E (§1) attaches to each non-trivial class a binary quartic
*class form* F_cl. The Thue equations actually solved are those of a monic *field form* G, the norm form
N(n − mβ) of a generator β whose order O satisfies O ∩ (Q + Qβ) = Z + Zβ (§2.2). For two classes G = F_cl. For class
(1,0), (1+α)/2 is integral for a root α of F_cl(1,t), so F_cl itself does not have this property (§2.3).

| class | μ | F_cl(m,n) | G(m,n) | β | relation |
|---|---|---|---|---|---|
| (1,0) | 3 − θ | 9m⁴ − 72m³n − 18m²n² + n⁴ | 4m⁴ − 5m³n − 3m²n² − 2mn³ + n⁴ | (1+α)/2 | F_cl(m,n) = G(2m, m+n) |
| (0,1) | −3 − θ | 9m⁴ − 72m³n + 18m²n² + n⁴ | the same | α | G = F_cl |
| (1,1) | −θ | 36m⁴ − 72m³n + n⁴ | the same | α | G = F_cl |

**Theorem A′ (the twelve Thue equations of the field forms).** For each class and each c ∈ {1, 4, 9, 36}, the
integer solutions of G(m,n) = c are exactly ±(m,n) for the following pairs (m,n):

| class | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|
| (1,0) | (0,1) | (1,0) | (1,−1), (2,1) | none |
| (0,1) | (0,1) | none | (1,0) | (1,3) |
| (1,1) | (0,1), (2,1) | (1,4) | none | (1,0) |

**Theorem A (the class forms).** For each class and each c ∈ {1, 4, 9, 36}, the integer solutions of F_cl(m,n) = c
are exactly ±(m,n) for:

| class | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|
| (1,0) | (0,1) | none | (1,0) | none |
| (0,1) | (0,1) | none | (1,0) | (1,3) |
| (1,1) | (0,1), (2,1) | (1,4) | none | (1,0) |

**Theorem B (integral points).** The following hold. (Hypothesis (H) of 1.2, on which this rests, is Theorem H of
[hypothesis_H.md](hypothesis_H.md).)
- the integral points of E are the 24 points with x ∈ {−5, −3, 0, 3, 7, 9, 24, 33, 39, 513, 1099, 5112}, with both
  signs of y;
- T(n) = P(m) with n, m ≥ 0 holds only for the common values 0, 1, 55, 91 and 208335, that is
  (m, n) ∈ {(0, 0), (1, 1), (5, 10), (6, 13), (85, 645)}.

**0.3 What is used.**
- **Theorems A′ and A** use:
  - exact computation with integers and rationals;
  - real numbers only through enclosures in rational intervals, obtained from exact sign changes, with the number of
    real roots fixed by Sturm's theorem;
  - elementary real-variable facts in Lemma 3.2: positive definiteness of the size form E of §3.1 (not the curve),
    that a bounded polygon is the convex hull of its vertices, and concavity of the logarithm of a positive affine
    function (or the root-free AM–GM of Remark 3.2(ii)); Remark 3.2(i) also uses e^y ≥ 1 + y, but Theorem A′ does not
    need that remark (§4);
  - Proposition 17 of the finite-polynomial lifting lemma
    [lifting.binomial.multi](../../methods/lifting.binomial.multi.md);
  - congruences.

  They use no completion, no p-adic power series, no ideals, no class numbers, no unit theorem and no linear forms in
  logarithms.
- **Theorem B** adds the 2-descent of §1 and [hypothesis_H.md](hypothesis_H.md).
  - **(H)**, the list of four classes, is Theorem H there. It is proved by:
    - a Euclidean division in Z[λ] (λ = θ/3), certified by 116 exact boxes (Z[λ] is norm-Euclidean);
    - the divisibility theory that follows from that division;
    - the units ±εᵏ, by a one-unit descent;
    - a congruence modulo 4.

    The rank of E(Q) and the Mordell–Weil theorem (N5) are no longer used.
  - **N1 and N2** of the descent script (the index–discriminant formula; integral elements form a ring, and the
    characteristic-polynomial criterion) are replaced by hypothesis_H.md (3.6): Z[λ] is integrally closed, and every
    element of it is integral. So O_K = Z[λ], as §1.3 uses (hypothesis_H.md §6).

  The rest of §1 (conics, parametrization, gcd bounds) is elementary.

**0.4 Rank.** LMFDB lists rank 2 for 8604.a1 with generators (3, 9) and (−3, 9). The proof does not use the rank:
rank ≤ 2 was the original route to (H), which hypothesis_H.md proves directly.

## 1. The reduction to Thue equations (the 2-descent; prior work)
Source: [weierstrass_2descent.py](../../targets/rank2_signatures_2026-10/weierstrass_2descent.py) →
[.out.txt](../../targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt), section `CURVE 8604.a1`. Every item
used below is re-derived by exact arithmetic in [03_reconstruction.out.txt](scripts/03_reconstruction.out.txt).

**1.1 The cubic field and the descent map.**
- Put K = Q(θ) with θ³ − 9θ + 81 = 0. Its signature is (1,1). O_K = Z + Z·θ/3 + Z·θ²/9 = Z[λ] with λ = θ/3 and
  λ³ − λ + 3 = 0, which is the cubic of the model 3Y² = X³ − X + 3 used by Uchiyama (§5), and [O_K : Z[θ]] = 27.
  That Z[λ] is the whole ring of integers is [hypothesis_H.md](hypothesis_H.md) (3.6); the descent script's N1
  (disc Z[λ] = −239 is squarefree) is not used.
- For a point P ≠ O put φ(P) = X(P) − θ in K*/K*². φ is a homomorphism whose kernel contains 2E(Q), by the identity
  (x₁−θ)(x₂−θ)(x₃−θ) = (λ'θ+ν)² for three collinear points on Y = λ'X + ν, together with the vertical-line case
  (checked again in 03 on a sample).
- E(Q)[2] = 0, because the cubic is irreducible (rational-root test).

**1.2 Hypothesis (H).** Every integral point P of E has φ(P) ∈ {1, 3−θ, −3−θ, −θ}.
- These four classes are φ of O, G₁ = (3, 9), G₂ = (−3, 9) and (0, 9) = −G₁ − G₂, for the LMFDB generators G₁, G₂.
  They are pairwise distinct: 3−θ, −3−θ, −θ and their pairwise products are non-squares, by a Legendre-symbol
  witness at a root of the cubic modulo 7 or 13 (checked exactly). Such a witness is constant on a class modulo
  squares because a square root in K of an element of Z[λ] lies in Z[λ] ([hypothesis_H.md](hypothesis_H.md) (5.8) and
  (3.6)). Theorem B does not need the distinctness.
- (H) would follow from rank E(Q) ≤ 2, by N5, since then |E(Q)/2E(Q)| ≤ 4. This was the original route.
- **(H) is proved for integral points in [hypothesis_H.md](hypothesis_H.md), Theorem H**, without the rank:
  X − θ = μ₀S² with S ∈ Z[λ] and μ₀ ∈ {1, π₃ₐπ₃꜀, π₃ₐπ₃ᵦ, π₃ᵦπ₃꜀} (π₃ₐ = λ, π₃ᵦ = λ − 1, π₃꜀ = λ + 1), and each μ₀ is
  one of the four classes times a square.
  - For the points with 3 | x, Uchiyama obtains the analogue of (H) from class number 1 and a fundamental unit of K
    (§6.3). hypothesis_H.md replaces both: its Sections 2–3 give the divisibility theory that class number 1
    provided, and its Section 4 proves that ε is a fundamental unit.

**1.3 From a class to a Thue equation.** Let μ be the class and X − θ = μξ².
- N(μ) = 81 for the three non-trivial classes. Ξ_μ = {ξ : μξ² ∈ O_K}, with O_K = Z[λ], is a lattice with
  O_K ⊆ Ξ_μ ⊆ (1/9)O_K, because Z[λ] is integrally closed ([hypothesis_H.md](hypothesis_H.md) (3.6) and §6); an
  enumeration of (1/9)O_K/O_K shows that Ξ_μ contains O_K with index 3.
- In a basis of Ξ_μ, ξ has coordinates (u,v,w) with A₂ = 0 (a smooth conic, Gram determinant −9), A₁ = −1 and A₀ = X.
- A common factor k of (u,v,w) satisfies k² | 3, because θ-coordinates of elements of O_K lie in (1/3)Z; so k = 1.
- The conic has the rational point ξ = 1. So (u,v,w) = ±𝐗(m,n)/g for coprime (m,n), by Lemma 1.3′ below. Here 𝐗(m,n)
  is the parametrization and g = gcd 𝐗(m,n) divides 36; that divisibility follows from the polynomial identities
  that the descent script checks.
- Then A₁ = −1 becomes F_cl(m,n) = g², with g ∈ {1, 2, 3, 6}.
- *g is forced.* In each non-trivial class u² − F_cl is 12 times an integral form, so u² ≡ F_cl (mod 12), and g is
  gcd(u, 6). Since r² ≡ 1, 4, 3, 0 (mod 6) exactly when gcd(r, 6) = 1, 2, 3, 6, every coprime solution of
  F_cl(m,n) = c ∈ {1, 4, 9, 36} has g = √c and gives an integral point of the class.
- X = X₄(m,n)/g² and Y = ±Y₆(m,n)/g³; the polynomials are in [core.py](scripts/core.py), `descent_maps`, and
  Y₆² = X₄³ − 9X₄F_cl² + 81F_cl³ identically.
- *The trivial class* (μ = 1) gives H = (4/3)F₀ = g² with F₀ = m³n − mn³ + 3n⁴ = n(m³ − mn² + 3n³) and g | 2. So
  F₀ = 3 and g = 2. Then n | 3, and for each n the rational-root test leaves the solutions ±(1,1), ±(1,−1), ±(5,−3)
  and ±(0,1). The last has g = 1 and gives only the rational point (1/4, ±71/8). The trivial class gives
  X ∈ {−5, 7, 1099} (§5).

So every integral point outside the trivial class gives a coprime solution of F_cl(m,n) = c ∈ {1, 4, 9, 36} in one
of the three classes. Theorem A lists all of them.

**1.3′ Lemma (parametrization of a conic).** Let q be a ternary quadratic form with integer coefficients and nonzero
Gram determinant, with polar form B(P, P′) = q(P + P′) − q(P) − q(P′). Let P₀ ∈ Z³ be primitive with q(P₀) = 0, and let
(P₀, e, f) be a basis of Z³. For d = me + nf put 𝐗(m,n) = q(d)P₀ − B(P₀,d)d. Then q(𝐗(m,n)) = 0. Every primitive
integer solution P of q(P) = 0 equals ±𝐗(m,n)/g for some coprime (m,n), where g is the gcd of the coordinates of
𝐗(m,n).

*Proof* (that of [triangular_tetrahedral/proof.md](../triangular_tetrahedral/proof.md), Lemma 1.3′).
1. *q(𝐗) = 0.* q(𝐗) = q(d)²q(P₀) − q(d)B(P₀,d)² + B(P₀,d)²q(d) = 0.
2. *No isotropic d ≠ 0 in span(e, f) has B(P₀, d) = 0.* Otherwise q would vanish on span(P₀, d), a plane. In a basis
   (P₀, d, w) the Gram matrix would then be [[0,0,x],[0,0,y],[x,y,z]], whose determinant is 0, a contradiction.
3. *Write the solution in the basis.* Write P = aP₀ + d with a ∈ Z and d = me + nf.
4. *If d = 0,* then P = ±P₀. The linear form (m,n) ↦ B(P₀, me + nf) is not identically 0 (otherwise B(P₀, ·) would
   vanish on all of Q³, since B(P₀, P₀) = 2q(P₀) = 0, against the nonzero determinant). So it has a primitive zero
   (m₀, n₀); by step 2, q(m₀e + n₀f) ≠ 0, so 𝐗(m₀, n₀) is a nonzero multiple of P₀, and P = ±𝐗(m₀, n₀)/g.
5. *If d ≠ 0,* then 0 = q(P) = aB(P₀,d) + q(d). By step 2, B(P₀,d) ≠ 0, so B(P₀,d)·P = −𝐗(m,n). If h = gcd(m,n), then
   𝐗(m,n) = h²𝐗(m/h, n/h), so P is proportional to 𝐗(m′,n′) with (m′,n′) coprime; 𝐗(m′,n′) ≠ 0, because its
   (e, f)-coordinates are −B(P₀,d′)(m′,n′) ≠ 0. P is primitive, so P = ±𝐗(m′,n′)/g. ∎

The descent script applies the lemma to the conic C_μ = (A₂ = 0), with P₀ = ξ = 1 and the bases it prints.

## 2. The Thue fields
Script: [00_fields.py](scripts/00_fields.py) → [.out.txt](scripts/00_fields.out.txt).

**2.1 The fields.** Write G(m,n) = a₀m⁴ + a₁m³n + a₂m²n² + a₃mn³ + n⁴, so that g(t) := G(1,t) =
t⁴ + a₃t³ + a₂t² + a₁t + a₀. For each class, g is irreducible with exactly two real roots (Sturm) and Galois
group S4. Write β for a root and L = Q(β). Then G(m,n) = m⁴g(n/m) = N(n − mβ). The order O used below has the
integral basis ω₀..ω₃ of [core.py](scripts/core.py) (`CLASSES[key]['BASIS']`):

| class | g(t) | ω₂ | ω₃ | [O : Z[β]] | disc O |
|---|---|---|---|---|---|
| (1,0) | t⁴ − 2t³ − 3t² − 5t + 4 | (1 + 2β + β²)/3 | (1 + β³)/3 | 9 | −3²·239 |
| (0,1) | t⁴ + 18t² − 72t + 9 | (3 + β²)/6 | (9 + 9β + β² + β³)/12 | 72 | −2⁶·3²·239 |
| (1,1) | t⁴ − 72t + 36 | β²/6 | (6β + β³)/12 | 72 | −2⁶·3²·239 |

In each case ω₀ = 1 and ω₁ = β. O is the maximal order (p-maximality at 3, and at 2 for classes (0,1) and (1,1)),
but maximality is never used. What
is used is that O is a ring containing Z[β] and the units below. The three fields are pairwise non-isomorphic: g
factors mod 13 (classes (1,0), (0,1)) and mod 7 (the other two pairs) with different degree patterns. Classes (0,1)
and (1,1) have the same discriminant.

**2.2 Markers.** Write x = (x₀, x₁, x₂, x₃) for the coordinates in the basis ω.
- The β²- and β³-coordinates of x are (x₂, x₃)·D for an invertible 2×2 matrix D: diag(1/3, 1/3) for class (1,0),
  [[1/6, 0], [1/12, 1/12]] for class (0,1) and diag(1/6, 1/12) for class (1,1).
- Since 1, β, β², β³ are linearly independent (g is irreducible), x ∈ Q + Qβ if and only if x₂ = x₃ = 0, and then
  x = x₀ + x₁β ∈ Z + Zβ.
- n − mβ has coordinates (n, −m, 0, 0).
- So **the solutions of G(m,n) = c are the elements z ∈ O with N(z) = c and x₂(z) = x₃(z) = 0.**

**2.3 From the class forms to the field forms.**
- *Why class (1,0) needs a substitution.* Let α be a root of f_cl(t) = F_cl(1,t) = t⁴ − 18t² − 72t + 9. Then
  β := (1+α)/2 is a root of g(t) = t⁴ − 2t³ − 3t² − 5t + 4, so it is integral: 16·g((1+t)/2) = f_cl(t). In the
  maximal order of Q(α), the element (1+α)/2 lies in Q + Qα but not in Z + Zα, so the marker lemma 2.2 fails for
  F_cl itself. With β as generator it holds (2.1).
- *The substitution.* N(n − mα) = N((n+m) − 2mβ), that is F_cl(m,n) = G(2m, m+n), an identity of forms
  ([00_fields.out.txt](scripts/00_fields.out.txt), §4).
- *Lemma.* For every integer c, (m,n) ↦ (2m, m+n) is a bijection from the integer solutions of F_cl(m,n) = c onto
  the integer solutions (m', n') of G(m',n') = c with m' even. Its inverse is (m', n') ↦ (m'/2, n' − m'/2).
  *Proof.* (m,n) ↦ (2m, m+n) maps Z² bijectively onto {(m', n') : m' even}, with that inverse, and the identity of
  forms carries one equation into the other. ∎
- Hence Theorem A follows from Theorem A′. For class (1,0) the solutions of G with m even are ±(0,1) (c = 1) and
  ±(2,1) (c = 9); they give F_cl(0,1) = 1 and F_cl(1,0) = 9. The solutions (1,0) (c = 4) and (1,−1) (c = 9), with m
  odd, give no solution of F_cl. For the other two classes G = F_cl.

**2.4 Units.** Each field gets two units of norm +1 (in ω-coordinates):

| class | ε₁ | ε₂ |
|---|---|---|
| (1,0) | (0, 1, −1, 0) = (−1 + β − β²)/3 | (0, −2, 1, 0) = (1 − 4β + β²)/3 |
| (0,1) | (1, −1, 1, 0) = (β − 3)²/6 | (1, 0, −2, 0) = −β²/3 |
| (1,1) | (0, 0, −1, 0) = −β²/6 | (8, −4, 3, 0) = 8 − 4β + β²/2 |

- They were found by a lattice search (LLL on weighted embeddings, Lagrange-reduced in log space), which is not part of
  the proof. Their norms are checked exactly.
- They are multiplicatively independent by Remark 3.2(i).
- As a by-product of §3, O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩ in all three orders (Remark 3.2(iii)). The pair is fundamental, but that
  is never used.

## 3. Coverage: finitely many unit orbits (the two-unit descent)
Script: [01_coverage.py](scripts/01_coverage.py) → [.out.txt](scripts/01_coverage.out.txt). This section fills the
slot `ghost.unit.multi_unit` of [ghost.unit](../../methods/ghost.unit.md) (see
[ghost.unit.multi_unit](../../methods/ghost.unit.multi_unit.md)). Sections 3.1 and 3.2 are those of
[triangular_tetrahedral/proof.md](../triangular_tetrahedral/proof.md); they use only that L is a quartic field of
signature (2,1), that O is an order of L, and that the units have norm 1. They are repeated so that this file stands
alone.

**3.1 Sizes.** Let σ₁ < σ₂ be the real embeddings and σ₃ a complex one. For x ∈ O ⊗ R put s₁ = σ₁(x)², s₂ = σ₂(x)²,
s₃ = |σ₃(x)|² and E(x) = s₁ + s₂ + s₃. Then:
- (a) N(x)² = s₁s₂s₃².
- (b) E is a positive definite quadratic form. It is the squared length of x ↦ (σ₁x, σ₂x, Re σ₃x, Im σ₃x), which is an
  R-linear isomorphism O⊗R → R⁴, because g has four distinct roots (a Vandermonde determinant). So only finitely many
  lattice points have E ≤ T.
- (c) For a unit u put f_v(u) = s_v(u), the *size factors* of u. Then s_v(ux) = f_v(u)s_v(x).
- *Rigour.*
  - σ₁ and σ₂ are enclosed in rational intervals of width < 2⁻²⁰⁰ by exact sign changes of g, and Sturm's count
    confirms that g has exactly two real roots. 01 prints both, for each class.
  - The complex pair β₃, β̄₃ enters only through the real factor t² − t₃t + n₃ of g, where t₃ = β₃ + β̄₃ = −a₃ − σ₁ − σ₂
    and n₃ = β₃β̄₃ = a₀/(σ₁σ₂).
  - If x ≡ A + Bt mod (t² − t₃t + n₃), then s₃ = A² + t₃AB + n₃B².
  - Every quantity below is an interval with rational endpoints ([core.py](scripts/core.py), `Iv`, `Embeddings`).

**3.2 Lemma (reduced representatives).** Let Γ be a group of units of O of norm 1 and U ⊂ Γ finite. Let f_v⁺(u) ≥ f_v(u)
be rational, and let P⁺ be the polygon
{r ∈ R³ : r ≥ 0, r₁+r₂+r₃ = 1, Σ_v r_v(f_v⁺(u) − 1) ≥ 0 for all u ∈ U}.
Suppose every vertex of P⁺ has all coordinates positive, and let b_min = min r₁r₂r₃² over the vertices. Then every
z ∈ O with N(z) = c ≠ 0 is z = γw with γ ∈ Γ, N(w) = c and E(w)⁴ ≤ c²/b_min.

*Proof.*
1. By 3.1(b), the set Σ_z = {γz : γ ∈ Γ, E(γz) ≤ E(z)} is finite. It contains z. Let w be an element of Σ_z with least
   E.
2. For u ∈ U, uw ∈ Γz. If E(uw) < E(w), then uw ∈ Σ_z has smaller E than w, which is impossible. So
   0 ≤ E(uw) − E(w) = Σ_v s_v(w)(f_v(u) − 1), by 3.1(c).
3. Since s_v(w) ≥ 0 and f_v⁺ ≥ f_v, also Σ_v s_v(w)(f_v⁺(u) − 1) ≥ 0. So r := s(w)/E(w) lies in P⁺. Note E(w) > 0,
   because w ≠ 0.
4. *π(r) := r₁r₂r₃² ≥ b_min on P⁺.*
   - P⁺ is a bounded polytope, so it is the convex hull of its vertices (Minkowski). The vertices are the feasible
     solutions of the nonsingular systems formed by r₁+r₂+r₃ = 1 and two constraints taken as equalities; the script
     solves all such pairs.
   - The vertices have positive coordinates, so every point of P⁺ does.
   - On P⁺, log π = log r₁ + log r₂ + 2 log r₃ is a sum of logarithms of positive affine functions, hence concave.
   - So for r = Σ λ_i v_i, a convex combination of vertices, log π(r) ≥ Σ λ_i log π(v_i) ≥ log b_min.
5. By 3.1(a), c² = N(w)² = s₁s₂s₃² = E(w)⁴π(r) ≥ E(w)⁴b_min. ∎

*Remarks.*
- **(i) Independence.** If every vertex of P⁺ is positive, then the log size vectors ℓ(u) = (log f_v(u))_v, u ∈ U, span
  the plane H = {y : y₁ + y₂ + 2y₃ = 0}. In particular ⟨ε₁, ε₂⟩ ≅ Z², so ε₁, ε₂ are multiplicatively independent and −1
  ∉ ⟨ε₁, ε₂⟩. *Proof:* each ℓ(u) lies in H, since Σ d_v log f_v(u) = log N(u)² = 0 with d = (1,1,2). If all ℓ(u) lay on
  one line Rd ⊂ H (d = 0 allowed), note that p = (1,1,2) is orthogonal to d, since d ∈ H = p^⊥. Choose q ⊥ d not
  parallel to p, with a negative coordinate, and move from p along q to the first point r = p + tq with a zero
  coordinate. Then r ≥ 0, r ≠ 0, r ⊥ d, and after normalizing Σ r_v = 1, for every u ∈ U, Σ r_v(f_v⁺(u) − 1) ≥ Σ
  r_v(f_v(u) − 1) ≥ ⟨r, ℓ(u)⟩ = 0 by e^y − 1 ≥ y. So r ∈ P⁺ has a zero coordinate, against step 4, which shows that
  every point of P⁺ is positive. A relation ε₁^aε₂^b = ±1 with (a,b) ≠ (0,0) would put all ℓ(u) on a line.
- **(ii) Variants.** A weaker bound needs only linearity: s_v(w) ≥ κ_vE(w), with κ_v the least v-coordinate of a
  vertex, gives E(w)⁴ ≤ c²/(κ₁κ₂κ₃²); [verify_coverage.py](scripts/verify_coverage.py) uses it as an independent
  cross-check. Step 4 also has a proof without logarithms. Along a segment whose endpoints have π ≥ b_min, expand
  the product binomially and apply AM–GM in the root-free form "y₁⋯y_N ≥ μ^N implies Σy_i ≥ Nμ" (proved by
  induction); so π ≥ b_min on the segment. Apply this first to the edges of P⁺, then to segments from a vertex to a
  point of an edge: P⁺ is a convex polygon, so these segments cover it.
- **(iii) The unit group.** Let T₁ be the least integer with T₁⁴ ≥ 1/b_min (T₁ = 6, 102, 201 in the three classes:
  the T of 3.3 for c = 1). With c = ±1 (the bound depends only on c²), every unit is γw with E(w) ≤ T₁. In each
  field the elements of norm 1 with E ≤ T₁ are all ±ε₁^aε₂^b, and none has norm −1
  ([01_coverage.out.txt](scripts/01_coverage.out.txt)). With (i), O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩. Nothing below uses this.
- **(iv) Provenance.** Koecher's reduction theory of unary forms; see
  [triangular_tetrahedral/proof.md](../triangular_tetrahedral/proof.md), Remark 3.2(iv), and its literature check.

**3.3 The certificates.** Take Γ = ⟨ε₁, ε₂⟩ and U = {ε₁^aε₂^b : a, b ∈ {−1,0,1}, (a,b) ≠ (0,0)}. Each f_v⁺ is the
upper endpoint of an interval of width < 10⁻⁵².

| class | vertices of P⁺ | b_min (rounded) | E(w) ≤ (c²/b_min)^{1/4} (rounded up) | T for c = 1, 4, 9, 36 |
|---|---|---|---|---|
| (1,0) | 6, all positive | 1.26·10⁻³ | 5.31·√c | 6, 11, 16, 32 |
| (0,1) | 6, all positive | 9.41·10⁻⁹ | 101.6·√c | 102, 204, 305, 610 |
| (1,1) | 6, all positive | 6.16·10⁻¹⁰ | 200.8·√c | 201, 402, 603, 1205 |

- T is the least integer with T⁴ ≥ c²/b_min.
- The least vertex coordinates are 1.69·10⁻², 2.59·10⁻⁵ and 5.71·10⁻⁸ in the three classes; positivity is checked
  exactly.
- *Enumeration.*
  - Let G* be the exact Gram matrix of E in the basis ω, Ĝ the symmetric rational matrix formed from the upper
    triangle of the midpoints of its interval enclosure, and ρ the largest radius.
  - Every pivot of the exact LDLᵀ elimination of Ĝ − 4ρI is positive, so it is positive definite (01 prints the
    four pivots of each class).
  - xᵀ(Ĝ − 4ρI)x ≤ xᵀG*x = E(x) for all real x, because |xᵀ(G* − Ĝ)x| ≤ ρ(Σ|x_i|)² ≤ 4ρΣx_i².
  - So the exact Fincke–Pohst enumeration of xᵀ(Ĝ − 4ρI)x ≤ T over Q contains every x with E(x) ≤ T. Each float-guided
    range of the enumeration is certified exactly: both endpoints fail the test, on either side of the centre. 01
    prints the number of ranges of each enumeration (1836 for the largest).
- The elements of norm exactly c are then merged into orbits of ⟨−1, ε₁, ε₂⟩. Two elements are merged only when their
  quotient is verified exactly to be ±ε₁^aε₂^b.
- The largest enumeration is class (1,1), c = 36: 38486 nonzero lattice points with xᵀ(Ĝ − 4ρI)x ≤ 1205. They are
  exactly those with E(x) ≤ 1205: [verify_coverage.py](scripts/verify_coverage.py) enumerates E(x) ≤ T independently
  and finds the same counts at every bound of 01.

**3.4 Proposition (coverage).** Every z ∈ O with N(z) = c ∈ {1,4,9,36} is ±ε₁^aε₂^b·w for one of 12 seeds w, one per
class and c:

| class | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|
| (1,0) | (−1,0,0,0) | (1,0,2,−1) | (1,−2,3,−1) | (−2,−2,3,−1) |
| (0,1) | (−1,0,0,0) | (1,1,−2,0) | (0,−1,0,0) | (3,−1,0,0) |
| (1,1) | (−1,0,0,0) | (4,−1,0,0) | (0,2,−3,0) | (0,−1,0,0) |

*Proof.* Lemma 3.2 with the certificates of 3.3, and the enumeration of 3.3. ∎

## 4. Exclusion: the zeros of the markers on each orbit
Script: [02_exclusion.py](scripts/02_exclusion.py) → [.out.txt](scripts/02_exclusion.out.txt), with the helpers in
[lifting2.py](scripts/lifting2.py) (a verbatim copy of the template's). Write ε^{(a,b)} := ε₁^aε₂^b. Lemmas 4.1–4.4 are
those of [triangular_tetrahedral/proof.md](../triangular_tetrahedral/proof.md) §4, repeated here; the basis remark of
4.3 is sharpened, and 4.4 has a paragraph on the steps of 02.

Fix a seed w.
- *−1 ∉ ⟨ε₁, ε₂⟩.* This follows from Remark 3.2(i). 02 also checks it directly: −1 is not a residue of ⟨ε₁, ε₂⟩
  modulo 8 for class (1,0), and modulo 4 for the classes (0,1) and (1,1).
- *Injectivity.* (a,b) ↦ ε^{(a,b)}w is injective, by Remark 3.2(i).
- So the solution pairs ±z of the orbit correspond one-to-one to Z_w := {(a,b) ∈ Z² : x₂(ε^{(a,b)}w) = x₃(ε^{(a,b)}w) = 0}.
  Theorem A′ needs only the inclusion "every solution is ±ε^{(a,b)}w with (a,b) ∈ Z_w"; together with Z_w = Z₀ below
  that accounts for every solution even without injectivity.

**4.1 Lemma (kernel lattice).** Let q ≥ 3 and Λ_q = {(a,b) : ε^{(a,b)} ≡ 1 (mod qO)}. Let h be the order of ε₁ mod q,
k the least integer > 0 with ε₂^k ∈ ⟨ε₁⟩ mod q, and j an integer with ε₂^k ≡ ε₁^j (mod q). Then
Λ_q = Z(h,0) + Z(−j,k).

*Proof.* The b with ε₂^b ∈ ⟨ε₁ mod q⟩ form the subgroup kZ. If (a,b) ∈ Λ_q, then b = sk, so ε₁^{a+sj} ≡ 1 and
h | a + sj; hence (a,b) = ((a+sj)/h)·(h,0) + s·(−j,k). Conversely both generators lie in Λ_q. ∎

**4.2 Classes and coset systems.**
- (a) *Residues.* On an exponent class (a₀,b₀) + Λ_q every element ε^{(a,b)}w is ≡ ε^{(a₀,b₀)}w (mod qO). So if
  either marker of ε^{(a₀,b₀)}w is ≢ 0 (mod q), the class contains no point of Z_w. The other classes are the
  *candidates*.
- (b) *The coset system.* Put η₁ = ε₁^h and η₂ = ε₁^{−j}ε₂^k: units of norm 1 with η₁ ≡ η₂ ≡ 1 (mod qO), whose
  multiplication matrices are T₁ = I + qK₁ and T₂ = I + qK₂, integer matrices that commute and have determinant 1.
  For z₀ = ε^{(a₀,b₀)}w, the map (s,s') ↦ (a₀,b₀) + s(h,0) + s'(−j,k) is a bijection from Z² onto the exponent class,
  the corresponding elements are T₁^sT₂^{s'}z₀, and the markers along the class are f_r(s,s') = l_rT₁^sT₂^{s'}z₀
  (r = 1, 2) with l₁ = (0,0,1,0), l₂ = (0,0,0,1). This is the setting of
  [lifting.binomial.multi](../../methods/lifting.binomial.multi/proof.md), §1 and Theorem 10, with q = p^e, and e ≥ 2
  when p = 2.

**4.3 Proposition (Jacobian criterion; Proposition 17 of lifting.binomial.multi).** Let z₀ be a common zero. Put
x_r = l_rK₁z₀, y_r = l_rK₂z₀ (r = 1, 2), λ = min ν_p(x₁,x₂,y₁,y₂), J = x₁y₂ − x₂y₁, δ = ν_p(J) − λ and ι = [p = 2]. If
λ ≤ e−1−ι and δ ≤ e−1−ι, then z₀ is the only common zero in Z², that is the only point of Z_w in its class mod Λ_q.
- The proof uses only finite integer polynomials: after a relabelling with ν(x₁) = λ, the stable reductions of f₁ and
  of f₂ − (x₂/x₁)f₁ are a linear form āX + b̄Y with ā ≠ 0 and ūY, because e + λ < 2e − ι and e + δ < 2e − ι; their
  common zero locus mod p has colength one.
- *On the basis.* A change of basis of Λ_q replaces K₁, K₂ by GL₂(Z)-combinations of K₁, K₂ plus q times integral
  matrices, so (x, y) changes by GL₂(Z) modulo q, and J changes by a sign modulo p^{e+min(λ,e)}. So λ is unchanged when
  λ < e. Since ν(J) ≥ 2λ, ν(J) < e + λ forces λ < e, so ν(J) is unchanged when ν(J) < e + λ. The hypotheses of 4.3 give
  ν(J) ≤ e − 1 − ι + λ < e + λ, so a positive verdict of 4.3 does not depend on the basis. Here λ = 0 and ν(J) < e at
  every zero used: ν(J) ≤ 1 with e = 3 at q = 8, and ν(J) = 0 with e = 1 or 2 at q = 23, 43 and 49.
  [verify_exclusion.py](scripts/verify_exclusion.py) recomputes λ and ν(J) in four bases
  ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), negative control (d) and the table of Proposition 17
  instances): they agree at every zero used, and differ only at two failing pairs (q = 7, class (1,1)), where ν(J) ≥ e.

**4.4 Lemma (sieve).** Let q₁, …, q_s be prime powers ≥ 3 (as in 4.2), L = ∩Λ_{q_i}, and Z₀ ⊂ Z_w the set of exponents
of the known solutions. Every point of Z_w lies in a class v mod L whose image mod each Λ_{q_i} is a candidate.
- (a) If v contains a point z ∈ Z₀, and 4.3 holds at z for some q_i, then v ∩ Z_w = {z}.
- (b) If v contains no point of Z₀, but lies inside the Λ_{q_i}-class of some z ∈ Z₀ at which 4.3 holds for that same
  q_i, then v ∩ Z_w = ∅.
- So if every candidate class mod L is of kind (a) or (b), then Z_w = Z₀.

*Proof.* The first sentence is 4.2(a). For (a) and (b), the Λ_{q_i}-class of z contains v, and its only point of Z_w
is z, by 4.3. ∎

*Steps.* (a) and (b) hold verbatim, with the same proof, for a class modulo any lattice L′ ⊆ Λ_{q_i}. 02 applies them
at the step lattices L_t = Λ_{q₁} ∩ ⋯ ∩ Λ_{q_t}: at step t it removes the candidate classes mod L_t of kind (a) or
(b), and refines only the others (the *live* classes) at the next step. A point of Z_w lies in a candidate class at
every step, so if no live class is left at the end, its class was removed at some step t, as kind (a); hence it is a
known zero, and Z_w = Z₀.

**4.5 The certificates.** Each known zero has a modulus where 4.3 holds. "Half" marks a solution of G with m odd,
which gives no solution of F_cl (2.3); it is certified like any other zero.

| class | c | known zeros (m,n) of G | moduli used | Prop. 17 instances |
|---|---|---|---|---|
| (1,0) | 1 | (0,1) | 8, 16 | q=8: λ=0, ν(J)=1 |
| (1,0) | 4 | (1,0) half | 8 | q=8: λ=0, ν(J)=0 |
| (1,0) | 9 | (1,−1) half, (2,1) | 8, 16 | q=8: ν(J)=0 for (1,−1), ν(J)=1 for (2,1); λ=0 |
| (1,0) | 36 | none | 4 | (no candidates mod 4) |
| (0,1) | 1 | (0,1) | 8, 16 | q=8: λ=0, ν(J)=1 |
| (0,1) | 4 | none | 4 | (no candidates mod 4) |
| (0,1) | 9 | (1,0) | 8, 16 | q=8: λ=0, ν(J)=1 |
| (0,1) | 36 | (1,3) | 43, 13 | q=43: λ=0, ν(J)=0 |
| (1,1) | 1 | (0,1), (2,1) | 4, 7, 23 | q=23: λ=0, ν(J)=0, 0 |
| (1,1) | 4 | (1,4) | 4, 8, 7, 23 | q=23: λ=0, ν(J)=0 (also behind one sibling class) |
| (1,1) | 9 | none | 8 | (no candidates mod 8) |
| (1,1) | 36 | (1,0) | 49 | q=49: λ=0, ν(J)=0 |

- **What the output prints.** For each modulus used, |G_q| = [Z² : Λ_q] and the basis (h,0), (−j,k) of 4.1; after
  each step, the Hermite form of the current lattice L and the candidate counts; every certified class with its
  lattice; every sibling class with its witness; and every Proposition 17 instance used. The largest lattice index
  reached is 77,616 (class (0,1), c = 36). The only sibling class removed at 02's steps is in class (1,1), c = 4: the class
  (88, 44) modulo the step lattice (264, 176, 88), emptied by Proposition 17 at (q, zero) = (23, (1,4)).
- **Margins.** The proof of Proposition 17 uses e + λ < 2e − ι and e + δ < 2e − ι. Eight of the eleven instances have
  δ = e − 1 − ι: the four at q = 8 with ν(J) = 1 (3 + δ = 4 < 5), and the four at q = 23 and 43 (e = 1: 1 < 2). Three
  have slack: (8, (1,0)) and (8, (1,−1)) in class (1,0), and (49, (1,0)) in class (1,1) (2 < 4). The inequality
  e + λ < 2e − ι (λ = 0 throughout) holds with margin 1 at q = 23 and 43, and with margin 2 at q = 8 and 49.
- **Unnecessary moduli.** In class (1,1), c = 4, the moduli 4 and 8 add nothing: [7, 23] alone certifies the orbit
  ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), other modulus sets). The certificate is left as the
  greedy search found it.
- **A 2-adic zero that is not an integer zero.** In the orbit (0,1), c = 36, the markers have a common zero at a 2-adic
  exponent α ∈ Z₂² with α ≡ (46,6) mod (64,32) ([check_2adic_zero.py](scripts/check_2adic_zero.py) →
  [.out.txt](scripts/check_2adic_zero.out.txt)). The certificate lifts the candidates to 2¹²⁰ and checks Hensel's
  condition ν(F) > 2ν(det DF) for the markers F as power series in the exponents, with ν(F) ≥ 120 and ν(det DF) = 7; the
  multivariable Hensel lemma, which is p-adic analysis, is used only in this remark, and nothing else uses the remark.
  - α is not an integer, because the only integer zero of the orbit is the exponent (0,0) of (1,3) (Theorem A′).
  - α lies outside the Λ₄-class of (0,0). So for every k the class of α mod Λ_{2^k} is a candidate that is neither
    certified nor a sibling: no sieve with moduli 2^k alone can finish this orbit. The odd moduli 43 and 13 do
    ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), other modulus sets, which also follows the candidates
    to 2⁴⁸).
  - This is the obstruction foreseen in [triangular_tetrahedral](../triangular_tetrahedral/proof.md) §6.4, met here
    for one prime and avoided with others.
- **Phantoms.** −1 is a residue of ⟨ε₁, ε₂⟩ modulo 13 in class (0,1): −1 ≡ ε^{(38,2)} (mod 13). So the exponent
  class (38,2) mod Λ₁₃, whose element is ≡ −w (mod 13), the class of −(1,3), is a candidate modulo 13. Of its 231
  subclasses modulo the final lattice (924, 168, 84), 230 fail the candidate test modulo 43, and one lies in the
  Λ₄₃-class of the exponent (0,0) of (1,3), certified at the first step. So in the order [43, 13] it never becomes
  live ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), negative control (c)).

**Proof of Theorem A′.**
- By 2.2, a solution of G(m,n) = c is an element z ∈ O with N(z) = c and vanishing markers.
- By 3.4, z or −z is ε^{(a,b)}w for one of the 12 seeds w, with (a,b) ∈ Z_w.
- By 4.5 and 4.4, Z_w = Z₀, the exponents of the known solutions.
- So z is ± a known solution, and (m,n) is read off from z = n − mβ. ∎

**Proof of Theorem A.** Theorem A′ and Lemma 2.3. ∎

## 5. Reconstruction
Script: [03_reconstruction.py](scripts/03_reconstruction.py) → [.out.txt](scripts/03_reconstruction.out.txt).

- *Points of the non-trivial classes.* Each solution of Theorem A has g = √c (1.3) and gives the point
  (X₄/g², ±Y₆/g³), and ±(m,n) ↔ (X, ±Y) is one-to-one:

  | class | c | ±(m,n) of F_cl | g | x | \|y\| |
  |---|---|---|---|---|---|
  | (1,0) | 1 | (0,1) | 1 | 3 | 9 |
  | (1,0) | 9 | (1,0) | 3 | 39 | 243 |
  | (0,1) | 1 | (0,1) | 1 | −3 | 9 |
  | (0,1) | 9 | (1,0) | 3 | 33 | 189 |
  | (0,1) | 36 | (1,3) | 6 | 24 | 117 |
  | (1,1) | 1 | (0,1) | 1 | 0 | 9 |
  | (1,1) | 1 | (2,1) | 1 | 5112 | 365499 |
  | (1,1) | 4 | (1,4) | 2 | 513 | 11619 |
  | (1,1) | 36 | (1,0) | 6 | 9 | 27 |

- *The extra solutions of the field form.* The solutions (1,0) (c = 4) and (1,−1) (c = 9) of G in class (1,0) have m
  odd. They give F_cl(1, −1) = 64 and F_cl(1, −3) = 144, values 16c outside {1, 4, 9, 36}, and only the non-integral
  rational points (57/16, ±621/64) and (−15/4, ±63/8).
- *The trivial class* gives x ∈ {−5, 7, 1099} (1.3).
- *Theorem B.* By (H), now Theorem H of [hypothesis_H.md](hypothesis_H.md), the integral points are the 24 points
  with x ∈ {−5, −3, 0, 3, 7, 9, 24, 33, 39, 513, 1099,
  5112}. They agree with LMFDB's list and with the Mordell–Weil decompositions of the descent file, and the class of
  each point is (i mod 2, j mod 2) for P = iG₁ + jG₂.
- *Transfer.* T(n) = P(m) holds exactly when (x, y) = (6m+3, 18n+9) ∈ E (0.1). For x = 6t + 3,
  x³ − 9x + 81 = 81(8P(t) + 1), so 9 | y and y/9 is odd: every integral point with x ≡ 3 (mod 6) gives a solution.
  These are x ∈ {−3, 3, 9, 33, 39, 513}, so m ∈ {−1, 0, 1, 5, 6, 85}, with the common values 0, 0, 1, 55, 91 and
  208335. The involution y ↦ −y is n ↦ −1 − n. With m, n ≥ 0 the solutions are (m,n) = (0,0), (1,1), (5,10), (6,13)
  and (85,645).
- *Uchiyama's model.* 3Y² = X³ − X + 3 is E through (x, y) = (3X, 9Y), and (X, Y) = (2m+1, 2n+1) turns it into
  T(n) = P(m). By Theorem B, its integral points are the 18 points with X ∈ {−1, 0, 1, 3, 8, 11, 13, 171, 1704}.

## 6. Status, comparison and novelty

**6.1 Status.**
- Theorems A′ and A are proved here in full, by finite elementary computations. They were verified in-session by
  independent scripts and refereed in-session by agents ([review.md](review.md), rounds 1 and 2); there is no outside
  audit.
- Theorem B is proved in full: §1 with Theorem H of [hypothesis_H.md](hypothesis_H.md), which replaces the rank and
  the facts N1, N2, and §5. hypothesis_H.md, Lemma 1.3′ and the other changes for it were verified by independent
  scripts and refereed in-session too (review.md, round 3); there is no outside audit.

**6.2 What R1 needed, and what it did not.**
- *Coverage* needed no unit theorem and no fundamental units: two units, a polygon of eight linear inequalities, and
  an exact enumeration of at most 38486 lattice points.
- *Exclusion* needed only Proposition 17 with a congruence sieve over q ≤ 49. The frame criterion (Theorem 18) and
  saturation (Theorem 11) were not needed.
- *New compared with 1611.a1.*
  - One class form is not the norm form of a generator with O ∩ (Q + Qβ) = Z + Zβ, so a field form and the
    substitution lemma 2.3 are needed. The two solutions of the field form with m odd are certified like the others.
  - One orbit has a common 2-adic zero of the markers that is not an integer zero (4.5, with a Hensel certificate).
    No sieve by powers of 2 alone can finish it, and odd moduli do. This is the obstruction foreseen in
    [triangular_tetrahedral/proof.md](../triangular_tetrahedral/proof.md) §6.4, met here for one prime and avoided
    with others.

**6.3 Uchiyama (1984), read in full** ([literature_check.md](literature_check.md) §1, with the exact checks in
[lit_uchiyama.out.txt](scripts/lit_uchiyama.out.txt)).
- *His model and reduction.* He works on 3Y² = X³ − X + 3, whose integral points are the integral points (3X, 9Y)
  of E with 3 | x. His cubic ring Z[λ], λ³ − λ + 3 = 0, is our O_K (λ = θ/3).
  - With the class number 1 of Q(λ), stated without a source, and its fundamental unit, taken from the
    Delone–Faddeev table (his reference [7, p. 113]), he factorizes X − λ.
  - His cases X ≡ 0, 1, −1 (mod 3) are our classes (1,1), (1,0) and (0,1). So his §1 proves (H) for the points with
    3 | x, which are all that T(n) = P(m) needs; it uses ideals and a class number. hypothesis_H.md replaces his class
    number (stated without a source) and his tabulated unit: a Euclidean algorithm in Z[λ] gives the divisibility
    theory, and a descent gives the unit.
- *His equations are ours.* His twelve sub-case equations are exactly our twelve equations F_cl(m,n) = c, through
  explicit GL₂(Q) substitutions. He excludes four of them by congruences; 02 removes three of these by a congruence
  alone, and the fourth, ((1,0), 4), has only the extra field-form solution G(1,0) = 4. Eight Thue equations remain.
- *His fields and units are ours.* His three quartic fields are our three. In each, his fundamental pair (Berwick's
  criterion) and ours generate the same group modulo ±1, so by Remark 3.2(iii) his unit groups are right.
- *His exclusion.* It is a case analysis on the 2-adic valuations of the exponents, written as finite congruences
  modulo powers of 2 on the two marker coordinates, mostly one marker at a time. There are no p-adic series. This is the
  completion-free form of the two-index argument, as already noted for 1611.a1.
- *Errors in print* (all checked exactly; 117 of 155 tested statements are true):
  - his integral bases of two of the three fields span orders of index 2 in the maximal order, so his discriminant
    −2⁸·3²·239 should be −2⁶·3²·239;
  - the printed units of Props. 3, 7 and 8 are not units (wrong matrices M in Props. 3 and 8; in Prop. 7 M is right,
    but the units were not transported correctly); corrected units are given;
  - Prop. 1, step 4: the printed congruence modulo 2^{h+5} is false in most tested cases, and its conclusion fails
    near the solutions (2,0) and (0,0) of the exponents. What is missing there is a Jacobian-type criterion, which
    Proposition 17 supplies;
  - Prop. 2's odd case and Prop. 4's lemma are not proved as printed, and Props. 3, 5–8 are sketched "on the same
    lines";
  - Prop. 8 cannot be completed by congruences modulo powers of 2 at all, although it is said to go "just as in the
    proof of Proposition 2": its equation is our orbit (0,1), c = 36, whose markers have a 2-adic zero that is not an
    integer zero (4.5). Odd moduli are needed. Prop. 3, by contrast (our orbit (1,1), c = 4), is repaired 2-adically by
    Proposition 17 at q = 16 ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), other modulus sets).

  The statements of all eight propositions are true (Theorem A here). So the theorem is right: Theorem A here with
  his §1 (given h_K = 1 and the unit of K) proves it; this card proves it without that input (Theorem H). His
  reduction is a sound
  proof given the class number and the unit, and his solution of the Thue equations is incomplete in print.

**6.4 Other proofs** ([literature_check.md](literature_check.md) §§2–3).
- Avanesov, Volzh. Mat. Sb. 8 (1971) 3–6, and Finkelstein–London, J. Number Theory 4 (1972) 455–462, both solved the
  problem; neither was read. Finkelstein–London's references (Avanesov 1968, Billevich's units, Delone–Faddeev) suggest
  the same toolkit.
- Avanesov–Gusev (2012) describe a height bound (Dem'yanenko) followed by a computer search.
- Integral-point software (elliptic logarithms) reproduces the 24 points; no paper treating this curve that way was
  found.
- No elementary proof was found.

**6.5 Novelty.** None is claimed. By the project's criterion this is a reconstruction of Uchiyama's proof: the same
fields, equations and units, and the same kind of completion-free two-index exclusion, with odd moduli where powers
of 2 cannot suffice. What this card adds:
- a complete proof of all twelve Thue equations, with exact, re-runnable certificates, where Uchiyama's printed
  argument has gaps;
- the errors in Uchiyama's printed argument, recorded in [literature_check.md](literature_check.md) §1.7;
- the units are not assumed fundamental: coverage comes from Lemma 3.2, and fundamentality is a by-product, proved
  without the unit theorem;
- the reduction uses no ideals and no class number: (H) is proved by a Euclidean algorithm in Z[λ] with an exact
  covering certificate, and a one-unit descent ([hypothesis_H.md](hypothesis_H.md)).

**6.6 What remains open.**
- *Done (2026-10-15, reviewed 2026-10-26):* an elementary proof of (H) for integral points, in
  [hypothesis_H.md](hypothesis_H.md). Z[λ] is norm-Euclidean (Taylor's theorem), certified by 116 exact boxes (the
  point fixed by ε needs 8/9, the value in Cavallar–Lemmermeyer's table); with the units ±εᵏ by a one-unit descent and
  a congruence modulo 4, this gives (H).
- The exclusion sieve has no a priori termination guarantee; one orbit here shows that the powers of a single prime
  can fail.

## 7. Files
- [core.py](scripts/core.py): fields, orders, exact ring arithmetic, rigorous intervals, units, descent maps, the
  maps between field forms and class forms.
- [00_fields.py](scripts/00_fields.py): fields, maximal orders, markers, the substitution lemma, units, known
  solutions.
- [01_coverage.py](scripts/01_coverage.py): Lemma 3.2 certificates, enumeration, orbits, unit group.
- [lifting2.py](scripts/lifting2.py) and [02_exclusion.py](scripts/02_exclusion.py): kernel lattices, coset systems,
  Proposition 17, sieve.
- [03_reconstruction.py](scripts/03_reconstruction.py): back to E and to T(n) = P(m).
- Independent checks: [verify_coverage.py](scripts/verify_coverage.py),
  [verify_exclusion.py](scripts/verify_exclusion.py), [verify_bruteforce.py](scripts/verify_bruteforce.py).
- [lit_uchiyama.py](scripts/lit_uchiyama.py): exact checks of Uchiyama (1984): model, reduction, the twelve sub-case
  equations, fields, bases, units, printed congruences.
- [check_2adic_zero.py](scripts/check_2adic_zero.py): the Hensel certificate for the 2-adic zero of 4.5 (outside the
  proof).
- Hypothesis (H), in [hypothesis_H.md](hypothesis_H.md):
  - [hring.py](scripts/hring.py): the ring Z[λ]; [hlib.py](scripts/hlib.py): shared helpers;
  - [04_ring.py](scripts/04_ring.py): the ring of triples;
  - [05_cover.py](scripts/05_cover.py) with [euclid_boxes.txt](scripts/euclid_boxes.txt): the Euclidean division, and
    the lower bound 8/9 at the point fixed by ε;
  - [06_units.py](scripts/06_units.py): the units;
  - [07_classes.py](scripts/07_classes.py): the class of X − θ;
  - [cover_search.py](scripts/cover_search.py): how the table was found (not used by the proof);
  - independent checks: [verify_H_units.py](scripts/verify_H_units.py), [verify_H_cover.py](scripts/verify_H_cover.py),
    [verify_H_classes.py](scripts/verify_H_classes.py) (review.md, round 3);
  - [lit_H_lemmermeyer.py](scripts/lit_H_lemmermeyer.py), [lit_H_lmfdb.py](scripts/lit_H_lmfdb.py): the identifications
    with the tables cited in [literature_check.md](literature_check.md) §6 (not used by the proof).
