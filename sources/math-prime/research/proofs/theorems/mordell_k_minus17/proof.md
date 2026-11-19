# y² = x³ + 17 (10404.f2) by strategy R1: two-unit ghosts and finite-polynomial exclusion

Written 2026-10-15 by Claude (Opus 5.5) as the third run of strategy R1 of
[GENERALIZATIONS.md](../../GENERALIZATIONS.md) §3.2, after [triangular_tetrahedral](../triangular_tetrahedral/proof.md)
(1611.a1), together with [triangular_square_pyramidal](../triangular_square_pyramidal/proof.md) (8604.a1). The
machinery is the same; the data are new. Scripts: [scripts/](scripts/). Each script has its output
beside it as `*.out.txt`. Every check is exact, and floating point only seeds searches. Hypothesis (H) of §1.2, which
this proof first took from the rank (or from Cassels's theorem), is proved in [hypothesis_H.md](hypothesis_H.md) (added
2026-10-15, reviewed in round 3 of [review.md](review.md)).

## 0. Statements

**0.1 The curve.** E: y² = x³ + 17, that is y² = x³ + k with k = 17, or x² − 17 = y³ after renaming (x² + k = y³
with k = −17 in the catalogue's convention, hence the slug). This is
LMFDB 10404.f2, whose model is already monic, so X = x and Y = y throughout. Its discriminant is
−124848 = −2⁴·3³·17², and j = 0. LMFDB gives rank 2, an mwrank value (see 0.4).

**0.2 Class forms and field forms.** The 2-descent of E (§1) attaches to each non-trivial class a binary quartic
*class form* F_cl. The Thue equations actually solved are those of a monic *field form* G, the norm form
N(n − mβ) of a generator β whose order O satisfies O ∩ (Q + Qβ) = Z + Zβ (§2.2). For class (1,0), G = F_cl. For
class (0,1), α/2 is integral for a root α of F_cl(1,t). For class (1,1), F_cl is not monic (F_cl(0,1) = 4). See §2.3.

| class | μ | F_cl(m,n) | right sides | G(m,n) | β | relation |
|---|---|---|---|---|---|---|
| (1,0) | −2 − θ | −12m⁴ − 24m³n + 12m²n² + n⁴ | 1, 4 | the same | α | G = F_cl |
| (0,1) | 4 − θ | 16m⁴ − 8m³n − 8mn³ + n⁴ | 1, 4 | m⁴ − m³n − 4mn³ + n⁴ | α/2 | F_cl(m,n) = G(2m, n) |
| (1,1) | −1 − θ | −3m⁴ + 20m³n − 12m²n² + 4n⁴ | 1, 4, 9, 36 | −12m⁴ + 40m³n − 12m²n² + n⁴ | 2α | 4F_cl(m,n) = G(m, 2n), F_cl(2m, n) = 4G(m,n) |

**Theorem A′ (the eight Thue equations of the field forms).** For each class and each c among its right sides, the
integer solutions of G(m,n) = c are exactly ±(m,n) for the following pairs (m,n):

| class | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|
| (1,0) | (0,1) | (1,2) | — | — |
| (0,1) | (0,1), (1,0), (2,1) | none | — | — |
| (1,1) | (0,1), (3,1) | none | none | (1,2) |

**Theorem A (the class forms).** For each class and each c among its right sides, the integer solutions of
F_cl(m,n) = c are exactly ±(m,n) for:

| class | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|
| (1,0) | (0,1) | (1,2) | — | — |
| (0,1) | (0,1), (1,1) | none | — | — |
| (1,1) | none | (0,1), (6,1) | (1,1) | none |

**Theorem B (integral points).** (Hypothesis (H) of 1.2, on which this rests, is Theorem H of
[hypothesis_H.md](hypothesis_H.md).) The integral points of E are the 16 points with x ∈ {−2, −1, 2, 4, 8, 43, 52, 5234}, with both signs of y. Equivalently, x² − 17 = y³ holds
in integers only for y ∈ {−2, −1, 2, 4, 8, 43, 52, 5234}, with |x| = 3, 4, 5, 9, 23, 282, 375, 378661.

**0.3 What is used.**
- **Theorems A′ and A** use:
  - exact computation with integers and rationals;
  - real numbers only through enclosures in rational intervals, obtained from exact sign changes, with the number of
    real roots fixed by Sturm's theorem;
  - elementary real-variable facts in Lemma 3.2: positive definiteness of the size form E of §3.1 (not the curve),
    that a bounded polygon is the convex hull of its vertices, and concavity of the logarithm of a positive affine
    function (or the root-free arguments of Remark 3.2(ii)); Remark 3.2(i) also uses e^y ≥ 1 + y, but Theorem A′ does
    not need that remark (§4);
  - Proposition 17 of the finite-polynomial lifting lemma
    [lifting.binomial.multi](../../methods/lifting.binomial.multi.md);
  - congruences.

  They use no completion, no p-adic power series, no ideals, no class numbers, no unit theorem and no linear forms in
  logarithms.
- **Theorem B** adds the 2-descent of §1 and [hypothesis_H.md](hypothesis_H.md).
  - **(H)**, the list of four classes, is Theorem H there. It is proved by:
    - a division with multipliers {1, 2} (Dedekind–Hasse) in R = Z + Zθ + Zω (ω = (1 + θ + θ²)/3), certified by 795
      exact boxes; R is not norm-Euclidean;
    - the divisibility theory that follows from that division;
    - the units ±ηᵏ, η = 18 + 7θ, by a one-unit descent;
    - the norm condition, which alone leaves the four classes.

    Neither the rank of E(Q) nor Cassels's theorem (§6.1) is used.
  - **N1 and N2** of the descent script (the index–discriminant formula and p-maximality; integral elements form a
    ring, and the characteristic-polynomial criterion) are replaced by hypothesis_H.md (3.6): R is integrally closed,
    and every element of it is integral. So O_K = R, as §1.3 uses (hypothesis_H.md §6).

  The rest of §1 (conics, parametrization, gcd bounds) is elementary.

**0.4 Rank.** LMFDB lists rank 2 for 10404.f2 with generators (−2, 3) and (4, 9). The proof does not use the rank:
rank ≤ 2 was the original route to (H), which hypothesis_H.md proves directly.

## 1. The reduction to Thue equations (the 2-descent; prior work)
Source: [weierstrass_2descent.py](../../targets/rank2_signatures_2026-10/weierstrass_2descent.py) →
[.out.txt](../../targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt), section `CURVE 10404.f2`. Every item
used below is re-derived by exact arithmetic in [03_reconstruction.out.txt](scripts/03_reconstruction.out.txt).

**1.1 The cubic field and the descent map.**
- Put K = Q(θ) with θ³ + 17 = 0. Its signature is (1,1). O_K = Z + Zθ + Z·(1 + θ + θ²)/3, with disc(O_K) = −867 =
  −3·17² and [O_K : Z[θ]] = 3. That it is the whole ring of integers is [hypothesis_H.md](hypothesis_H.md) (3.6); the
  descent script's N1 and N2 are not used.
- For a point P ≠ O put φ(P) = X(P) − θ in K*/K*². φ is a homomorphism whose kernel contains 2E(Q), by the identity
  (x₁−θ)(x₂−θ)(x₃−θ) = (λθ+ν)² for three collinear points on Y = λX + ν, together with the vertical-line case
  (checked again in 03 on a sample).
- E(Q)[2] = 0, because t³ + 17 is irreducible.

**1.2 Hypothesis (H).** Every integral point P of E has φ(P) ∈ {1, −2−θ, 4−θ, −1−θ}.
- These four classes are φ of O, G₁ = (−2, 3), G₂ = (4, 9) and (−1, 4) = −G₁ − G₂, for the LMFDB generators G₁, G₂.
  They are pairwise distinct: −2−θ, 4−θ, −1−θ and their pairwise products are non-squares, by a Legendre-symbol witness
  at a root of t³ + 17 modulo a prime (11 or 5 for the three classes; checked exactly). Such a witness is constant on a
  class modulo squares because a square root in K of an element of O_K lies in O_K ([hypothesis_H.md](hypothesis_H.md)
  (5.7) and (3.6)). Theorem B does not need the distinctness.
- (H) would follow from rank E(Q) ≤ 2, by N5, since then |E(Q)/2E(Q)| ≤ 4; it is also Cassels's theorem (§6.1).
  These were the original routes.
- **(H) is proved for integral points in [hypothesis_H.md](hypothesis_H.md), Theorem H**, without the rank and without
  ideals: X − θ = μ₀S² with S ∈ O_K and μ₀ ∈ {1, ρσ, η, ηρσ} (ρ = 1 − θ + ω, σ = 1 + θ + ω), and each μ₀ is one of
  the four classes times a square.

**1.3 From a class to a Thue equation.** Let μ be the class and X − θ = μξ².
- N(μ) = 9, 81 and 16 for the classes (1,0), (0,1) and (1,1). Ξ_μ = {ξ : μξ² ∈ O_K} is a lattice with
  O_K ⊆ Ξ_μ ⊆ (1/k)O_K for k = 3, 9, 4, because O_K is integrally closed ([hypothesis_H.md](hypothesis_H.md) (3.6) and
  §6); it contains O_K with index 1, 3 and 4.
- In a basis of Ξ_μ, ξ has coordinates (u,v,w) with A₂ = 0 (a smooth conic), A₁ = −1 and A₀ = X.
- A common factor k of (u,v,w) satisfies k² | 3, because θ-coordinates of elements of O_K lie in (1/3)Z; so k = 1.
- The conic has the rational point ξ = 1. So (u,v,w) = ±𝐗(m,n)/g for coprime (m,n), by Lemma 1.3′ below. Here 𝐗(m,n)
  is the parametrization and g = gcd 𝐗(m,n) divides 4, 4 and 12; that divisibility follows from the polynomial
  identities that the descent script checks.
- Then A₁ = −1 becomes F_cl(m,n) = g², with g ∈ {1, 2} for classes (1,0) and (0,1), and g ∈ {1, 2, 3, 6} for class
  (1,1).
- *g is forced.*
  - Classes (1,0) and (0,1): u² − F_cl is 4 times an integral form and g = gcd(u, 2). So every coprime solution of
    F_cl = c ∈ {1, 4} has g = √c.
  - Class (1,1): 2 | g exactly when F_cl is even, and 3 | g exactly when m ≡ n (mod 3). A coprime solution of
    F_cl = c ∈ {1, 4, 9, 36} has g = √c, unless 3 | c and 3 | n; then g = √c/3 and the solution gives only a rational
    point. The one solution of Theorem A with 3 | c is (1,1), with 3 ∤ n, so every solution of Theorem A gives an
    integral point.
- X = X₄(m,n)/g² and Y = ±Y₆(m,n)/g³; the polynomials are in [core.py](scripts/core.py), `descent_maps`, and
  Y₆² = X₄³ + 17F_cl³ identically.
- *X is an integer.* A₀ has denominators 3 for class (1,1) and the trivial class, so the template's argument does
  not apply there; the argument below covers all classes. An element of O_K
  with θ²-coordinate 0 lies in Z + Zθ, and μξ² = X₄/g² − (F_cl/g²)θ is such an element when F_cl = g².
- *The trivial class* (μ = 1) gives H = 4F₀ = g² with F₀ = n(3m³ + 3m²n + mn² + 2n³) and g | 6, g even, so
  F₀ ∈ {1, 9}. Then n | F₀, and the rational-root test leaves the solutions ±(1,−1) (F₀ = 1) and ±(1,1) (F₀ = 9). The
  second has g = 2, not 6, and gives only the rational point (−8/9, ±109/27). The trivial class gives X = 8 (§5).

So every integral point outside the trivial class gives a coprime solution of F_cl(m,n) = c, c a right side of its
class, in one of the three classes. Theorem A lists all of them.

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
| (1,0) | t⁴ + 12t² − 24t − 12 | β²/2 | (2β + β³)/4 | 8 | −2⁶·3³·17² |
| (0,1) | t⁴ − 4t³ − t + 1 | β² | β³ | 1 | −3³·17² |
| (1,1) | t⁴ − 12t² + 40t − 12 | β²/2 | (10β + 2β² + β³)/12 | 24 | −2⁶·3·17² |

In each case ω₀ = 1 and ω₁ = β. O is the maximal order (p-maximality is tested at every p with p² | disc O: at 2, 3
and 17 for class (1,0); at 3 and 17 for class (0,1), where O = Z[β]; at 2 and 17 for class (1,1), where 3 ∥ disc O),
but maximality is never used. What
is used is that O is a ring containing Z[β] and the units below. The three fields are pairwise non-isomorphic: g
factors mod 5 (class (1,0) against the other two) and mod 11 (classes (0,1), (1,1)) with different degree patterns,
and the three discriminants differ.

**2.2 Markers.** Write x = (x₀, x₁, x₂, x₃) for the coordinates in the basis ω.
- The β²- and β³-coordinates of x are (x₂, x₃)·D for an invertible 2×2 matrix D (lower triangular with diagonal
  (1/2, 1/4), (1, 1) and (1/2, 1/12)).
- Since 1, β, β², β³ are linearly independent (g is irreducible), x ∈ Q + Qβ if and only if x₂ = x₃ = 0, and then
  x = x₀ + x₁β ∈ Z + Zβ.
- n − mβ has coordinates (n, −m, 0, 0).
- So **the solutions of G(m,n) = c are the elements z ∈ O with N(z) = c and x₂(z) = x₃(z) = 0.**

**2.3 From the class forms to the field forms.**
- *Class (0,1).* Let α be a root of f_cl(t) = F_cl(1,t) = t⁴ − 8t³ − 8t + 16. Then β := α/2 is a root of
  g(t) = t⁴ − 4t³ − t + 1, since 16·g(t/2) = f_cl(t). So α/2 is integral, it lies in Q + Qα but not in Z + Zα, and the
  marker lemma 2.2 fails for F_cl itself. With β it holds: here O = Z[β].
  - *Identity.* N(n − mα) = N(n − 2mβ), that is F_cl(m,n) = G(2m, n).
  - *Lemma.* For every integer c, (m,n) ↦ (2m, n) is a bijection from the integer solutions of F_cl(m,n) = c onto the
    integer solutions (m', n') of G(m',n') = c with m' even. *Proof.* It maps Z² bijectively onto
    {(m',n') : m' even}, and the identity carries one equation into the other. ∎
  - The solutions of G with m even are ±(0,1) and ±(2,1) (c = 1). They give F_cl(0,1) = F_cl(1,1) = 1. The solution
    (1,0), with m odd, gives no solution of F_cl.
- *Class (1,1).* F_cl(1,t) = 4t⁴ − 12t² + 20t − 3 is not monic. For β = 2α, with α a root of it,
  g(t) = 4F_cl(1, t/2) = t⁴ − 12t² + 40t − 12 is monic, and G(m,n) = 4F_cl(m, n/2). The identities of forms
  4F_cl(m,n) = G(m, 2n) and F_cl(2m, n) = 4G(m,n) are checked in [00_fields.out.txt](scripts/00_fields.out.txt), §4.
  - *Congruences.* F_cl(m,n) ≡ m⁴ (mod 4) and G(m,n) ≡ n⁴ (mod 4): −3 ≡ 1 (mod 4), and the other coefficients
    (20, −12, 4 for F_cl; −12, 40, −12 for G) are divisible by 4.
  - *Lemma.* Let S = {1, 4, 9, 36}. The maps
    - (m,n) ↦ (m, 2n) for m odd, and (m,n) ↦ (m/2, n) for m even,
    - with inverses (m',n') ↦ (m', n'/2) for n' even, and (m',n') ↦ (2m', n') for n' odd,

    are mutually inverse bijections between the integer solutions of F_cl(m,n) ∈ S and those of G(m',n') ∈ S. The
    value changes by the factor 4 or 1/4.
  - *Proof.* Integer solutions of either equation are primitive, since no fourth power > 1 divides an element of S.
    1. If F_cl(m,n) = c ∈ S with m odd, then c ≡ 1 (mod 4), so c ∈ {1, 9}, and G(m, 2n) = 4c ∈ {4, 36} ⊂ S, with
       n' = 2n even.
    2. If m = 2k is even, then c = F_cl(2k, n) = 4G(k, n), so c ∈ {4, 36} and G(k, n) = c/4 ∈ {1, 9}. Then n is odd,
       because G(k, n) ≡ n⁴ (mod 4) is odd.
    3. Conversely, if G(m',n') = c' ∈ S with n' even, then c' ≡ 0 (mod 4), so c' ∈ {4, 36}, and
       F_cl(m', n'/2) = c'/4 ∈ {1, 9}. Here m' is odd, since (m', n') is primitive.
    4. If n' is odd, then c' ≡ 1 (mod 4), so c' ∈ {1, 9}, and F_cl(2m', n') = 4c' ∈ {4, 36}.
    5. The four maps preserve the parity conditions and compose to the identity in both orders. ∎
  - The solutions of G give: (0,1) and (3,1) (c' = 1, n' odd) → F_cl(0,1) = F_cl(6,1) = 4; (1,2) (c' = 36, n' even) →
    F_cl(1,1) = 9. G = 4 and G = 9 have no solutions, so F_cl = 1 and F_cl = 36 have none.
- Hence Theorem A follows from Theorem A′.

**2.4 Units.** Each field gets two units of norm +1 (in ω-coordinates):

| class | ε₁ | ε₂ |
|---|---|---|
| (1,0) | (2, −2, 1, 0) = 2 − 2β + β²/2 | (1, 2, −2, 0) = 1 + 2β − β² |
| (0,1) | (0, −1, 0, 0) = −β | (1, −4, 5, −1) = 1 − 4β + 5β² − β³ |
| (1,1) | (0, 1, 0, −1) = (2β − 2β² − β³)/12 | (10, −47, −7, 18) = 10 − 32β − β²/2 + 3β³/2 |

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
- **(ii) Variants.**
  - A weaker bound needs only linearity: s_v(w) ≥ κ_vE(w), with κ_v the least v-coordinate of a vertex, gives
    E(w)⁴ ≤ c²/(κ₁κ₂κ₃²). [verify_coverage.py](scripts/verify_coverage.py) uses it where the enumeration is feasible:
    classes (1,0) and (0,1), and c = 1 in class (1,1).
  - For class (1,1) with c = 4, 9, 36 that bound is too weak (κ₁ = 3.4·10⁻⁹). verify_coverage.py then applies the same
    linearity argument on a subdivision of P⁺ into triangles, on each of which κ₁κ₂κ₃² ≥ 0.9·b_min. This gives a
    logarithm-free weak form of step 4, with b_min replaced by β ≥ 0.9·b_min (verify_coverage.out.txt prints
    β = 0.9036, 0.9074 and 0.9058·b_min), and its enumeration, to E ≤ 138179, covers the resulting bounds (at most
    2315).
  - Step 4 also has a proof by binomial expansion. Along a segment whose endpoints have π ≥ b_min, expand the
    product binomially and apply AM–GM in the root-free form "y₁⋯y_N ≥ μ^N implies Σy_i ≥ Nμ" (proved by
    induction); so π ≥ b_min on the segment. Apply this first to the edges of P⁺, then to segments from a vertex to a
    point of an edge: P⁺ is a convex polygon, so these segments cover it.
- **(iii) The unit group.** Let T₁ be the least integer with T₁⁴ ≥ 1/b_min (T₁ = 370, 8, 377 in the three classes:
  the T of 3.3 for c = 1). With c = ±1 (the bound depends only on c²), every unit is γw with E(w) ≤ T₁. In each
  field the elements of norm 1 with E ≤ T₁ are all ±ε₁^aε₂^b, and none has norm −1
  ([01_coverage.out.txt](scripts/01_coverage.out.txt)). With (i), O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩. Nothing below uses this.
- **(iv) Provenance.** Koecher's reduction theory of unary forms; see
  [triangular_tetrahedral/proof.md](../triangular_tetrahedral/proof.md), Remark 3.2(iv), and its literature check.

**3.3 The certificates.** Take Γ = ⟨ε₁, ε₂⟩ and U = {ε₁^aε₂^b : a, b ∈ {−1,0,1}, (a,b) ≠ (0,0)}. Each f_v⁺ is the
upper endpoint of an interval of width < 10⁻⁵¹.

| class | vertices of P⁺ | b_min (rounded) | E(w) ≤ (c²/b_min)^{1/4} (rounded up) | T for each right side c |
|---|---|---|---|---|
| (1,0) | 6, all positive | 5.39·10⁻¹¹ | 369.1·√c | 370 (c = 1), 739 (c = 4) |
| (0,1) | 6, all positive | 3.91·10⁻⁴ | 7.12·√c | 8 (c = 1), 15 (c = 4) |
| (1,1) | 6, all positive | 4.99·10⁻¹¹ | 376.3·√c | 377, 753, 1129, 2258 (c = 1, 4, 9, 36) |

- T is the least integer with T⁴ ≥ c²/b_min.
- The least vertex coordinates are 1.44·10⁻⁶, 3.51·10⁻³ and 3.37·10⁻⁹ in the three classes; positivity is checked
  exactly.
- *Enumeration.*
  - Let G* be the exact Gram matrix of E in the basis ω, Ĝ the symmetric rational matrix formed from the upper
    triangle of the midpoints of its interval enclosure, and ρ the largest radius.
  - Every pivot of the exact LDLᵀ elimination of Ĝ − 4ρI is positive, so it is positive definite (01 prints the
    four pivots of each class).
  - xᵀ(Ĝ − 4ρI)x ≤ xᵀG*x = E(x) for all real x, because |xᵀ(G* − Ĝ)x| ≤ ρ(Σ|x_i|)² ≤ 4ρΣx_i².
  - So the exact Fincke–Pohst enumeration of xᵀ(Ĝ − 4ρI)x ≤ T over Q contains every x with E(x) ≤ T. Each float-guided
    range of the enumeration is certified exactly: both endpoints fail the test, on either side of the centre. 01
    prints the number of ranges of each enumeration (7214 for the largest).
- The elements of norm exactly c are then merged into orbits of ⟨−1, ε₁, ε₂⟩. Two elements are merged only when their
  quotient is verified exactly to be ±ε₁^aε₂^b.
- The largest enumeration is class (1,1), c = 36: 213798 nonzero lattice points with xᵀ(Ĝ − 4ρI)x ≤ 2258.
  [verify_coverage.py](scripts/verify_coverage.py), with a different lower form of E, finds the same counts at every
  bound of 01.

**3.4 Proposition (coverage).** Every z ∈ O with N(z) = c, c a right side of its class, is ±ε₁^aε₂^b·w for one of
17 seeds w:

| class | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|
| (1,0) | (−1,0,0,0) | (2,−1,0,0) | — | — |
| (0,1) | (−1,0,0,0) | none: no element of norm 4 | — | — |
| (1,1) | (−1,0,0,0) | (−3,9,1,−4) | (−1,3,−1,−3), (1,0,1,−1), (−4,1,0,−1), (2,5,−2,−5), (−9,10,2,−4), (−2,2,0,−1) | (−4,3,0,−2), (−2,7,0,−4), (0,11,0,−6), (0,3,0,−2), (−8,5,0,−2), (2,−1,0,0) |

*Proof.* Lemma 3.2 with the certificates of 3.3, and the enumeration of 3.3. ∎

*Remark (class (0,1) has no element of norm 4).* The enumeration shows it. An independent reason: O = Z[β] and
g ≡ t⁴ + t + 1 (mod 2) is irreducible, so O/2O is the field F₁₆. The norm of x reduces mod 2 to the determinant of
multiplication by x̄ on F₁₆, which vanishes only for x̄ = 0. So an even norm forces x ∈ 2O and 16 | N(x).

## 4. Exclusion: the zeros of the markers on each orbit
Script: [02_exclusion.py](scripts/02_exclusion.py) → [.out.txt](scripts/02_exclusion.out.txt), with the helpers in
[lifting2.py](scripts/lifting2.py) (a verbatim copy of the template's). Write ε^{(a,b)} := ε₁^aε₂^b. Lemmas 4.1–4.4 are
those of [triangular_tetrahedral/proof.md](../triangular_tetrahedral/proof.md) §4, repeated here; the basis remark of
4.3 is sharpened, and 4.4 has a paragraph on the steps of 02.

Fix a seed w.
- *−1 ∉ ⟨ε₁, ε₂⟩.* This follows from Remark 3.2(i). 02 also checks it directly: −1 is not a residue of ⟨ε₁, ε₂⟩
  modulo 8 for the classes (1,0) and (1,1), and modulo 4 for class (0,1).
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
  ν(J) ≤ e − 1 − ι + λ < e + λ, so a positive verdict of 4.3 does not depend on the basis. Here λ = 0 at every zero used
  except (8, (1,2)) in class (1,0), where λ = 1 and ν(J) = 2 < e = 3. [verify_exclusion.py](scripts/verify_exclusion.py)
  recomputes λ and ν(J) in three bases ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), part 1); they
  agree whenever they are < e.

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

**4.5 The certificates.** Each known zero has a modulus where 4.3 holds. "Half" marks a solution of G that gives no
solution of F_cl (2.3); it is certified like any other zero.

| class | c | known zeros (m,n) of G | moduli used | Prop. 17 instances |
|---|---|---|---|---|
| (1,0) | 1 | (0,1) | 8, 16 | q=8: λ=0, ν(J)=1 |
| (1,0) | 4 | (1,2) | 8, 16 | q=8: λ=1, ν(J)=2, δ=1 |
| (0,1) | 1 | (0,1), (1,0) half, (2,1) | 3, 9, 27, 7, 19 | q=3: λ=0, ν(J)=0 at all three; q=7 at (0,1), (1,0) and q=19 at (2,1) (siblings): λ=0, ν(J)=0 |
| (1,1) | 1 | (0,1), (3,1) | 3, 9, 5, 4 | q=5 at (0,1): λ=0, ν(J)=0; q=9 at (3,1): λ=0, ν(J)=1 |
| (1,1) | 4 | none | 5 | (no candidates mod 5) |
| (1,1) | 9 | none (6 orbits) | 8, 3, 3, 3, 3, 9 (in seed order) | (no candidates) |
| (1,1) | 36 | (1,2) in the orbit of (2,−1,0,0); none in the other five | 3, 3, 9, 3, 3, 9 (in seed order) | q=9 at (1,2): λ=0, ν(J)=1 |

- **What the output prints.** For each modulus used, |G_q| = [Z² : Λ_q] and the basis (h,0), (−j,k) of 4.1; after
  each step, the Hermite form of the current lattice L and the candidate counts; every certified class with its
  lattice; every sibling class with its witness; and every Proposition 17 instance used. The largest lattice index
  reached is 10,554,462 (class (0,1), c = 1).
- **Phantoms.** In the orbit (0,1), c = 1, −1 ≡ ε₁^{81} (mod 27), so the three classes of −z, z a known zero, are
  candidates modulo 3, 9 and 27. They are the three classes left after q = 27. Two of their subclasses die as siblings
  at q = 7 and three at q = 19, moduli at which −1 is not a residue. In the orbit (1,1), c = 1, −1 is a residue modulo 5
  and 4; the classes of −z there are likewise all accounted for
  ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), part 1: the lines on −1 modulo the used moduli, the
  PHANTOM lines and the sibling-empty classes).
- **Unnecessary moduli.** In that orbit 3 and 9 add nothing once 27 is used (Proposition 17 holds at 27 as well), so
  [27, 7, 19] suffices; in the orbit (1,1), c = 1, 3 adds nothing once 9 is used
  ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), part 2: the minimal sub-lists [27, 7, 19] and [9, 5,
  4]). The certificate is left as the greedy search found it.
- **The same-modulus sibling rule is needed.** In class (1,1), ε₁² = 1 − 3β, so the zeros (0,1) and (3,1) of the orbit c
  = 1 lie in one Λ₃-class, and Proposition 17 fails at q = 3 for both
  ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), part 0 for ε₁², part 1 for q = 3). A sibling rule that
  let Proposition 17 at another modulus vouch for a Λ₃-class would accept an incomplete certificate if a zero were
  dropped (same output, control (a6)).
- **Margins.** The proof of Proposition 17 uses e + λ < 2e − ι and e + δ < 2e − ι. All eleven instances have
  δ = e − 1 − ι, so the second holds with margin exactly 1. The first holds with margin 1 at the eight instances with
  λ = e − 1 − ι (the seven with e = 1, and (8, (1,2))), and with margin 2 at (8, (0,1)), (9, (3,1)) and (9, (1,2)):
  - at q = 8 (e = 3, ι = 1): 3 + δ = 4 < 5, and 3 + λ ≤ 4 < 5 (λ = 1 at (8, (1,2)), at its bound);
  - at q = 9 (e = 2, ι = 0): 2 + δ = 3 < 4, and 2 + λ = 2 < 4;
  - at q = 3, 5, 7, 19 (e = 1): 1 < 2.

**Proof of Theorem A′.**
- By 2.2, a solution of G(m,n) = c is an element z ∈ O with N(z) = c and vanishing markers.
- By 3.4, z or −z is ε^{(a,b)}w for one of the 17 seeds w, with (a,b) ∈ Z_w.
- By 4.5 and 4.4, Z_w = Z₀, the exponents of the known solutions.
- So z is ± a known solution, and (m,n) is read off from z = n − mβ. ∎

**Proof of Theorem A.** Theorem A′ and Lemma 2.3. ∎

## 5. Reconstruction
Script: [03_reconstruction.py](scripts/03_reconstruction.py) → [.out.txt](scripts/03_reconstruction.out.txt).

- *Points of the non-trivial classes.* Each solution of Theorem A has g = √c (1.3) and gives the point
  (X₄/g², ±Y₆/g³), and ±(m,n) ↔ (X, ±Y) is one-to-one:

  | class | c | ±(m,n) of F_cl | g | x | \|y\| |
  |---|---|---|---|---|---|
  | (1,0) | 1 | (0,1) | 1 | −2 | 3 |
  | (1,0) | 4 | (1,2) | 2 | 43 | 282 |
  | (0,1) | 1 | (0,1) | 1 | 4 | 9 |
  | (0,1) | 1 | (1,1) | 1 | 52 | 375 |
  | (1,1) | 4 | (0,1) | 2 | 2 | 5 |
  | (1,1) | 4 | (6,1) | 2 | 5234 | 378661 |
  | (1,1) | 9 | (1,1) | 3 | −1 | 4 |

- *The extra solution of the field form.* The solution (1,0) of G in class (0,1) has m odd. It corresponds to the
  projective point (m:n) = (1:0) of F_cl, that is to (1/2, 0), with F_cl(1,0) = 16 outside {1, 4}, and gives only the
  rational point (1/4, ±33/8).
- *The trivial class* gives x = 8 (1.3).
- *Theorem B.* By (H), now Theorem H of [hypothesis_H.md](hypothesis_H.md), the integral points are the 16 points with x ∈ {−2, −1, 2, 4, 8, 43, 52, 5234}. They agree
  with LMFDB's list and with the Mordell–Weil decompositions of the descent file, and the class of each point is
  (i mod 2, j mod 2) for P = iG₁ + jG₂.
- *Mordell's form.* x² − 17 = y³ is y² = x³ + 17 with the variables exchanged, so its integer solutions are
  y ∈ {−2, −1, 2, 4, 8, 43, 52, 5234} with x = ±3, ±4, ±5, ±9, ±23, ±282, ±375, ±378661.

## 6. Status, comparison and novelty

**6.1 Status.**
- Theorems A′ and A are proved here in full, by finite elementary computations. They were verified in-session by
  independent scripts and refereed in-session by agents ([review.md](review.md), rounds 1 and 2); there is no outside
  audit.
- Theorem B is proved in full: §1 with Theorem H of [hypothesis_H.md](hypothesis_H.md), which replaces the rank, the
  facts N1, N2 and Cassels's theorem, and §5. hypothesis_H.md, Lemma 1.3′ and the other changes for it were verified
  by independent scripts and refereed in-session too (review.md, round 3); there is no outside audit.
- *(H) is a published theorem.* Cassels (Acta Math. 82, 1950, §22 and Tables 1–2) shows, for Y² = X³ − D with
  D = −17, that x − δ (δ³ = D) lies in one of four classes modulo squares; they match our four
  ([literature_check.md](literature_check.md) §3, checked in [lit_checks_17.out.txt](scripts/lit_checks_17.out.txt)).
  His proof (§22, by the argument of §§11–12) uses ideal factorization in Q(∛17), the class number h = 1, the unit
  theorem and the unit η = 18 + 7θ. His Lemma 6, which rests on class field theory, serves other D in the paper (for
  example D = +17, p. 262), not the list at D = −17 (pages 249–251 and 261–262, checked in round 2 of the review).
  That input is not elementary in the project's sense. hypothesis_H.md replaces it: the same factorization, with the
  class number replaced by a Dedekind–Hasse division and the unit by a one-unit descent.

**6.2 What R1 needed, and what it did not.**
- *Coverage* needed no unit theorem: two units per field, and an exact enumeration of at most 213798 lattice points.
- *Exclusion* needed only Proposition 17 with a congruence sieve over q ≤ 27, once with λ = 1. The frame criterion
  and saturation were not needed.
- *New compared with 1611.a1.*
  - Two class forms need field forms: α/2 is integral in class (0,1), and the form of class (1,1) is not monic. The
    substitution lemmas of 2.3 return to the class forms; for class (1,1) the correspondence changes the right side by
    a factor 4.
  - Phantom classes of −z survive three moduli (3, 9 and 27) and die as siblings at 7 and 19 (4.5); phantoms and
    the same-modulus sibling rule already occurred for 1611.a1.
  - An explicit case (class (1,1), c = 1, where ε₁² = 1 − 3β) shows that the sibling rule must use the same modulus:
    a cross-modulus rule would accept an incomplete certificate (4.5).

**6.3 Prior proofs** ([literature_check.md](literature_check.md) §§1–2).
- *The theorem is classical.* The first complete solution is Nagell's (Avh. Norske Vid.-Akad. Oslo 1930, no. 7),
  read through its review only. Hemer (thesis 1952; Ark. Mat. 1954) and London–Finkelstein (1973) solved all
  0 < k ≤ 100 or |k| ≤ 100; Gebel–Pethő–Zimmer (1998) all |k| ≤ 10⁴ by elliptic logarithms; Bennett–Ghadermarzi (2015)
  all 0 < |k| ≤ 10⁷ by cubic Thue equations and Baker's method.
- *The algebraic proofs whose method is known take the cubic route* (Hemer by its review; Nagell by inference from
  Gauthier–Lê and Hilliker): Mordell's correspondence with binary cubic forms of discriminant −108·17, or the
  factorization ±y + √17 = ε(a + b√17)³ in Q(√17). Under Mordell's correspondence the 16 points fall into four cubic
  Thue equations, in four cubic fields of discriminants −1836, −459, −1836 and −204, each with one unit: one exponent
  per equation (computed here, [literature_check.md](literature_check.md) Appendix A). Hemer's finite congruence lemma
  (1954; 1956, Lemma 3) excludes single exponents; it is the one-unit analogue of Proposition 17.
- *Related card.* [mordell_negative_k](../mordell_negative_k/README.md) (LEAN) proves that x² + k = y³ has finitely
  many solutions for each k < 0, so for x² − 17 = y³, by the one-unit cubic route (`ghost.unit` +
  `lifting.binomial`). It does not give the list of solutions.
- *R1 takes the quartic route*: the 2-descent over K = Q(θ), θ³ = −17, eight quartic Thue equations in three S4
  fields with cubic resolvent field K, two units, two exponents. The two routes share no field and no equation.

**6.4 Novelty.**
- In the sources checked, no proof of this theorem by the quartic route, or by any two-exponent Skolem-type (unit)
  argument, was found;
  nor for x² − 17 = y³, nor for a family of Mordell equations containing k = 17.
- **No novelty is claimed yet.** Primary sources that could contain such a proof were not read: Nagell 1930, Hemer
  1952 and 1954, London–Finkelstein 1973 and Ljunggren (Acta Arith. 8, 1963). The cubic route and R1 also share an
  outline (descent, Thue equations, units, congruences in the exponents); whether that makes R1 a rewording is the
  owner's decision.
- *What is not new:* the quartic route itself (the 2-descent over Q(∛k) to quartic Thue equations), named by
  Bennett–Ghadermarzi and implemented for general curves in KANT with Baker-type solving; the two-index mechanism
  (Skolem 1933–35, Ljunggren 1942; on the sister curves, Avanesov 1966/67 and Uchiyama 1984, completion-free); the
  one-exponent congruence lemma for Mordell's equation (Delone; Hemer); the descent and its four classes (Cassels
  1950); the coverage lemma (Koecher's reduction theory).

**6.5 What remains open.**
- *Done (2026-10-15, reviewed 2026-10-26):* an elementary proof of (H) for integral points, in
  [hypothesis_H.md](hypothesis_H.md): Z + Zθ + Zω is not norm-Euclidean (Euclidean minimum at least 1115/1028, the
  value Lemmermeyer and Cavallar–Lemmermeyer record), but a division with multipliers {1, 2}, certified by 795 exact
  boxes, gives the divisibility theory that class number 1 provided (Bézout, prime elements, square allocation); the
  units are ±ηᵏ by a one-unit descent, and the norm condition alone leaves the four classes.
- A complete elementary proof by the cubic route (four one-unit cubic Thue equations), as an independent second
  proof.
- The primary sources listed in 6.4.
- The exclusion sieve has no a priori termination guarantee.

## 7. Files
- [core.py](scripts/core.py): fields, orders, exact ring arithmetic, rigorous intervals, units, descent maps, the
  maps between field forms and class forms.
- [00_fields.py](scripts/00_fields.py): fields, maximal orders, markers, the substitution lemma, units, known
  solutions.
- [01_coverage.py](scripts/01_coverage.py): Lemma 3.2 certificates, enumeration, orbits, unit group.
- [lifting2.py](scripts/lifting2.py) and [02_exclusion.py](scripts/02_exclusion.py): kernel lattices, coset systems,
  Proposition 17, sieve.
- [03_reconstruction.py](scripts/03_reconstruction.py): back to E and to x² − 17 = y³. The 2-descent itself is
  [weierstrass_2descent.py](../../targets/rank2_signatures_2026-10/weierstrass_2descent.py) (section CURVE 10404.f2 of
  its output).
- [lit_checks_17.py](scripts/lit_checks_17.py) → [.out.txt](scripts/lit_checks_17.out.txt): the computations of the
  literature check (the cubic route, the resolvent fields, Cassels's classes).
- Independent checks: [verify_coverage.py](scripts/verify_coverage.py),
  [verify_exclusion.py](scripts/verify_exclusion.py), [verify_bruteforce.py](scripts/verify_bruteforce.py).
- Hypothesis (H), in [hypothesis_H.md](hypothesis_H.md):
  - [hring.py](scripts/hring.py): the ring Z + Zθ + Zω; [hlib.py](scripts/hlib.py): shared helpers;
  - [04_ring.py](scripts/04_ring.py): the ring of triples;
  - [05_cover.py](scripts/05_cover.py) with [dh_boxes.txt](scripts/dh_boxes.txt): the division with multipliers, and
    the proof that the ring is not norm-Euclidean;
  - [06_units.py](scripts/06_units.py): the units;
  - [07_classes.py](scripts/07_classes.py): the class of X − θ;
  - [cover_search.py](scripts/cover_search.py): how the table was found (not used by the proof);
  - independent checks: [verify_H_units.py](scripts/verify_H_units.py), [verify_H_cover.py](scripts/verify_H_cover.py),
    [verify_H_classes.py](scripts/verify_H_classes.py) (review.md, round 3);
  - [lit_H_lemmermeyer.py](scripts/lit_H_lemmermeyer.py), [lit_H_lmfdb.py](scripts/lit_H_lmfdb.py): the identifications
    with the tables cited in [literature_check.md](literature_check.md) §7 (not used by the proof).
