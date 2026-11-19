# Triangular = tetrahedral (1611.a1) by strategy R1: two-unit ghosts and finite-polynomial exclusion

Written 2026-10-15 by Claude (Opus 5.5) as the test of strategy R1 of
[GENERALIZATIONS.md](../../GENERALIZATIONS.md) §3.2 on the first rank-2 curve whose Thue fields have a complex place.
It was revised after the in-session review recorded in [review.md](review.md). Scripts: [scripts/](scripts/). Each
script has its output beside it as `*.out.txt`. Every check is exact; floating point only seeds some searches, whose
results are then verified exactly. Hypothesis (H) of §1.2, which this proof
first took from the rank, is proved in [hypothesis_H.md](hypothesis_H.md) (added 2026-10-15).

## 0. Statements

**0.1 The curve.** E: Y² = X³ − 144X + 1296. It is LMFDB 1611.a1, y² + y = x³ − 9x + 20, through X = 4x, Y = 8y + 4.
LMFDB gives rank 2, an mwrank value. The proof does not use it (see 0.3).
- C(n,2) = C(m,3) holds exactly when (x, y) = (3(m−1), 9n−5) is a point of 1611.a1.
- The identity is 162·(C(n,2) − C(m,3)) = y² + y − (x³ − 9x + 20), checked in
  [weierstrass_2descent.out.txt](../../targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt).

**0.2 The three Thue forms.** They come from the 2-descent of E (§1). Write F(m,n) = a₀m⁴ + a₁m³n + a₂m²n² + a₃mn³ + n⁴,
so that f(t) := F(1,t) = t⁴ + a₃t³ + a₂t² + a₁t + a₀ is monic.

| class | F(m,n) | f(t) = F(1,t) |
|---|---|---|
| (1,0) | 36m⁴ − 6mn³ + n⁴ | t⁴ − 6t³ + 36 |
| (0,1) | −5m⁴ + 14m³n + 24m²n² + 2mn³ + n⁴ | t⁴ + 2t³ + 24t² + 14t − 5 |
| (1,1) | −8m⁴ + 4m³n + 24m²n² + 4mn³ + n⁴ | t⁴ + 4t³ + 24t² + 4t − 8 |

**Theorem A (the twelve Thue equations).** For each class and each c ∈ {1, 4, 9, 36}, the integer solutions of
F(m,n) = c are exactly ±(m,n) for the following pairs (m,n), where "none" means no solutions:

| class | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|
| (1,0) | (0,1) | (1,2) | none | (1,0), (1,6) |
| (0,1) | (0,1) | (−1,1) | (4,1) | (1,1) |
| (1,1) | (0,1) | none | (−1,1), (2,1) | none |

**Theorem B (integral points).** The following hold. (Hypothesis (H) of 1.2, on which this rests, is Theorem H of
[hypothesis_H.md](hypothesis_H.md).)
- the integral points of E are the 28 points with X ∈ {−15, −12, −8, 0, 4, 9, 12, 24, 40, 48, 108, 252, 420, 4953},
  with both signs of Y;
- C(n,2) = C(m,3) with n, m ≥ 0 only for the values 0, 1, 10, 120, 1540, 7140.

**0.3 What is used.**
- **Theorem A** uses:
  - exact computation with integers and rationals;
  - real numbers only through enclosures in rational intervals, obtained from exact sign changes, with the number of
    real roots fixed by Sturm's theorem;
  - three elementary real-variable facts in Lemma 3.2: positive definiteness of E, concavity of the logarithm of a
    positive affine function (or the root-free AM–GM of Remark 3.2(ii)), and e^y ≥ 1 + y;
  - Proposition 17 of the finite-polynomial lifting lemma
    [lifting.binomial.multi](../../methods/lifting.binomial.multi.md);
  - congruences.

  It uses no completion, no p-adic power series, no ideals, no class numbers, no unit theorem and no linear forms in
  logarithms.
- **Theorem B** adds the 2-descent of §1 and [hypothesis_H.md](hypothesis_H.md).
  - **(H)**, the list of four classes, is Theorem H there. It is proved by:
    - a division with multipliers {1, 2} in Z[λ] (λ = θ/6), certified by 398 exact boxes;
    - the divisibility theory that follows from that division;
    - the units ±εᵏ, by a one-unit descent;
    - a congruence modulo 4.

    Z[λ] is not norm-Euclidean, so a plain Euclidean algorithm would not do. The rank of E(Q) and the Mordell–Weil
    theorem (N5) are no longer used.
  - **N1 and N2** of the descent script (the index–discriminant formula; integral elements form a ring, and the
    characteristic-polynomial criterion) are replaced by hypothesis_H.md (3.6): Z[λ] is integrally closed, and every
    element of it is integral. So O_K = Z[λ], as §1.3 uses (hypothesis_H.md §7).

  The rest of §1 (conics, parametrization, gcd bounds) is elementary.

## 1. The reduction to Thue equations (the 2-descent; prior work)
Source: [weierstrass_2descent.py](../../targets/rank2_signatures_2026-10/weierstrass_2descent.py) →
[.out.txt](../../targets/rank2_signatures_2026-10/weierstrass_2descent.out.txt), section `CURVE 1611.a1`.
A verifier re-derived it independently, with other class representatives and conic points
([verify_weierstrass.out.txt](../../targets/rank2_signatures_2026-10/verify_weierstrass.out.txt)), and §1 is also
re-checked in [03_reconstruction.out.txt](scripts/03_reconstruction.out.txt).

**1.1 The cubic field and the descent map.**
- Put K = Q(θ) with θ³ − 144θ + 1296 = 0. Its signature is (1,1), and O_K = Z + Z·θ/6 + Z·θ²/36 = Z[θ/6]. That this
  ring is the whole ring of integers is [hypothesis_H.md](hypothesis_H.md) (3.6); the descent script's N1 is not used.
- For a point P ≠ O put φ(P) = X(P) − θ in K*/K*². φ is a homomorphism whose kernel contains 2E(Q). The proof is
  the one-line identity (x₁−θ)(x₂−θ)(x₃−θ) = (λθ+ν)² for three collinear points, given in the header of the descent
  script, together with the vertical-line case. The script also checks the identity on 109 pairs.
- E(Q)[2] = 0, because the cubic is irreducible.

**1.2 Hypothesis (H).** Every integral point P of E has φ(P) ∈ {1, 9−θ, −15−θ, −12−θ}.
- These four classes are φ of 0, G₁, G₂ and G₁+G₂, for the LMFDB generators G₁, G₂. They are pairwise distinct
  (non-squares modulo 7 and 23, checked exactly).
- (H) would follow from rank E(Q) ≤ 2, by N5, since then |E(Q)/2E(Q)| ≤ 4. This was the original route.
- **(H) is proved for integral points in [hypothesis_H.md](hypothesis_H.md), Theorem H**, without the rank:
  X − θ = μ₀S² with S ∈ Z[λ] and μ₀ ∈ {1, επ₃ₐπ₃ᵦ, −π₃ₐπ₃꜀, −επ₃ᵦπ₃꜀}, and each μ₀ is one of the four classes times a
  square.
  - For the points of his own model (6 | X), Avanesov obtains the analogue of (H) from class number 1 and a
    fundamental unit of K (§6.3). hypothesis_H.md proves both facts instead of citing them.
  - The lattice step 1.3 uses that Z[λ] is integrally closed (hypothesis_H.md (3.6) and §7), in place of N1 and N2.

**1.3 From a class to a Thue equation.** Let μ be the class and X − θ = μξ².
- Then ξ lies in the lattice Ξ_μ = {ξ : μξ² ∈ O_K}, where O_K = Z[λ] (hypothesis_H.md (3.6)); it is a lattice, with
  bounded denominators, by hypothesis_H.md §7. Its coordinates (u,v,w) satisfy A₂ = 0 (a smooth conic),
  A₁ = −1 and A₀ = X.
- A common factor k of (u,v,w) satisfies k² | 6, so k = 1.
- The conic has the rational point ξ = 1. So (u,v,w) = ±X(m,n)/g for coprime (m,n), by Lemma 1.3′ below. Here X(m,n)
  is the parametrization and g = gcd X(m,n) divides 36. The divisibility g | R = 36 follows from the polynomial
  identities that the descent script checks: α²q − n²R and β²q − m²R are l times integral forms, in its notation.
- Then A₁ = −1 becomes F(m,n) = g², with F the class form of 0.2 and g ∈ {1, 2, 3, 6}.
- X = X₄(m,n)/g² and Y = ±Y₆(m,n)/g³; the polynomials are in [core.py](scripts/core.py), `descent_maps`.
- The trivial class gives n(m³ + 3m²n − mn² + 3n³) = 6. Here n | 6, so (m,n) lies on finitely many lines, and it
  gives X ∈ {−8, 4, 40} (§5).

So every integral point outside the trivial class gives coprime (m,n) with F(m,n) ∈ {1, 4, 9, 36} in one of the three
classes. Theorem A lists all of them.

**1.3′ Lemma (parametrization of a conic).** Let q be a ternary quadratic form with integer coefficients and nonzero
Gram determinant, with polar form B(P, P′) = q(P + P′) − q(P) − q(P′). Let P₀ ∈ Z³ be primitive with q(P₀) = 0, and let
(P₀, e, f) be a basis of Z³. For d = me + nf put X(m,n) = q(d)P₀ − B(P₀,d)d. Then q(X(m,n)) = 0. Every primitive
integer solution P of q(P) = 0 equals ±X(m,n)/g for some coprime (m,n), where g is the gcd of the coordinates of
X(m,n).

*Proof.*
1. *q(X) = 0.* q(X) = q(d)²q(P₀) − q(d)B(P₀,d)² + B(P₀,d)²q(d) = 0.
2. *No isotropic d ≠ 0 in span(e, f) has B(P₀, d) = 0.* Otherwise q would vanish on span(P₀, d), a plane. In a basis
   (P₀, d, w) the Gram matrix would then be [[0,0,x],[0,0,y],[x,y,z]], whose determinant is 0, a contradiction.
3. *Write the solution in the basis.* Write P = aP₀ + d with a ∈ Z and d = me + nf.
4. *If d = 0,* then P = ±P₀.
   - The linear form (m,n) ↦ B(P₀, me + nf) is not identically 0. Otherwise B(P₀, ·) would vanish on all of Q³, since
     B(P₀, P₀) = 2q(P₀) = 0, against the nonzero determinant.
   - So it has a primitive zero (m₀, n₀). By step 2, q(m₀e + n₀f) ≠ 0.
   - So X(m₀, n₀) is a nonzero multiple of P₀, and P = ±X(m₀, n₀)/g.
5. *If d ≠ 0,* then 0 = q(P) = aB(P₀,d) + q(d).
   - By step 2, B(P₀,d) ≠ 0. So B(P₀,d)·P = −X(m,n), and P is proportional to X(m,n).
   - If h = gcd(m,n), then X(m,n) = h²X(m/h, n/h), so P is proportional to X(m′,n′) with (m′,n′) coprime.
   - X(m′,n′) ≠ 0, because its (e, f)-coordinates are −B(P₀,d′)(m′,n′) ≠ 0.
   - P is primitive, so P = ±X(m′,n′)/g. ∎

The descent script applies the lemma to the conic C_μ = (A₂ = 0), with P₀ = ξ = 1 and the bases it prints. This
argument was checked by script, but written nowhere until the referee round of 2026-10-15 (review.md, round 4).

## 2. The Thue fields
Script: [00_fields.py](scripts/00_fields.py) → [.out.txt](scripts/00_fields.out.txt).

**2.1 The fields and the norm form.** For each class, f is irreducible with exactly two real roots (Sturm) and Galois
group S4. Write α for a root and L = Q(α).
- F(m,n) = m⁴f(n/m) = N(n − mα).
- The order O used below has the integral basis ω₀..ω₃ of [core.py](scripts/core.py) (`CLASSES[key]['BASIS']`):

| class | ω₂ | ω₃ | [O : Z[α]] | disc O |
|---|---|---|---|---|
| (1,0) | α²/6 | α³/6 | 36 | −2⁴·3²·179 |
| (0,1) | (α²+4α+1)/6 | (α³+3α+2)/6 | 36 | −2⁴·3²·179 |
| (1,1) | (α²+2α+4)/6 | (α³+4)/12 | 72 | −2²·3²·179 |

  In each case ω₀ = 1 and ω₁ = α. O is the maximal order (p-maximality at 2 and 3), but maximality is never used.
  What is used is that O is a ring containing Z[α] and the units below.

**2.2 Markers.** Write x = (x₀, x₁, x₂, x₃) for the coordinates in the basis ω.
- The α²- and α³-coordinates of x are (x₂, x₃)·D for an invertible 2×2 matrix D: diag(1/6, 1/6), or diag(1/6, 1/12)
  for class (1,1).
- Since 1, α, α², α³ are linearly independent (f is irreducible), x ∈ Q + Qα if and only if x₂ = x₃ = 0, and then
  x = x₀ + x₁α ∈ Z + Zα.
- n − mα has coordinates (n, −m, 0, 0).
- So **the solutions of F(m,n) = c are the elements z ∈ O with N(z) = c and x₂(z) = x₃(z) = 0.** The two coordinates
  x₂, x₃ are the *markers*.

**2.3 Units.** Each field gets two units of norm +1 (in ω-coordinates):

| class | ε₁ | ε₂ |
|---|---|---|
| (1,0) | (0, 1, −1, 0) = α − α²/6 | (−86, −41, −117, 30) |
| (0,1) | (0, −1, 1, 0) | (−21, 36, 16, 14) |
| (1,1) | (−1, 0, 1, 0) | (1, −3, −5, −2) |

- They were found by a lattice search, which is not part of the proof. Their norms are checked exactly.
- They are multiplicatively independent by Remark 3.2(i).
- As a by-product of §3, O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩ in all three orders (Remark 3.2(iii)). The pair is fundamental, but
  that is never used.

## 3. Coverage: finitely many unit orbits (the two-unit descent)
Script: [01_coverage.py](scripts/01_coverage.py) → [.out.txt](scripts/01_coverage.out.txt). This section fills the
slot `ghost.unit.multi_unit` of [ghost.unit](../../methods/ghost.unit.md) (see
[ghost.unit.multi_unit](../../methods/ghost.unit.multi_unit.md)).

**3.1 Sizes.** Let σ₁ < σ₂ be the real embeddings and σ₃ a complex one. For x ∈ O ⊗ R put s₁ = σ₁(x)², s₂ = σ₂(x)²,
s₃ = |σ₃(x)|² and E(x) = s₁ + s₂ + s₃. Then:
- (a) N(x)² = s₁s₂s₃².
- (b) E is a positive definite quadratic form. It is the squared length of x ↦ (σ₁x, σ₂x, Re σ₃x, Im σ₃x), which is
  an R-linear isomorphism O⊗R → R⁴, because f has four distinct roots (a Vandermonde determinant). So only finitely
  many lattice points have E ≤ T.
- (c) For a unit u put f_v(u) = s_v(u), the *size factors* of u. Then s_v(ux) = f_v(u)s_v(x).
- *Rigour.*
  - σ₁ and σ₂ are enclosed in rational intervals of width < 2⁻²⁰⁰ by exact sign changes of f, and Sturm's count
    confirms that f has exactly two real roots.
  - The complex pair α₃, ᾱ₃ enters only through the real factor t² − t₃t + n₃ of f, where
    t₃ = α₃ + ᾱ₃ = −a₃ − σ₁ − σ₂ and n₃ = α₃ᾱ₃ = a₀/(σ₁σ₂).
  - If x ≡ A + Bt mod (t² − t₃t + n₃), then s₃ = A² + t₃AB + n₃B².
  - Every quantity below is an interval with rational endpoints ([core.py](scripts/core.py), `Iv`, `Embeddings`).

**3.2 Lemma (reduced representatives).** Let G be a group of units of norm 1 and U ⊂ G finite. Let f_v⁺(u) ≥ f_v(u)
be rational, and let P⁺ be the polygon
{r ∈ R³ : r ≥ 0, r₁+r₂+r₃ = 1, Σ_v r_v(f_v⁺(u) − 1) ≥ 0 for all u ∈ U}.
Suppose every vertex of P⁺ has all coordinates positive, and let b_min = min r₁r₂r₃² over the vertices. Then every
z ∈ O with N(z) = c ≠ 0 is z = g·w with g ∈ G, N(w) = c and E(w)⁴ ≤ c²/b_min.

*Proof.*
1. By 3.1(b), the set Σ_z = {gz : g ∈ G, E(gz) ≤ E(z)} is finite. It contains z. Let w be an element of Σ_z with least
   E.
2. For u ∈ U, uw ∈ Gz. If E(uw) < E(w), then uw ∈ Σ_z has smaller E than w, which is impossible. So
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
- **(i) Independence.** If every vertex of P⁺ is positive, then the log size vectors ℓ(u) = (log f_v(u))_v, u ∈ U,
  span the plane H = {y : y₁ + y₂ + 2y₃ = 0}. In particular ⟨ε₁, ε₂⟩ ≅ Z², so ε₁, ε₂ are multiplicatively
  independent and −1 ∉ ⟨ε₁, ε₂⟩.
  - Each ℓ(u) lies in H, since Σ d_v log f_v(u) = log N(u)² = 0 with d = (1,1,2).
  - Suppose instead that all ℓ(u) lie on one line Rd ⊂ H, with d = 0 allowed. The plane d^⊥ contains p = (1,1,2),
    because d ∈ H = p^⊥. Choose q ∈ d^⊥ not parallel to p, with q having a negative coordinate (replace q by −q if
    necessary).
  - Move from p along q until a coordinate first vanishes: r = p + tq. Then r ≥ 0, r ≠ 0, r ⊥ d, and some r_v = 0.
    Normalize Σ r_v = 1.
  - For u ∈ U, by e^y − 1 ≥ y and r ≥ 0: Σ r_v(f_v⁺(u) − 1) ≥ Σ r_v(f_v(u) − 1) ≥ Σ r_v ℓ_v(u) = ⟨r, ℓ(u)⟩ = 0.
  - So r ∈ P⁺ and has a zero coordinate, against step 4. Hence the ℓ(u) span H.
  - A relation ε₁^aε₂^b = ±1 with (a,b) ≠ (0,0) would put every ℓ(u), u ∈ U ⊂ ⟨ε₁,ε₂⟩, on the line Rℓ(ε₁) or
    Rℓ(ε₂), which is impossible.
- **(ii) Variants.**
  - A weaker bound needs only linearity: s_v(w) ≥ κ_vE(w), with κ_v the least v-coordinate of a vertex, gives
    E(w)⁴ ≤ c²/(κ₁κ₂κ₃²). [verify_coverage.py](scripts/verify_coverage.py) uses it as an independent cross-check.
  - Step 4 also has a proof without logarithms: expand the product binomially along a segment, and apply AM–GM in
    the root-free form "y₁⋯y_N ≥ μ^N implies Σy_i ≥ Nμ", proved by induction.
- **(iii) The unit group.** Let T(1) be the bound for c = 1.
  - Take c = 1 and c = −1. The bound depends only on c², so the lemma puts every unit u as u = gw with E(w) ≤ T(1).
  - In each field the elements of norm 1 with E ≤ T(1) are all ±ε₁^aε₂^b, and none has norm −1
    ([01_coverage.out.txt](scripts/01_coverage.out.txt)).
  - With (i), O^× = {±1}×⟨ε₁⟩×⟨ε₂⟩. So the pair is fundamental, proved without the unit theorem. Nothing below
    uses this.
- **(iv) Provenance.** The step "take an element of a unit orbit with minimal trace form" is Koecher's reduction
  theory of unary forms. In it, the minimizers form a polyhedral cone cut out by finitely many units.
  - Koecher's paper (Math. Ann. 141, 1960) was not accessed; the attribution follows Porter–Bali–Leibak
    (arXiv:2302.00505), with Leibak–Porter–Ling (JTNB 35, 2023).
  - Bounds for a reduced associate are classical: Levesque–Waldschmidt (arXiv:1505.06653), Lemma 6, after
    Shorey–Tijdeman, Lemma A.15; Fieker's thesis (1997), Algorithm 4.8.
  - Lemma 3.2 is an explicit certificate form of this, using the E-minimal element (weights 1 rather than the T₂
    weights 1, 1, 2), a relaxed finite U and the vertex-product bound.
  - The literature check did not find that form, and no novelty is claimed for it.

**3.3 The certificates.** Take G = ⟨ε₁, ε₂⟩ and U = {ε₁^aε₂^b : a, b ∈ {−1,0,1}, (a,b) ≠ (0,0)}. Each f_v⁺ is
the upper endpoint of an interval of width < 10⁻⁵⁵.

| class | vertices of P⁺ | b_min (rounded) | E(w) ≤ (c²/b_min)^{1/4} | T for c = 1, 4, 9, 36 |
|---|---|---|---|---|
| (1,0) | 4, all positive | 7.68·10⁻⁷ | 33.8·√c | 34, 68, 102, 203 |
| (0,1) | 6, all positive | 2.01·10⁻⁶ | 26.6·√c | 27, 54, 80, 160 |
| (1,1) | 6, all positive | 3.81·10⁻⁵ | 12.7·√c | 13, 26, 39, 77 |

- T is the least integer with T⁴ ≥ c²/b_min.
- *Enumeration.*
  - Let G* be the exact Gram matrix of E in the basis ω, Ĝ the symmetric rational matrix formed from the upper
    triangle of the midpoints of its interval enclosure, and ρ the largest radius.
  - Every pivot of Ĝ − 4ρI is positive (checked exactly), so it is positive definite.
  - xᵀ(Ĝ − 4ρI)x ≤ xᵀG*x = E(x) for all real x, because |xᵀ(G* − Ĝ)x| ≤ ρ(Σ|x_i|)² ≤ 4ρΣx_i².
  - So the exact Fincke–Pohst enumeration of xᵀ(Ĝ − 4ρI)x ≤ T over Q contains every x with E(x) ≤ T. Each
    float-guided range of the enumeration is certified exactly: both endpoints fail the test, on either side of the
    centre.
- The elements of norm exactly c are then merged into orbits of ⟨−1, ε₁, ε₂⟩. Two elements are merged only when their
  quotient is verified exactly to be ±ε₁^aε₂^b.

**3.4 Proposition (coverage).** Every z ∈ O with N(z) = c ∈ {1,4,9,36} is ±ε₁^aε₂^b·w for one of 16 seeds w:

| class | c = 1 | c = 4 | c = 9 | c = 36 |
|---|---|---|---|---|
| (1,0) | (−1,0,0,0) | (2,−1,0,0) | (−3,0,6,−1) | (0,0,6,−1) |
| (0,1) | (−1,0,0,0) | (−1,−1,0,0) | (−3,−2,−3,−1) | (1,−1,0,0) |
| (1,1) | (−1,0,0,0) | (3,−1,−2,−1), (2,−1,−1,−1), (0,0,−1,0) | (−1,−1,0,0) | (6,−2,−3,−2), (−1,−2,−3,−1), (0,−1,−3,−1) |

*Proof.* Lemma 3.2 with the certificates of 3.3, and the enumeration of 3.3. ∎

The class (1,1) orbits with c = 4 and c = 36 also die at once: F ≡ n⁴ (mod 4), so an even value of F forces 2 | n,
and then 8 | F. They are treated uniformly below anyway.

## 4. Exclusion: the zeros of the markers on each orbit
Script: [02_exclusion.py](scripts/02_exclusion.py) → [.out.txt](scripts/02_exclusion.out.txt), with the helpers in
[lifting2.py](scripts/lifting2.py). Write ε^{(a,b)} := ε₁^aε₂^b.

Fix a seed w.
- *−1 ∉ ⟨ε₁, ε₂⟩.* This follows from Remark 3.2(i). 02 also checks it directly: −1 is not a residue of ⟨ε₁, ε₂⟩
  modulo 4 for the classes (1,0) and (0,1), and modulo 8 for (1,1). Moduli 19, 11 and 47 give an independent
  cross-check ([verify_bruteforce.out.txt](scripts/verify_bruteforce.out.txt)).
- *Injectivity.* (a,b) ↦ ε^{(a,b)}w is injective, by Remark 3.2(i).
- So each solution pair ±z of the orbit has exactly one element in the family ε^{(a,b)}w, (a,b) ∈ Z², and the
  solution pairs ±z of the orbit correspond one-to-one to Z_w := {(a,b) ∈ Z² : x₂(ε^{(a,b)}w) = x₃(ε^{(a,b)}w) = 0}.
- Theorem A needs only the inclusion "every solution is ±ε^{(a,b)}w with (a,b) ∈ Z_w". Together with Z_w = Z₀
  below, that accounts for every solution even without injectivity.

**4.1 Lemma (kernel lattice).** Let q ≥ 3 and Λ_q = {(a,b) : ε^{(a,b)} ≡ 1 (mod qO)}. Let:
- h be the order of ε₁ mod q;
- k the least integer > 0 with ε₂^k ∈ ⟨ε₁⟩ mod q;
- j an integer with ε₂^k ≡ ε₁^j (mod q).

Then Λ_q = Z(h,0) + Z(−j,k).

*Proof.* The b with ε₂^b ∈ ⟨ε₁ mod q⟩ form a subgroup of Z, namely kZ.
- If (a,b) ∈ Λ_q, then b = sk for some s, so ε₁^{a+sj} ≡ 1 and h | a + sj. Hence (a,b) = ((a+sj)/h)·(h,0) + s·(−j,k).
- Conversely, both generators lie in Λ_q. ∎

**4.2 Classes and coset systems.**
- (a) *Residues.* On an exponent class (a₀,b₀) + Λ_q every element ε^{(a,b)}w is ≡ ε^{(a₀,b₀)}w (mod qO). So if
  either marker of ε^{(a₀,b₀)}w is ≢ 0 (mod q), the class contains no point of Z_w. The other classes are the
  *candidates*.
- (b) *The coset system.* Put η₁ = ε₁^h and η₂ = ε₁^{−j}ε₂^k. They are units of norm 1 with η_j ≡ 1 (mod qO). So their
  multiplication matrices are T_j = I + qK_j: integer matrices that commute and have determinant 1.
  - For z₀ = ε^{(a₀,b₀)}w, the map (i,j') ↦ (a₀,b₀) + i(h,0) + j'(−j,k) is a bijection from Z² onto the exponent class
    (a₀,b₀) + Λ_q, by 4.1.
  - The corresponding elements are T₁^iT₂^{j'}z₀.
  - The markers along the class are f_k(i,j') = l_kT₁^iT₂^{j'}z₀, with l₁ = (0,0,1,0) and l₂ = (0,0,0,1).
- This is the setting of [lifting.binomial.multi](../../methods/lifting.binomial.multi/proof.md), §1 and Theorem 10,
  with q = p^e, and e ≥ 2 when p = 2.

**4.3 Proposition (Jacobian criterion; Proposition 17 of lifting.binomial.multi).** Let z₀ be a common zero. Put:
- x_k = l_kK₁z₀ and y_k = l_kK₂z₀;
- λ = min ν_p(x₁,x₂,y₁,y₂);
- J = x₁y₂ − x₂y₁ and δ = ν_p(J) − λ;
- ι = [p = 2].

If λ ≤ e−1−ι and δ ≤ e−1−ι, then z₀ is the only common zero in Z². That is, it is the only point of Z_w in its class
mod Λ_q.
- The proof uses only finite integer polynomials: the cutoff is chosen after the zeros, and a colength-one argument
  modulo p does the rest.
- The conclusion is about the class, not the basis of Λ_q, so the hypotheses need to hold in one basis only. The
  basis used is that of 4.1.
- *On the basis.* A change of basis replaces K_j by a GL₂(Z)-combination of K₁, K₂ plus q times an integer matrix.
  So (x, y) changes by GL₂(Z) modulo q.
  - Hence λ is unchanged when λ < e, and ν(J) is unchanged when ν(J) < e.
  - At every zero used here λ = 0, so Proposition 17 needs ν(J) ≤ e − 1 − ι < e, and basis-dependent values of ν(J)
    never certify anything.
  - [verify_exclusion.py](scripts/verify_exclusion.py) recomputes λ and ν(J) in three bases. They agree at every
    zero used.

**4.4 Lemma (sieve).** Let q₁, …, q_k be prime powers other than 2 (as in 4.2), L = ∩Λ_{q_i}, and Z₀ ⊂ Z_w the set
of exponents of the known solutions. Every point of Z_w lies in a class v mod L whose image mod each Λ_{q_i} is a
candidate.
- (a) If v contains a point z ∈ Z₀, and 4.3 holds at z for some q_i, then v ∩ Z_w = {z}.
- (b) If v contains no point of Z₀, but lies inside the Λ_{q_i}-class of some z ∈ Z₀ at which 4.3 holds for that
  same q_i, then v ∩ Z_w = ∅.
- So if every candidate class mod L is of kind (a) or (b), then Z_w = Z₀.

*Proof.* The first sentence is 4.2(a). For (a) and (b), the Λ_{q_i}-class of z contains v, and its only point of Z_w
is z, by 4.3. ∎

**4.5 The certificates.** The script finds, for each orbit, a short sequence of moduli for which 4.4 applies. Each
known zero has its own modulus where 4.3 holds:

| class | c | known zeros (m,n) | moduli used | Prop. 17 instances (λ = 0 in all) |
|---|---|---|---|---|
| (1,0) | 1 | (0,1) | 4 | q=4: ν(J)=0 |
| (1,0) | 4 | (1,2) | 4, 31, 7 | q=31: ν(J)=0 |
| (1,0) | 9 | none | 4 | (no candidates mod 4) |
| (1,0) | 36 | (1,0), (1,6) | 4, 31 | q=31: ν(J)=0, 0 |
| (0,1) | 1 | (0,1) | 4, 7, 3 | q=7: ν(J)=0 |
| (0,1) | 4 | (−1,1) | 4, 31, 8 | q=31: ν(J)=0 |
| (0,1) | 9 | (4,1) | 4, 8, 31, 3 | q=31: ν(J)=0 |
| (0,1) | 36 | (1,1) | 4, 31, 8 | q=31: ν(J)=0 |
| (1,1) | 1 | (0,1) | 9, 3*, 27, 13 | q=9: ν(J)=0; q=13 (sibling): ν(J)=0 |
| (1,1) | 4 | none (3 orbits) | 4 | (no candidates mod 4) |
| (1,1) | 9 | (−1,1), (2,1) | 4, 8, 16, 41, 25 | q=16: ν(J)=2 (e=4) for (−1,1); q=41: ν(J)=0 for (2,1); q=25 (siblings): ν(J)=0, 1 (e=2) |
| (1,1) | 36 | none (3 orbits) | 4 | (no candidates mod 4) |

(*) After q = 9, the modulus 3 adds nothing: Λ₉ ⊆ Λ₃, and the candidate test mod 3 is implied by the one mod 9. The
greedy search picked it on a tie.

- **What the output prints.** The group orders |G_q| = [Z² : Λ_q], the candidate counts at each step, and every
  Proposition 17 instance used, whether for certification or for a sibling class. The largest lattice index reached
  is 4,166,400, for class (1,1), c = 9.
- **Margins.** Every instance has λ = 0, and 12 of the 14 have δ = e − 1 − ι exactly.
  - The proof of Proposition 17 uses the strict inequalities e + δ < Λ(2) = 2e − ι and e + λ < Λ(2).
  - For e = 1 (q = 7, 13, 31, 41) these read 1 < 2, and for q = 4 they read 2 < 3.
  - The two instances with δ > 0 give 6 < 7 at (16, (−1,1)) and 3 < 4 at (25, (2,1)).
  - Only (9, (0,1)) and (25, (−1,1)) have slack.
  - The instance (16, (−1,1)) could be dropped. The instance (25, (2,1)) is needed for two sibling classes.

**Proof of Theorem A.**
- By 2.2, a solution of F(m,n) = c is an element z ∈ O with N(z) = c and vanishing markers.
- By 3.4, z or −z is ε^{(a,b)}w for one of the 16 seeds w, with (a,b) ∈ Z_w.
- By 4.5 and 4.4, Z_w = Z₀, the exponents of the known solutions.
- So z is ± a known solution, and (m,n) is read off from z = n − mα. ∎

## 5. Reconstruction
Script: [03_reconstruction.py](scripts/03_reconstruction.py) → [.out.txt](scripts/03_reconstruction.out.txt).
- For each class, the polynomial identity Y₆² = X₄³ − 144X₄F² + 1296F³ holds. So every solution with F(m,n) = g²,
  g = gcd(u,v,w), gives a point of E.
- In the three non-trivial classes u² ≡ F (mod 6) identically. So every coprime solution of F = c has g = √c and gives
  an integral point, and ±(m,n) ↔ (X, ±Y) is one-to-one.
- The trivial-class equation n(m³ + 3m²n − mn² + 3n³) = 6 has the solutions ±(1,1), ±(1,−1), ±(3,−1), ±(7,−2).
  The last has g = 1 and gives only the non-integral point (3601/4, ±216071/8). The trivial class therefore gives
  X ∈ {−8, 4, 40}.
- The solutions of Theorem A give the remaining X ∈ {−15, −12, 0, 9, 12, 24, 48, 108, 252, 420, 4953}.
- By (H), now Theorem H of [hypothesis_H.md](hypothesis_H.md), these are all the integral points: Theorem B.
- On the LMFDB model they are 22 points. Through (x,y) = (3(m−1), 9n−5) they give C(n,2) = C(m,3) ∈
  {0, 1, 10, 120, 1540, 7140}.

## 6. Status, comparison and novelty

**6.1 Status.**
- Theorem A is proved here in full, by finite elementary computations. It is refereed in-session only (see
  [review.md](review.md)), with no outside audit.
- Theorem B is proved in full: §1 with Theorem H of [hypothesis_H.md](hypothesis_H.md), which replaces the rank
  and the facts N1, N2. That file is refereed in-session only, too.

**6.2 What R1 needed, and what it did not.**
- *Coverage* needed no unit theorem and no fundamental units: two units, a polygon of eight linear inequalities, and
  a small exact enumeration.
- *Exclusion* needed only the Jacobian criterion (Proposition 17) with a congruence sieve over small prime powers.
  The frame criterion (Theorem 18) and saturation (Theorem 11) were not needed.
- *No spurious p-adic zero blocked the sieve.* The classes where ε^{(a,b)} ≡ −1 (mod q) carry the "phantom" zero −z.
  These phantoms were removed in one of two ways:
  - by a modulus at which −1 is not a residue (q = 4 for the classes (1,0) and (0,1); q = 8 or 16 for (1,1));
  - in the orbit (1,1), c = 1, as a sibling of the certified zero at q = 13. That works because −1 is not a
    residue modulo 27·13, although it is one modulo 27 and modulo 13 separately.

**6.3 Avanesov (1966/67), read in full** ([literature_check.md](literature_check.md)).
- *His reduction.* He works on 6Y² = X³ − 4X + 6. In terms of his (X, Y), the coordinates of E are X_E = 6X and
  Y_E = 36Y. He uses:
  - class number 1 of Q(λ), λ³ − 4λ + 6 = 0 (Delone–Faddeev);
  - its unit 40λ² − 101λ + 95 (computed by Billevich's method);
  - ideal factorization, with Mordell (1923) cited for the finiteness of the common ideal factors.

  This gives eight cases and six Thue equations with right side 1. His reduction covers only the points of his
  model, those with 6 | X_E.
- *His fields and equations.* They are ours ([00_fields.out.txt](scripts/00_fields.out.txt), §6):
  - the fields correspond through η = 6/α (class (1,0)), η = 3/(α+1) (class (1,1)) and η = 3(1+α)/(1−α)
    (class (0,1));
  - each of his equations L_i = 1 is equivalent to exactly one of ours, by an explicit GL₂(Q) substitution:

    | his lemma | L2 | L3 | L4 | L5 | L6 | L7 |
    |---|---|---|---|---|---|---|
    | our (class, c) | ((1,0), 36) | ((1,0), 4) | ((1,1), 9) | ((1,1), 1) | ((0,1), 4) | ((0,1), 36) |

    For example F₍₁,₀₎(N, −6M) = 36·L2(N, M). These six are exactly the (class, c) whose solutions can give points
    with 6 | X_E, so they are exactly what the tetrahedral problem needs. The other six equations of Theorem A carry
    the non-tetrahedral points X_E = 9, −15, 4953, or no solution;
  - his units of the Lemma 2 field generate the same group as ours for class (1,0), modulo ±1. By Remark 3.2(iii),
    his units are therefore fundamental, which confirms his p. 411 claim (made by Billevich's method);
  - his remark (p. 413) that all three fields have discriminant −2⁴·3²·179 is false for the field of Lemmas 4–5:
    its discriminant is −2²·3²·179 = −6444. The slip has no visible consequence, since Lemmas 3–7 are printed without
    proof.
- *His exclusion.* He proves Lemma 2 by Billevich's fundamental units and Skolem's p-adic method at p = 3: two
  classes, each with a determinant ≢ 0 mod 3, which is Proposition 17. Lemmas 3–7 are stated "similarly", without
  proof.
- *So the mechanism is his.* Cubic-field descent, then (2,1) quartic Thue equations, then Skolem's method in two
  indices.
- *What differs here:*
  - the reduction avoids ideals and class numbers. (H) is proved by a division with multipliers {1, 2} in Z[λ],
    certified by a box table (hypothesis_H.md). Avanesov's field is not norm-Euclidean (M ≥ 121/109, and = 121/109 by
    Cavallar–Lemmermeyer's computation), so his class
    number 1 cannot be replaced by a plain Euclidean algorithm;
  - the units are not assumed fundamental. Coverage comes from Lemma 3.2, and fundamentality appears only as a
    by-product (Remark 3.2(iii)) that is never used;
  - Skolem's analytic step is replaced by the finite-polynomial Proposition 17 and a congruence sieve;
  - all twelve Thue equations are proved, where Avanesov printed the proof of one of six. Every computation is exact
    and re-runnable.
- **Other proofs** ([literature_check.md](literature_check.md) §2):
  - Stroeker–de Weger (Math. Comp. 68, 1999) use elliptic logarithms.
  - Gallegos-Ruiz–Katsipis–Tengely–Ulas (J. Number Theory 208, 2020) recompute the case with Magma's integral
    points.
  - No elementary proof was found.
- **The completion-free form is in print too.** Uchiyama (Tsukuba J. Math. 8, 1984, 131–157) solved the sister
  problem, triangular = square-pyramidal (8604.a1). He used eight (2,1) quartic Thue equations and a two-index
  argument written entirely as congruences modulo 2^{h+c}, with no p-adic series, and says his method is alike in
  principle to Avanesov's.
- **No novelty is claimed for the theorem or the mechanism.** By the project's criterion, this is a reconstruction of
  Avanesov's proof in finite form. Even its finite, congruence-based shape has a published precedent on a sister
  curve.
- **What may be new here is presentation only.** It consists of:
  - the explicit cone certificate for coverage;
  - the kernel-lattice sieve with Proposition 17;
  - the complete re-runnable certificates for all twelve equations.

**6.4 What remains open.**
- *Done (2026-10-15):* an elementary proof of (H) for integral points, in [hypothesis_H.md](hypothesis_H.md).
  - The suggested route fails: Z[λ] is not norm-Euclidean, because the point ξ* fixed by ε needs norm 121/109.
  - What works is a division with multipliers {1, 2} (Dedekind–Hasse), certified by 398 exact boxes, with a one-unit
    descent for the units.
  - A second proof of unique factorization uses the box principle.
- The same machinery should apply verbatim to 8604.a1 (triangular = square-pyramidal; Uchiyama 1984 has the
  classical proof) and to 10404.f2 (y² = x³ + 17). Their 2-descents are in the same output file.
- The exclusion sieve has no a priori termination guarantee. A profinite zero of the markers that is not an integer
  zero would block it; none occurred here.

## 7. Files
- [core.py](scripts/core.py): fields, orders, exact ring arithmetic, rigorous intervals, units, descent maps.
- [00_fields.py](scripts/00_fields.py): fields, maximal orders, markers, units, known solutions, Avanesov's fields.
- [01_coverage.py](scripts/01_coverage.py): Lemma 3.2 certificates, enumeration, orbits, unit group.
- [lifting2.py](scripts/lifting2.py) and [02_exclusion.py](scripts/02_exclusion.py): kernel lattices, coset systems,
  Proposition 17, sieve.
- [03_reconstruction.py](scripts/03_reconstruction.py): back to E, to 1611.a1 and to C(n,2) = C(m,3).
- Independent checks: [verify_coverage.py](scripts/verify_coverage.py), [verify_exclusion.py](scripts/verify_exclusion.py),
  [verify_bruteforce.py](scripts/verify_bruteforce.py).
- Hypothesis (H), in [hypothesis_H.md](hypothesis_H.md):
  - [cubic.py](scripts/cubic.py): shared definitions for Z[λ];
  - [04_ring.py](scripts/04_ring.py): the ring of triples;
  - [05_dh_cover.py](scripts/05_dh_cover.py) with [dh_boxes.txt](scripts/dh_boxes.txt): the division with
    multipliers, and the proof that Z[λ] is not norm-Euclidean;
  - [06_units.py](scripts/06_units.py): the units;
  - [07_classes.py](scripts/07_classes.py): the class of X − θ;
  - [08_ufd_box.py](scripts/08_ufd_box.py): unique factorization by the box principle;
  - [cover_search_K.py](scripts/cover_search_K.py): how the table was found (not used by the proof).
