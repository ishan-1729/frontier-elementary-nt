# Hypothesis (H) without the rank: arithmetic of Z[λ], λ³ = λ − 3

Written 2026-10-15 by Claude (Opus 5.5), revised 2026-10-26 after review round 3 ([review.md](review.md)). It proves
hypothesis (H) of [proof.md](proof.md) §1.2 for integral points,
with elementary means only. With it, Theorem B of proof.md no longer depends on the rank of E(Q). The method is that
of [norm.triples](../../methods/norm.triples.md), as in the proof for 1611.a1
([triangular_tetrahedral/hypothesis_H.md](../triangular_tetrahedral/hypothesis_H.md)), but simpler: Z[λ] is
norm-Euclidean, so a plain Euclidean division suffices. Scripts: [scripts/](scripts/), each with its output beside it
as `*.out.txt`. Every check in them is exact. Floating point only seeds some searches, and each result is then verified
in exact arithmetic.

**Theorem H.** Let (X, Y) be an integer solution of

    Y² = X³ − 9X + 81                                                  (E)

and let θ = 3λ, where λ³ = λ − 3, so that θ³ − 9θ + 81 = 0. Then

    X − θ = μ₀ S²   with S ∈ Z[λ] and μ₀ ∈ {1, π₃ₐπ₃꜀, π₃ₐπ₃ᵦ, π₃ᵦπ₃꜀},

where π₃ₐ = λ, π₃ᵦ = λ − 1 and π₃꜀ = λ + 1. Moreover

    3 − θ = π₃ₐπ₃꜀·π₃ᵦ²,   −3 − θ = π₃ₐπ₃ᵦ·π₃꜀²,   −θ = π₃ᵦπ₃꜀·π₃ₐ².

So X − θ ∈ μK*² with μ ∈ {1, 3 − θ, −3 − θ, −θ}, which is hypothesis (H).

**What the proof uses.**
- Integer and rational arithmetic, polynomial identities, congruences, induction (descent), and real inequalities.
- Real numbers enter only through rational interval enclosures of λ₁, the one real root of t³ − t + 3, obtained from
  exact sign changes. They are used for the units (Section 4) and for Remark 2.5, which Theorem H does not use. The
  division of Section 2 needs only rational triples and a table of 116 boxes, checked in exact rational arithmetic.
- No ideal theory, class number, reciprocity law, unit theorem, Mordell–Weil theorem or rank is used.
- Uchiyama (1984) took the class number 1 of Q(λ) without a source, and the fundamental unit from the Delone–Faddeev
  table. Here both are replaced by proofs. Sections 2–3 give the divisibility theory that class number 1 provided
  (Bézout, prime elements, square allocation), and Section 4 gives the unit.

**Where (H) was used, and what replaces it.** Every item below is proved here.

| input of proof.md, Theorem B | before | now |
|---|---|---|
| (H), the list of four classes (§1.2) | rank E(Q) ≤ 2 (LMFDB, mwrank) with Mordell–Weil (N5) | Theorem H, Sections 1–5 |
| N2: integral elements form a ring; α² integral implies α integral (§1.3) | standard facts, not proved | (3.6): Z[λ] is integrally closed |
| N1: O_K = Z[λ], via the index–discriminant formula (§1.1) | squarefree discriminant | (3.6): Z[λ] is the set of algebraic integers of K |

**Checks.** A reference such as [04: 1.5] means section 1.5 of `04_ring.py`, with the result in `04_ring.out.txt`.
The ring is defined in [hring.py](scripts/hring.py), and [hlib.py](scripts/hlib.py) holds helpers shared with the card
[mordell_k_minus17](../mordell_k_minus17/hypothesis_H.md) (intervals, Taylor and Bernstein bounds, lattice enumeration).

| Script | Checks |
|---|---|
| [04_ring.py](scripts/04_ring.py) ([output](scripts/04_ring.out.txt)) | ring laws, norm, adjugate, nonvanishing of the norm, the two factors, special elements, real values (Section 1) |
| [05_cover.py](scripts/05_cover.py) ([output](scripts/05_cover.out.txt)) | the box table of the division, its cover of D, the lower bound 8/9 at a point fixed by ε (Section 2) |
| [06_units.py](scripts/06_units.py) ([output](scripts/06_units.out.txt)) | root interval, descent constants, the units of small height (Section 4) |
| [07_classes.py](scripts/07_classes.py) ([output](scripts/07_classes.out.txt)) | the two factors, the norm condition, the step modulo 4, the identification, the integral points (Section 5) |
| [cover_search.py](scripts/cover_search.py) ([output](scripts/cover_search.out.txt)) | how the box table [euclid_boxes.txt](scripts/euclid_boxes.txt) was found (not used by the proof) |
| [verify_H_units.py](scripts/verify_H_units.py), [verify_H_cover.py](scripts/verify_H_cover.py), [verify_H_classes.py](scripts/verify_H_classes.py) (with outputs) | independent re-implementations of Sections 1 and 4, 2, and 5 (review.md, round 3) |
| [lit_H_lemmermeyer.py](scripts/lit_H_lemmermeyer.py), [lit_H_lmfdb.py](scripts/lit_H_lmfdb.py) (with outputs) | the identifications with the tables of the literature check, §6 (not used by the proof) |

---

## 1. The ring of integer triples

**(1.1) Multiplication.** A triple (a, b, c) stands for a + bλ + cλ², where λ³ = λ − 3, so λ⁴ = λ² − 3λ. Define

    (a,b,c)(d,e,g) = (ad − 3bg − 3ce,  ae + bd + bg + ce − 3cg,  ag + be + cd + cg).

This is the product of a + bt + ct² and d + et + gt², reduced modulo f(t) = t³ − t + 3 [04: 1.1]. With coordinatewise
addition the triples form a commutative ring with unit (1, 0, 0) [04: 1.1]. Write R for the ring of integer triples and
R_Q for the ring of rational triples. An integer r is the triple (r, 0, 0). Put θ = 3λ = (0, 3, 0); then
θ³ − 9θ + 81 = 0, and X − θ is the triple (X, −3, 0).

**(1.2) Norm.** Let M(α) be the matrix whose columns are α·1, α·λ, α·λ². For α = (a, b, c),

    M(α) = [[a, −3c, −3b], [b, a + c, b − 3c], [c, b, a + c]],

and M(αβ) = M(α)M(β) [04: 1.2]. The norm is the cubic form

    N(α) = det M(α) = a³ + 2a²c − ab² + 9abc + ac² − 3b³ + 3bc² + 9c³.

The proof uses only the polynomial identities N(αβ) = N(α)N(β), N(rα) = r³N(α) and N(−α) = −N(α) [04: 1.2]. Also
N(X − θ) = X³ − 9X + 81 identically in X [04: 1.5]. The basis (1, λ, λ²) has discriminant det(Tr(λ^{i+j})) = −239
[04: 1.2].

**(1.3) Adjugate triple.** For α = (a, b, c) put

    α* = (a² + 2ac − b² + 3bc + c²,  −ab − 3c²,  −ac + b² − c²).

Then α·α* = (N(α), 0, 0) [04: 1.3]. Consequences:
- **(1.3a)** α ∈ R is a unit of R exactly when N(α) = ±1.
- **(1.3b)** A rational triple α with N(α) ≠ 0 has the inverse α*/N(α) in R_Q.

**(1.4) Lemma: the norm vanishes only at 0.** If α ∈ R_Q and N(α) = 0, then α = 0.

*Proof.* Suppose α ≠ 0. Scale it to integers (a, b, c) with gcd 1 and N(a, b, c) = 0. Modulo 2, N(a, b, c) is odd for
each of the seven nonzero residue triples [04: 1.4] (it is the norm of F₂[t]/(t³ + t + 1), the field with eight
elements). So a, b and c are all even, against gcd 1. ∎

So R has cancellation. We write β | α when α = βγ with γ ∈ R; then N(β) | N(α). Associates differ by a unit.

**(1.5) Special elements** [04: 1.6]:

| element | triple | norm | identity |
|---|---|---|---|
| ε = λ² + λ − 1 | (−1, 1, 1) | 1 | ε⁻¹ = λ² − 2λ + 2 = (2, −2, 1) |
| π₃ₐ = λ | (0, 1, 0) | −3 | π₃ₐπ₃ᵦπ₃꜀ = λ³ − λ = −3 |
| π₃ᵦ = λ − 1 | (−1, 1, 0) | −3 | |
| π₃꜀ = λ + 1 | (1, 1, 0) | −3 | |
| π₂₃₉ = 3λ² − 1 | (−1, 0, 3) | 239 | F′(θ) = 3θ² − 9 = 9π₂₃₉ = (π₃ₐπ₃ᵦπ₃꜀)²π₂₃₉, where F(x) = x³ − 9x + 81 |

ε is Uchiyama's fundamental unit. No one of π₃ₐ, π₃ᵦ, π₃꜀, π₂₃₉ divides another [04: 1.6]; by hand:
- For a prime p and an integer r with f(r) ≡ 0 (mod p), a + bλ + cλ² ↦ a + br + cr² (mod p) is a ring map R → F_p,
  because the product (1.1) differs from the product of polynomials by a multiple of f.
- The maps with r = 0, 1, −1 (mod 3) send π₃ₐ, π₃ᵦ, π₃꜀ respectively to 0, and the map with r = 124 (mod 239) sends π₂₃₉
  to 0. Each of the four maps sends none of the other three elements to 0 [04: 1.6].
- If x | y, every ring map that sends x to 0 sends y to 0. So no one of the four divides another.

**(1.6) Real values.** f(−2) = −3 < 0 < 9/8 = f(−3/2), so f has a root λ₁ ∈ (−2, −3/2), and λ₁² > 9/4 > 4/3. Bisection
encloses λ₁ in a rational interval of width below 2⁻¹⁹⁹ [04: 1.7].
- f(t) = (t − λ₁)(t² + λ₁t + λ₁² − 1), and t² + λ₁t + λ₁² − 1 = (t + λ₁/2)² + (3λ₁² − 4)/4.
- 3λ₁² − 4 > 0, so the quadratic factor has no real root, and λ₁ is the only real root of f.
- The quadratic factor has the complex pair z, z̄ as roots, with z + z̄ = s := −λ₁ and zz̄ = p := λ₁² − 1.

For α = (a, b, c) put σ₁(α) = a + bλ₁ + cλ₁² and Q(α) = A² + sAB + pB², where A + Bt is α reduced modulo t² − st + p
(A = a − cp, B = b + cs) [04: 1.7]; Q(α) = |a + bz + cz²|². The following hold as polynomial identities modulo
f(λ₁) = 0 [04: 1.7]:
- σ₁(α)Q(α) = N(α);
- σ₁(αβ) = σ₁(α)σ₁(β);
- Q(αβ) = Q(α)Q(β).

*Signs.* Q = (A + sB/2)² + (3λ₁²/4 − 1)B², and 3λ₁²/4 − 1 > 0 since λ₁² > 4/3. So Q ≥ 0 on real triples, where it is
only semidefinite (it vanishes where A = B = 0). Q(α) > 0 for every nonzero rational triple α, since
σ₁(α)Q(α) = N(α) ≠ 0 by (1.4).

*Height.* T₂(α) = σ₁(α)² + 2Q(α) has, modulo f(λ₁), the Gram matrix [[3, 0, 2], [0, 3λ₁² − 2, λ₁], [2, λ₁,
−λ₁² − 9λ₁ + 2]], with leading minors 3, 9λ₁² − 6 ≈ 19.15 and 239 [04: 1.7]. So T₂ is a positive definite quadratic
form, and only finitely many integer triples have T₂ below any given bound.

## 2. Euclidean division

**(2.1) Lemma.** For every rational triple ξ in D = [−½, ½] × [−½, ½] × [0, ½] there is κ ∈ R with |N(ξ − κ)| < 1.

*Proof.* The table [euclid_boxes.txt](scripts/euclid_boxes.txt) lists 116 closed boxes in D. Each has ends that are
multiples of 1/128 and a translation κ ∈ R.

*The bound on a box.* For a box with centre ξ_c and half-widths h, put v = ξ_c − κ and w = ξ − ξ_c, so |wᵢ| ≤ hᵢ. The
exact Taylor expansion of the cubic N [05: 2.1],

    N(v + w) = N(v) + Σᵢ gᵢ(v)wᵢ + Σ_{i≤j} q_{ij}(v)wᵢwⱼ + N(w),

with, for v = (x, y, z),

    g₁ = 3x² + 4xz − y² + 9yz + z²,   g₂ = −2xy + 9xz − 9y² + 3z²,   g₃ = 2x² + 9xy + 2xz + 6yz + 27z²,
    q₁₁ = 3x + 2z,  q₁₂ = −2y + 9z,  q₁₃ = 4x + 9y + 2z,  q₂₂ = −x − 9y,  q₂₃ = 9x + 6z,  q₃₃ = x + 3y + 27z,

gives |N(ξ − κ)| ≤ B on the box, where

    B = |N(v)| + Σᵢ |gᵢ(v)|hᵢ + Σ_{i≤j} |q_{ij}(v)|hᵢhⱼ + Σ_m |n_m|h^m.

Here n_m is the coefficient of N at the exponent triple m, h^m = h₁^{m₁}h₂^{m₂}h₃^{m₃}, and Σ|n_m| = 29. B is computed
exactly, and B < 1 for every box; the largest is 130993/131072 [05: 2.2].

*The boxes cover D.* Each of the 128·128·64 = 1048576 cells of side 1/128 that tile D lies in exactly one box. So every
point of D lies in a closed cell, hence in a box [05: 2.3]. Two more checks confirm the cover:
- the boxes are exactly the leaves of a binary tree of midpoint splits of D, listed in depth-first order;
- their volumes add up to vol(D) = ½.

So the listed κ gives |N(ξ − κ)| ≤ B < 1 at every point of D. The argument works for real ξ as well; only rational ξ are
needed. ∎

**(2.2) Corollary.** For every rational triple ξ ∉ R there is κ ∈ R with 0 < |N(ξ − κ)| < 1.

*Proof.* Let r be the triple of nearest integers to the coordinates of ξ, and put ξ′ = ξ − r, or ξ′ = r − ξ if the third
coordinate of ξ − r is negative. Then ξ′ ∈ D. Take κ′ for ξ′ from (2.1). If ξ′ = ξ − r, then ξ − (κ′ + r) = ξ′ − κ′;
if ξ′ = r − ξ, then ξ − (r − κ′) = −(ξ′ − κ′). Either way we get κ with |N(ξ − κ)| = |N(ξ′ − κ′)| < 1. Since ξ ∉ R,
ξ − κ ≠ 0, and N(ξ − κ) ≠ 0 by (1.4). ∎

**(2.3) Euclidean division.** Let α, β ∈ R with β ≠ 0 and β ∤ α. Then there is κ ∈ R with

    0 < |N(α − κβ)| < |N(β)|.

*Proof.* N(β) ≠ 0 by (1.4). The rational triple ξ = α·β*/N(β) satisfies ξβ = α, and ξ ∉ R because β ∤ α. Take κ from
(2.2) and put ρ = α − κβ = (ξ − κ)β ∈ R. Then |N(ρ)| = |N(ξ − κ)|·|N(β)| and 0 < |N(ξ − κ)| < 1. ∎

**(2.4) Remark: the table.**
- *Shape.* 116 boxes; the deepest is reached by 17 midpoint splits, and every |κᵢ| ≤ 8 [05: 2.2].
- *How it was found.* A search ([cover_search.py](scripts/cover_search.py)) tried, on each box with centre ξ_c, the
  translates κ = ε^{−k}κ′ with κ′ an integer triple near ε^k·ξ_c and k ∈ {−1, 0, 1} (then N(ξ − κ) = N(ε^kξ − κ′)),
  with two split rules, and kept the smaller table. The proof uses only the table.
- *Sharp bound.* The Taylor bounds are not sharp. Expand N(ξ − κ) on each box in the 64 products of cubic Bernstein
  polynomials. The largest coefficient over all boxes is 32733/32768. It is a corner value, at ξ = (−3/8, −7/32, 7/16),
  in the box ZyxYZxzyXyZxY with κ = 0. So the largest value of |N(ξ − κ)| over the boxes, each on its own closed box,
  is exactly 32733/32768 [05: 2.6]. The corner lies in five other boxes, where the values are 22973/32768 (twice),
  26109/32768 (twice) and 24061/32768 [05: 2.6].
- *Confirmation.* 3000 random pairs (α, β), 2998 of them with β ∤ α, were divided by (2.3), with the expected
  inequality every time [05: 2.7].

**(2.5) Remark: the Euclidean minimum is at least 8/9.** Theorem H does not use this remark. Let ξ₀ = (−1/3, −2/9, 4/9).
Then εξ₀ − ξ₀ = (0, −1, −1) and ε⁻¹ξ₀ − ξ₀ = (3, −2, 1) are in R, so the coset C = ξ₀ + R satisfies εC = C. Its element
ξ₀ + 1 = (2/3, −2/9, 4/9) has |N| = 8/9, and no element of C has smaller |N| [05: 2.5].

*Proof of the bound.* Suppose some y ∈ C has |N(y)| < 8/9.
- The orbit {ε^k y : k ∈ Z} lies in C, since εC = C and ε⁻¹C = C. Its elements with T₂ ≤ T₂(y) are finitely many,
  since T₂ is positive definite and C is a translate of the lattice R. So the orbit has an element y₀ of least T₂. Also
  N(y₀) = N(y).
- Neither εy₀ nor ε⁻¹y₀ has smaller T₂, so T₂(y₀)³ ≤ K₃·N(y)², where K₃ = 4(E² + E + 1)³/(E²(E + 1)²) with E = Q(ε)
  as in (4.2), and K₃ ≤ 309.63 [05: 2.5.2].
  - This is the scaled lemma (4.3′) with its exact threshold (see after (4.3′)). Here the rounded C of (4.3′) would
    also do, since C³·(8/9)² ≤ T_e³ [05: 2.5.3].
  - The expression is unchanged under E ↦ 1/E. 05 writes it with e = σ₁(ε) = 1/E, and bounds K₃ above by its value at
    the lower end of the rational enclosure of e.
- So T₂(y₀)³ < K₃·(8/9)² ≤ T_e³, with T_e := 1251/200 = 6.255.
- [05: 2.5.3] uses T_low(y) = yᵀ(G′ − 3ρI)y. Here G′ is a rational matrix within ρ = 2⁻²⁰ of the Gram matrix G of T₂,
  entrywise; that follows from the rational enclosure of λ₁.
  - T_low ≤ T₂: T₂ − T_low = yᵀ(G − G′)y + 3ρΣyᵢ² ≥ −ρ(Σ|yᵢ|)² + 3ρΣyᵢ² ≥ 0, since (Σ|yᵢ|)² ≤ 3Σyᵢ².
  - T_low is positive definite, by an exact LDLᵀ.
  - 05 lists the five elements y ∈ C with T_low(y) ≤ T_e. The list contains every y ∈ C with T₂(y) ≤ T_e. The values of
    |N| on it are 8/9, 1, 1, 1 and 7/3, so its least |N| is 8/9, attained only at ξ₀ + 1.
- y₀ is on that list and has |N(y₀)| < 8/9, a contradiction. ∎

*The other fixed points.* The points fixed by ε modulo R are (ε − 1)⁻¹R; there are |N(ε − 1)| = 9 classes of them
modulo R, and |N(ε + 1)| = 9 classes fixed by −ε [05: 2.5.4].
- On the classes of ±ξ₀ the least |N| is exactly 8/9, by the bound above and N(−y) = −N(y).
- Each of the other seven ε-fixed classes, and each of the nine (−ε)-fixed classes, has an element with |N| ≤ 1/3
  [05: 2.5.4].

Cavallar and Lemmermeyer give the Euclidean minimum of this field as M = 8/9 (§7 below). The table proves
M ≤ 32733/32768 < 1, by the Bernstein coefficients of (2.4) (and M ≤ 130993/131072 by the Taylor bounds). That is all
the proof needs.

## 3. Common divisors, primes and square allocation

**(3.1) Gcd and Bézout.** Let α, β ∈ R, not both 0. Then there is δ ∈ R with δ | α, δ | β and δ = xα + yβ for some
x, y ∈ R. Every common divisor of α and β divides δ.

*Proof.* The set L = {xα + yβ : x, y ∈ R} contains a nonzero element. Choose a nonzero δ ∈ L with |N(δ)| least; |N| is a
positive integer on nonzero elements of R, by (1.4).
- Suppose δ ∤ α. Then (2.3) gives κ with 0 < |N(α − κδ)| < |N(δ)|.
- But α − κδ ∈ L, which contradicts the choice of δ. So δ | α, and likewise δ | β.
- A common divisor of α and β divides xα + yβ = δ. ∎

The rest of this section is that of [triangular_tetrahedral/hypothesis_H.md](../triangular_tetrahedral/hypothesis_H.md)
§3 (itself mordell_1963 §3), with R this ring. Those proofs use only (1.2) (N is multiplicative and integer-valued on
R), (1.3a), (1.3b), (1.4), (3.1) and, in (3.6), the characteristic polynomial [04: 1.2]; no Euclidean division enters
after (3.1). The statements and short proofs are repeated so that this file stands alone.

*Gcd and coprime.* We call δ a gcd of α and β. We call α and β coprime if every common divisor is a unit. This holds
exactly when xα + yβ = 1 for some x, y ∈ R.

**(3.2) Euclid's lemma.** If γ | αβ and γ, α are coprime, then γ | β.

*Proof.* From xγ + yα = 1 we get β = xγβ + yαβ, and γ divides both terms. ∎

**(3.2′) Products.** If γ is coprime to α and to β, then γ is coprime to αβ.

*Proof.* Multiply xγ + yα = 1 and x′γ + y′β = 1: the product is γ(…) + yy′αβ = 1. ∎

**(3.3) Prime elements.** Let ϖ ∈ R with N(ϖ) = ±q for a prime number q. Then ϖ is not a unit; every divisor of ϖ is a
unit or an associate of ϖ; and if ϖ | αβ then ϖ | α or ϖ | β. This applies to π₃ₐ, π₃ᵦ, π₃꜀ (q = 3) and π₂₃₉
(q = 239).

*Proof.* N(ϖ) ≠ ±1, so ϖ is not a unit (1.3a). If ϖ = δγ, then N(δ)N(γ) = ±q, so δ or γ is a unit. Suppose ϖ | αβ and
ϖ ∤ α. A gcd of ϖ and α divides ϖ and is not an associate of ϖ, so it is a unit; ϖ and α are coprime, and ϖ | β by
(3.2). ∎

Two prime elements are either associates or coprime: a gcd divides both, so it is a unit or an associate of each.

**(3.4) Exact powers.** Let ϖ be as in (3.3) and α ∈ R nonzero. If ϖⁿ | α then |N(ϖ)|ⁿ divides |N(α)| ≥ 1. So there is
a largest n ≥ 0 with ϖⁿ | α, and then α = ϖⁿα′ with ϖ ∤ α′.

**(3.5) Lemma: coprime square allocation.** Let α, β, W ∈ R with W ≠ 0, and let u be a unit, with αβ = uW² and α, β
coprime. Then α = u′R₀² for a unit u′ and some R₀ ∈ R.

*Proof.* Induction on |N(W)|. Let g be a gcd of α and W.
- *g a unit.* Then xα + yW = 1. Squaring, and using W² = u⁻¹αβ, α divides 1. So α is a unit, and α = α·1².
- *g not a unit.* Write α = gα₁ and W = gW₁. Cancelling g in gα₁β = ug²W₁² gives α₁β = ugW₁². g and β are coprime (a
  common divisor divides α and β), so g | α₁ by (3.2); write α₁ = gα₂. Then α₂β = uW₁², α₂ and β are coprime, and
  |N(W₁)| < |N(W)|. By induction α₂ = u′R₁², so α = g²α₂ = u′(gR₁)².

When |N(W)| = 1, W is a unit, every gcd of α and W is a unit, and the first case applies. ∎

**(3.6) Remark: integral closedness (replaces N1 and N2 of the descent).** By (1.3b) and (1.4), R_Q is a field; it is
K = Q(λ). R is integrally closed in it.
- *R is integrally closed.* Suppose ξ ∈ R_Q satisfies ξⁿ + c_{n−1}ξ^{n−1} + ⋯ + c₀ = 0 with integers cᵢ.
  - Write ξ = α/m with α ∈ R and m ≥ 1 an integer. Divide α and m by a gcd δ (3.1): ξ = α′/m′ with α′, m′ ∈ R coprime
    (m′ need not be an integer).
  - Multiplying the equation by m′ⁿ gives m′ | α′ⁿ. Raising xα′ + ym′ = 1 to the n-th power shows m′ | 1. So m′ is a
    unit and ξ = α′m′⁻¹ ∈ R.
- *Every element of R is integral over Z.* Each triple α is a root of the characteristic polynomial of M(α), which is
  monic with integer coefficients [04: 1.2] (Cayley–Hamilton, and α ↦ M(α) is an injective ring map).

So the elements of K that are roots of monic integer polynomials are exactly R; that is, O_K = R. In particular, if
ξ ∈ K and ξ² = c ∈ R, then ξ ∈ R: ξ is a root of p_c(t²), where p_c is the monic characteristic polynomial of c.
Section 6 uses this in place of the facts N1 and N2 of the descent script.

## 4. The units

**(4.1) Theorem.** The units of R are exactly ±ε^k (k ∈ Z), and k is determined by the unit. So every unit u can be
written u = σε^l v² with σ ∈ {1, −1}, l ∈ {0, 1} and v a unit.

**(4.2) Values of ε** [06: 4.2].
- e := σ₁(ε) ∈ (0.12288061, 0.12288062), so 0 < e < 1.
- E := Q(ε) = σ₁(ε*) = σ₁(ε⁻¹) = λ₁² − 2λ₁ + 2, enclosed directly in (8.137980, 8.137981).
- σ₁(ε)·E = N(ε) = 1. So σ₁(ε⁻¹) = E and Q(ε⁻¹) = e.

**(4.3) Descent lemma.** For a nonzero rational triple y put x = σ₁(y)² and q = Q(y). Both are positive: σ₁(y) ≠ 0
because σ₁(y)Q(y) = N(y) ≠ 0 by (1.4), and q > 0 by (1.6). Also xq² = N(y)² and T₂(y) = x + 2q. By (1.6) and (4.2)
[06: 4.3]:

    T₂(yε) − T₂(y) = (E − 1)(2q − x(E + 1)/E²),      T₂(yε⁻¹) − T₂(y) = (E − 1)(x(E + 1) − 2q/E).

So T₂(yε) < T₂(y) exactly when x/q > ρ_hi := 2E²/(E + 1), and T₂(yε⁻¹) < T₂(y) exactly when x/q < ρ_lo := 2/(E(E + 1)).
ρ_hi increases and ρ_lo decreases with E, so E < 8.137981 gives ρ_hi < 14.49483 and ρ_lo > 0.02689446 [06: 4.3].

*Lemma.* Let u be a unit with T₂(uε) ≥ T₂(u) and T₂(uε⁻¹) ≥ T₂(u). Then T₂(u) ≤ C := max(h(A), h(B)), where
h(q) := q⁻² + 2q, A = 410137/10⁶ and B = 417211/125000. Here C = h(B) = 6.765141… ≤ 6.76515 [06: 4.3].

*Proof.*
- Here xq² = 1, so x/q = q⁻³ ∈ [ρ_lo, ρ_hi]. Hence q³ ∈ (1/14.49483, 1/0.02689446).
- So A < q < B, since A³·14.49483 ≤ 1 ≤ B³·0.02689446 [06: 4.3].
- T₂(u) = h(q). By the identities h(q) − h(A) = (q − A)(2 − (q + A)/(q²A²)) and
  h(q) − h(B) = (q − B)(2 − (q + B)/(q²B²)), with (q + A)/A² ≥ (q + B)/B², we get h(q) ≤ max(h(A), h(B)) for
  A ≤ q ≤ B [06: 4.3]: if the first bracket is ≤ 0, then h(q) ≤ h(A); otherwise the second bracket is positive too,
  and h(q) ≤ h(B). (Equivalently, h is convex.) ∎

**(4.3′) Scaled form.** Let y be a nonzero rational triple with |N(y)| = n. If T₂(yε) ≥ T₂(y) and T₂(yε⁻¹) ≥ T₂(y),
then T₂(y) ≤ C·n^{2/3}.

*Proof.* Put w = q/n^{2/3}. Then x/q = w⁻³ and T₂(y) = n^{2/3}h(w); apply the argument of (4.3) to w [06: 4.3]. ∎

*The exact threshold.* Remark 2.5 uses (4.3′) in the cubed form T₂(y)³ ≤ K₃·N(y)², with the threshold
K₃ = 4(E² + E + 1)³/(E²(E + 1)²) = 309.6209… in place of C³. Indeed T₂(y)³/N(y)² = ψ(x/q) with ψ(t) = (t + 2)³/t, and
ψ decreases on (0, 1] and increases on [1, ∞); so on [ρ_lo, ρ_hi] its largest value is ψ(ρ_lo) = ψ(ρ_hi) = K₃
[05: 2.5.2]. K₃ is unchanged under E ↦ 1/E.

**(4.4) The units of small height.** 06 lists the integer triples with T_low ≤ C, where T_low is built as in (2.5) with
ρ = 2⁻³⁰ in place of 2⁻²⁰. T_low ≤ T₂ by the inequality of (2.5), so the list contains every triple with T₂ ≤ C. It
has five triples, 0, ±1 and ±λ (T₂(±λ) = 3λ₁² − 2 ≈ 6.384), and the only units among them are ±1 [06: 4.4].

**(4.5) Proof of (4.1).** Let u be a unit.
- The set {uε^k : T₂(uε^k) ≤ T₂(u)} is finite, since T₂ is positive definite. It is nonempty.
- Let v = uε^j be an element of it of least T₂. Then T₂(vε) ≥ T₂(v) and T₂(vε⁻¹) ≥ T₂(v), so T₂(v) ≤ C by (4.3).
- So v = ±1 by (4.4), and u = ±ε^{−j}.
- *Uniqueness.* 0 < σ₁(ε) < 1, so σ₁(ε^k) = σ₁(ε)^k is not ±1 unless k = 0.
- Finally write k = 2k₁ + l with l ∈ {0, 1}: then ±ε^k = ±ε^l(ε^{k₁})². ∎

**(4.6) Remark.** By (4.1), ε is a fundamental unit; this confirms Uchiyama's unit, which he took from the
Delone–Faddeev table. In addition, the four classes ±1 and ±ε are distinct modulo squares of units [06: 4.5]: a square
has σ₁ > 0, while σ₁(−1) < 0 < σ₁(ε); and a unit v with v² = ε would have T₂(v) = e + 2√E ≤ 5.83 < C, so v = ±1 by
(4.4), which is false.

## 5. The class of X − θ

**(5.1) Proposition (Theorem H).** Let (X, Y) be an integer solution of (E). Then X − θ = μ₀S² with S ∈ R and
μ₀ ∈ {1, π₃ₐπ₃꜀, π₃ₐπ₃ᵦ, π₃ᵦπ₃꜀}.

The proof takes (5.2)–(5.7).

**(5.2) Two factors.** Put α = X − θ = (X, −3, 0) and β = X² + θX + θ² − 9 = (X² − 9, 3X, 9). Then [07: 5.1]

    αβ = X³ − 9X + 81 = Y²,    N(α) = X³ − 9X + 81 = Y²,    β − (X + 2θ)α = 9π₂₃₉ = (π₃ₐπ₃ᵦπ₃꜀)²π₂₃₉.

Since α ≠ 0, N(α) ≠ 0 by (1.4), so Y ≠ 0 and β ≠ 0.

**(5.3) Common divisors.** Let 𝒮 = {π₃ₐ, π₃ᵦ, π₃꜀, π₂₃₉}. These are prime elements by (3.3), and they are pairwise
non-associate: no one of them divides another ((1.5); [07: 5.2]). A common divisor δ of α and β divides
β − (X + 2θ)α = (π₃ₐπ₃ᵦπ₃꜀)²π₂₃₉.

**(5.4) Separating 𝒮.**
- By (3.4), applied to each ϖ ∈ 𝒮 in turn, write
  - α = Π_{ϖ∈𝒮} ϖ^{e_ϖ}·α′,
  - β = Π_{ϖ∈𝒮} ϖ^{f_ϖ}·β′,
  - Y = Π_{ϖ∈𝒮} ϖ^{w_ϖ}·Y′,

  where no element of 𝒮 divides α′, β′ or Y′. Each later quotient divides the earlier ones, so a ϖ that was divided
  out to its exact power cannot divide the final quotient.
- *Matching exponents.* Fix ϖ ∈ 𝒮 and write αβ = Y² as ϖ^{e+f}A = ϖ^{2w}B, where A = α′β′·Π_{ϖ′≠ϖ}ϖ′^{e′+f′} and
  B = Y′²·Π_{ϖ′≠ϖ}ϖ′^{2w′}. Neither A nor B is divisible by ϖ, by (3.3): ϖ divides none of α′, β′, Y′, and none of the
  other elements of 𝒮 (5.3).
  - If e + f > 2w, cancelling ϖ^{2w} (1.4) gives ϖ | B.
  - If e + f < 2w, it gives ϖ | A.

  Both are impossible, so e_ϖ + f_ϖ = 2w_ϖ for every ϖ, and α′β′ = Y′².
- *Coprime.* Let δ be a common divisor of α′ and β′. Then δ divides α and β, hence (π₃ₐπ₃ᵦπ₃꜀)²π₂₃₉ (5.3).
  - For each ϖ ∈ 𝒮, a gcd of δ and ϖ is a unit or an associate of ϖ (3.3).
  - An associate would give ϖ | δ | α′, which is impossible.
  - So δ is coprime to every ϖ ∈ 𝒮, hence to (π₃ₐπ₃ᵦπ₃꜀)²π₂₃₉ by (3.2′). Since δ divides that product, δ | 1.

  So α′ and β′ are coprime.
- *Allocation.* (3.5) with u = 1 and W = Y′ gives α′ = u′R₀², with u′ a unit. Thus

      α = u′ · π₃ₐ^{e_a} π₃ᵦ^{e_b} π₃꜀^{e_c} π₂₃₉^{e₂₃₉} · R₀².

**(5.5) The norm.** Take norms in the last display. With N(π₃ₐ) = N(π₃ᵦ) = N(π₃꜀) = −3, N(π₂₃₉) = 239 and N(u′) = ±1:

    Y² = N(u′)·(−3)^{e_a+e_b+e_c}·239^{e₂₃₉}·N(R₀)².

Y ≠ 0 and N(R₀) ≠ 0 (since α′ ≠ 0) are integers.
- Compare the exponents of 239 and 3 on both sides: e₂₃₉ and e_a + e_b + e_c are even.
- Compare signs: N(u′) = 1.

**(5.6) Eight representatives.** By (4.1), u′ = σε^k with σ = ±1, and N(u′) = σ; so σ = 1 by (5.5). Write k = 2k₁ + l
and e_x = 2m_x + ν_x (x = a, b, c), with l, ν_x ∈ {0, 1}. Then

    X − θ = μ₀S²,    μ₀ = ε^l π₃ₐ^{ν_a} π₃ᵦ^{ν_b} π₃꜀^{ν_c},    S = ε^{k₁}π₃ₐ^{m_a}π₃ᵦ^{m_b}π₃꜀^{m_c}π₂₃₉^{e₂₃₉/2}R₀.

S ∈ R: ε^{k₁} ∈ R for every integer k₁ because ε⁻¹ = (2, −2, 1) ∈ R (1.5), and every other exponent is nonnegative.
ν_a + ν_b + ν_c is even, since e_a + e_b + e_c is. With l ∈ {0, 1} this leaves eight μ₀ [07: 5.3]:
- the four of (5.1), T := {1, π₃ₐπ₃꜀, π₃ₐπ₃ᵦ, π₃ᵦπ₃꜀};
- their products with ε: ε, επ₃ₐπ₃꜀, επ₃ₐπ₃ᵦ, επ₃ᵦπ₃꜀. This set is called εT below.

**(5.7) The step modulo 4.** S has integer coordinates, and (S + 2Z)² = S² + 4(SZ + Z²) for every Z ∈ R. So μ₀S² modulo
4R, that is coordinatewise modulo 4, depends only on S modulo 2, and takes at most eight values [07: 5.4].
X − θ = (X, −3, 0) needs its λ-coordinate ≡ 1 and its λ²-coordinate ≡ 0 (mod 4).

For the four μ₀ of εT this never happens [07: 5.4]. The eight values of μ₀S² modulo 4 are
- (0,0,0), (1,2,0), (2,2,1), (2,3,0), (2,3,3), (3,1,1), (3,2,1), (3,3,2) for μ₀ = ε and επ₃ₐπ₃꜀;
- (0,0,0), (0,1,2), (0,2,3), (1,0,1), (1,0,2), (1,3,2), (2,1,3), (3,1,3) for μ₀ = επ₃ᵦπ₃꜀ and επ₃ₐπ₃ᵦ.

None has the shape (∗, 1, 0).

For the four μ₀ of T, residues S with the required shape exist, as they must, since integral points exist in each
class [07: 5.4]. The step then also fixes X modulo 4: X ≡ 3 for μ₀ ∈ {1, π₃ₐπ₃꜀}, and X ≡ 0 or 1 for
μ₀ ∈ {π₃ₐπ₃ᵦ, π₃ᵦπ₃꜀}. All 12 integral X obey this [07: 5.6].

This proves (5.1). ∎

**(5.8) Hypothesis (H).** By [07: 5.5], 3 − θ = π₃ₐπ₃꜀·π₃ᵦ², −3 − θ = π₃ₐπ₃ᵦ·π₃꜀² and −θ = π₃ᵦπ₃꜀·π₃ₐ². So
μ₀ = μ·τ⁻² for the matching μ ∈ {1, 3 − θ, −3 − θ, −θ} and τ ∈ {1, π₃ᵦ, π₃꜀, π₃ₐ}, and X − θ = μ(S/τ)². This is
hypothesis (H) of proof.md §1.2.

The four μ are pairwise distinct modulo K*².
- *Characters.* For a prime p and a root r of f modulo p, a + bλ + cλ² ↦ a + br + cr² (mod p) is a ring map R → F_p
  (1.5). So χ_{p,r}(α) = ((a + br + cr²)/p), a Legendre symbol, is multiplicative on R.
- The characters at the roots 4 mod 7 and 8 mod 13 take the four different values (1, 1), (−1, −1), (−1, 1) and
  (1, −1) on 1, 3 − θ, −3 − θ and −θ [07: 5.5]. All these values are nonzero.
- A character with nonzero values on μ and μ′ takes the same value on both when μ = μ′ξ² with ξ ∈ K. Indeed ζ = μ′ξ
  has ζ² = μμ′ ∈ R, so ζ ∈ R by (3.6), and χ(μ)χ(μ′) = χ(ζ)². The left side is nonzero, so it is 1.

As a confirmation not used by the proof, each of the 24 integral points has X − θ = μ₀S² with S ∈ R for exactly one μ₀
of (5.1), in the class that the 2-descent assigns to it [07: 5.6].

**(5.9) Remark: what the step modulo 4 is.** It is a consequence, read modulo 4, of the local condition at 2 of the
2-descent, here for integral points with S integral. 8604.a1 has bad reduction at 2; 2 is inert in K, and R/2R is the
field with eight elements, where every element is a square. The step does not separate 1 from π₃ₐπ₃꜀, nor π₃ₐπ₃ᵦ from
π₃ᵦπ₃꜀: their values modulo 4 coincide [07: 5.4]. The norm condition of (5.5) is the condition "norm a square". The
local conditions at 3 and 239 are not needed.

## 6. What changes in proof.md

- *§1.2.* Hypothesis (H) is now Theorem H. The rank of E(Q), the Mordell–Weil theorem (N5) and the LMFDB value are no
  longer used.
- *§1.3.* The descent works in O_K = Z[λ], which is R. It uses these facts about Ξ_μ = {ξ ∈ K : μξ² ∈ R}, for the four
  μ (each in R):
  - *R ⊆ Ξ_μ = −Ξ_μ,* since μ ∈ R.
  - *Closed under addition.* μ(ξ + η)² = μξ² + 2μξη + μη², and μξη ∈ R: indeed (μξη)² = (μξ²)(μη²) ∈ R, so μξη ∈ R by
    (3.6).
  - *Bounded denominators.* If N(μ) | k², then kξ ∈ R for every ξ ∈ Ξ_μ. Indeed (kξ)² = (k²/N(μ))·μ*·(μξ²) ∈ R, so
    kξ ∈ R by (3.6).

  So R ⊆ Ξ_μ ⊆ (1/9)R (N(μ) = 1 or 81), Ξ_μ is a full lattice, and the descent script's enumeration of (1/9)R/R finds it.
  "O_K" in the script may be read as R; (3.6) replaces N2 and makes N1 unnecessary.
- *Theorem B.* It now holds without the rank, by Theorem A, proof.md §1 (with the parametrization lemma 1.3′ added
  there), proof.md §5 (the reconstruction of the points from the orbits) and Theorem H.

The 2-descent itself is unchanged: Theorem H only says in which class X − θ lies.

## 7. Provenance and novelty

Details and reading marks are in [literature_check.md](literature_check.md) §6.
- **The field.** K = Q(λ), λ³ − λ + 3 = 0, has discriminant −239 [04: 1.2]. Uchiyama (1984, p. 133) states class
  number 1 without a source and takes the unit λ² + λ − 1 from the Delone–Faddeev table. LMFDB 3.1.239.1 lists class
  number 1 and the fundamental unit ε.
- **Norm-Euclidean.** That Z[λ] is norm-Euclidean is E. M. Taylor's theorem (thesis 1975; J. London Math. Soc. (2) 14,
  1976), on Lemmermeyer's attribution; Taylor's text was not read. Cavallar and Lemmermeyer list −239 as "E 8/9",
  norm-Euclidean with Euclidean minimum 8/9 (*The Euclidean algorithm in cubic number fields*, 1998; arXiv:1202.6019v1,
  §7 table, p. 16, read in review round 3), and Lemmermeyer's 1989 table gives the same value. Their value comes from a
  floating-point covering with exact checks at the exceptional points; they trust it through agreement with earlier
  work and do not present it as a proof. Remark 2.5 proves the lower bound 8/9 exactly, at the point ξ₀ fixed by ε. The
  table of §2 proves M ≤ 32733/32768 < 1; M ≤ 8/9 is not proved here, and it is not needed.
- **Euclidean division by a box covering** is the standard method for Euclidean minima (Cavallar–Lemmermeyer; Lezowski,
  Math. Comp. 83, 2014). Here each box is checked in exact rational arithmetic.
- **The units.** The descent on T₂ is the one-unit version of the trace-height descent of mordell_1963 §4, as in
  [norm.triples.multiplier](../../methods/norm.triples.multiplier.md) for 1611.a1.
- **(H) itself.** Uchiyama's §1 proves the analogue of (H) for the points with 3 | x, with ideals, his class number
  (stated without a source) and his tabulated unit ([literature_check.md](literature_check.md) §1). No ideal-free or
  table-free proof of (H) for this curve was found in the sources checked (literature_check.md §6). Here the same
  factorization is made table-free.
- **Novelty.** None is claimed, by the project's criterion. The mechanism is classical: a Euclidean algorithm in the
  cubic field (Taylor's theorem, here certified in exact arithmetic), the divisibility theory it gives, units by
  descent, and the condition at 2 of the 2-descent. The table and the descent are presented as certificates.
