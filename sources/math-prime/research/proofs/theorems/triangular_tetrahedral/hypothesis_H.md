# Hypothesis (H) without the rank: arithmetic of Z[λ], λ³ = 4λ − 6

Written 2026-10-15 by Claude (Opus 5.5). It proves hypothesis (H) of [proof.md](proof.md) §1.2 for integral points,
with elementary means only. With it, Theorem B of proof.md no longer depends on the rank of E(Q). The method is that
of [norm.triples](../../methods/norm.triples.md) (Mordell's equation, [mordell_1963](../mordell_1963/proof.md)
§§1–5). It needs one new ingredient, because Z[λ] is **not** norm-Euclidean (Remark 2.5): a division with a
multiplier s ∈ {1, 2} (Dedekind–Hasse). Scripts: [scripts/](scripts/), each with its output beside it as `*.out.txt`.
Every check in them is exact. Floating point only seeds some searches, and each result is then verified in exact
arithmetic.

**Theorem H.** Let (X, Y) be an integer solution of

    Y² = X³ − 144X + 1296                                              (E)

and let θ = 6λ, where λ³ = 4λ − 6, so that θ³ − 144θ + 1296 = 0. Then

    X − θ = μ₀ S²   with S ∈ Z[λ] and μ₀ ∈ {1, επ₃ₐπ₃ᵦ, −π₃ₐπ₃꜀, −επ₃ᵦπ₃꜀},

where ε = λ² − 5λ − 19, π₃ₐ = λ − 1, π₃ᵦ = 2λ + 5 and π₃꜀ = λ² − 3λ + 3. Moreover

    9 − θ = επ₃ₐπ₃ᵦ·(−45 + 48λ − 19λ²)²,   −15 − θ = −π₃ₐπ₃꜀·(5 + 2λ)²,   −12 − θ = −επ₃ᵦπ₃꜀·(−8 + 8λ − 3λ²)².

So X − θ ∈ μK*² with μ ∈ {1, 9 − θ, −15 − θ, −12 − θ}, which is hypothesis (H).

**What the proof uses.**
- Integer and rational arithmetic, polynomial identities, congruences, induction (descent), and real inequalities.
- Real numbers enter only through rational interval enclosures of λ₁, the one real root of t³ − 4t + 6, obtained from
  exact sign changes. They are used for the units (Section 4) and for Remark 2.5, which Theorem H does not use. The
  division of Section 2 needs only rational triples and a table of 398 boxes, checked in exact rational arithmetic.
- No ideal theory, class number, reciprocity law, unit theorem, Mordell–Weil theorem or rank is used.
- Avanesov took the class number 1 of Q(λ) from the Delone–Faddeev table, and a fundamental unit computed by
  Billevich's method. Here both are replaced by proofs. Sections 2–3 give the divisibility theory that class number 1
  provided (Bézout, prime elements, square allocation), and Section 4 gives the unit.
- Section 6 proves unique factorization a second time, by the box principle, as an independent check. It is not part
  of the proof of Theorem H. It works with the kernels P(p, g), which are prime ideals in modern language (Remark after
  (6.6)), so only Sections 1–5 meet the no-ideal-theory standard.

**Where (H) was used, and what replaces it.** Every item below is proved here.

| input of proof.md, Theorem B | before | now |
|---|---|---|
| (H), the list of four classes (§1.2) | rank E(Q) ≤ 2 (LMFDB, mwrank) with Mordell–Weil (N5) | Theorem H, Sections 1–5 |
| N2: integral elements form a ring; α² integral implies α integral (§1.3) | standard facts, not proved | (3.6): Z[λ] is integrally closed |
| N1: O_K = Z-span(1, θ/6, θ²/36), via the index–discriminant formula (§1.3) | p-maximality tests | (3.6): Z[λ] is the set of algebraic integers of K |

**Checks.** A reference such as [04: 1.5] means section 1.5 of `04_ring.py`, with the result in `04_ring.out.txt`.
Shared definitions are in [cubic.py](scripts/cubic.py).

| Script | Checks |
|---|---|
| [04_ring.py](scripts/04_ring.py) ([output](scripts/04_ring.out.txt)) | ring laws, norm, adjugate, nonvanishing of the norm, special elements, real values (Section 1) |
| [05_dh_cover.py](scripts/05_dh_cover.py) ([output](scripts/05_dh_cover.out.txt)) | the box table of the division with multiplier, its cover of D, the half-lattice points, the point ξ* (Section 2) |
| [06_units.py](scripts/06_units.py) ([output](scripts/06_units.out.txt)) | root interval, descent constants, the units of small height (Section 4) |
| [07_classes.py](scripts/07_classes.py) ([output](scripts/07_classes.out.txt)) | the two factors, their common divisors, the norm condition, the steps modulo 4 and 8, the identification (Section 5) |
| [08_ufd_box.py](scripts/08_ufd_box.py) ([output](scripts/08_ufd_box.out.txt)) | a second proof of Section 3 by the box principle (Section 6) |
| [cover_search_K.py](scripts/cover_search_K.py) ([output](scripts/cover_search_K.out.txt)) | how the box table [dh_boxes.txt](scripts/dh_boxes.txt) was found (not used by the proof) |
| verify_H_units.py, verify_H_dh.py, verify_H_classes.py, verify_H_box.py (with outputs) | independent re-implementations of Sections 1 and 4, 2, 5 and 6 (review.md, round 4) |

---

## 1. The ring of integer triples

**(1.1) Multiplication.** A triple (a, b, c) stands for a + bλ + cλ², where λ³ = 4λ − 6, so λ⁴ = 4λ² − 6λ. Define

    (a,b,c)(d,e,f) = (ad − 6bf − 6ce,  ae + bd + 4bf + 4ce − 6cf,  af + be + cd + 4cf).

This is the product of a + bt + ct² and d + et + ft², reduced modulo f(t) = t³ − 4t + 6 [04: 1.1]. With
coordinatewise addition the triples form a commutative ring with unit (1, 0, 0) [04: 1.1]. Write R for the ring of
integer triples and R_Q for the ring of rational triples. An integer r is the triple (r, 0, 0). Put θ = 6λ = (0, 6, 0);
then θ³ − 144θ + 1296 = 0, and X − θ is the triple (X, −6, 0).

**(1.2) Norm.** Let M(α) be the matrix whose columns are α·1, α·λ, α·λ². For α = (a, b, c),

    M(α) = [[a, −6c, −6b], [b, a + 4c, 4b − 6c], [c, b, a + 4c]],

and M(αβ) = M(α)M(β) [04: 1.2]. The norm is the cubic form

    N(α) = det M(α) = a³ + 8a²c − 4ab² + 18abc + 16ac² − 6b³ + 24bc² + 36c³.

The proof uses only the polynomial identities N(αβ) = N(α)N(β), N(rα) = r³N(α) and N(−α) = −N(α) [04: 1.2]. Also
N(X − θ) = X³ − 144X + 1296 identically in X [04: 1.5].

**(1.3) Adjugate triple.** For α = (a, b, c) put

    α* = (a² + 8ac − 4b² + 6bc + 16c²,  −ab − 6c²,  b² − ac − 4c²).

Then α·α* = (N(α), 0, 0) [04: 1.3]. Consequences, as in mordell_1963 (1.3):
- **(1.3a)** α ∈ R is a unit of R exactly when N(α) = ±1.
- **(1.3b)** A rational triple α with N(α) ≠ 0 has the inverse α*/N(α) in R_Q.

**(1.4) Lemma: the norm vanishes only at 0.** If α ∈ R_Q and N(α) = 0, then α = 0.

*Proof.* Suppose α ≠ 0. Scale it to integers (a, b, c) with gcd 1 and N(a, b, c) = 0. Then f is Eisenstein at 2, and
the descent of mordell_1963 (1.4) goes through [04: 1.4]:
- N − a³ has even coefficients, so a is even: a = 2a′.
- N(2a′, b, c) = 2P, where P + 3b³ has even coefficients. So b³ is even, and b = 2b′.
- N(2a′, 2b′, c) = 4P′, where P′ − 9c³ has even coefficients. So c is even, a contradiction. ∎

So R has cancellation. We write β | α when α = βγ with γ ∈ R; then N(β) | N(α). Associates differ by a unit.

**(1.5) Special elements** [04: 1.5]:

| element | triple | norm | identity |
|---|---|---|---|
| ε = λ² − 5λ − 19 | (−19, −5, 1) | 1 | ε⁻¹ = 40λ² − 101λ + 95 = (95, −101, 40) |
| π₂ = λ² + λ − 4 | (−4, 1, 1) | −2 | π₂³ = −2ε, so 2 = −π₂³ε⁻¹ |
| π₃ₐ = λ − 1 | (−1, 1, 0) | −3 | π₃ₐπ₃ᵦπ₃꜀ = 3 |
| π₃ᵦ = 2λ + 5 | (5, 2, 0) | −3 | |
| π₃꜀ = λ² − 3λ + 3 | (3, −3, 1) | 3 | |
| π₁₇₉ = 45λ² − 114λ + 107 | (107, −114, 45) | 179 | π₂²π₁₇₉ = 3λ² − 4 |

Avanesov's unit is ε₀ = ε⁻¹, and his factorizations of 2, 3 and 3λ² − 4 are these.

**(1.6) Real values.** f(−3) = −9 < 0 < 6 = f(−2), so f has a root λ₁ ∈ (−3, −2). Bisection encloses λ₁ in a
rational interval of width below 2⁻²⁰⁰ [04: 1.6]; in particular λ₁ < −5/2.
- f(t) = (t − λ₁)(t² + λ₁t + λ₁² − 4), and t² + λ₁t + λ₁² − 4 = (t + λ₁/2)² + (3λ₁² − 16)/4.
- 3λ₁² − 16 > 3(5/2)² − 16 > 0, so the quadratic factor has no real root, and λ₁ is the only real root of f.
- The quadratic factor has the complex pair z, z̄ as roots, with z + z̄ = s := −λ₁ and zz̄ = p := λ₁² − 4.

For α = (a, b, c) put σ₁(α) = a + bλ₁ + cλ₁² and

    Q(α) = (a + bz + cz²)(a + bz̄ + cz̄²) = a² + sab + (s² − 2p)ac + pb² + spbc + p²c².

The right-hand side has real coefficients, so no complex number is needed. With s = −λ₁ and p = λ₁² − 4 it reduces,
modulo f(λ₁) = 0, to

    Q(α) = a² − λ₁ab + (8 − λ₁²)ac + (λ₁² − 4)b² + 6bc + (16 − 6λ₁ − 4λ₁²)c² = σ₁(α*).

The following hold as polynomial identities modulo f(λ₁) = 0 [04: 1.6]:
- σ₁(α)Q(α) = N(α);
- σ₁(αβ) = σ₁(α)σ₁(β);
- Q(αβ) = Q(α)Q(β).

*Signs.* Q is a sum of two squares with positive coefficients, Q = (a − (λ₁/2)b + ((8 − λ₁²)/2)c)² +
((3λ₁² − 16)/4)(b − λ₁c)², and 3λ₁² − 16 > 0 [04: 1.6]. So Q ≥ 0, but Q is only semidefinite: it vanishes on the
real line through (λ₁² − 4, λ₁, 1). A nonzero rational triple is not on that line, because λ₁ is irrational; and
anyway N = σ₁Q ≠ 0 on it by (1.4). So Q(α) > 0 for every nonzero rational triple α.

*Height.* Put T₂(α) = σ₁(α)² + 2Q(α). Modulo f(λ₁),

    T₂ = 3a² + 16ac + (3λ₁² − 8)b² + 8λ₁bc + (32 − 18λ₁ − 4λ₁²)c²,

with leading minors 3, 9λ₁² − 24 > 0 and 716 [04: 1.6]. So T₂ is a positive definite quadratic form, and only finitely
many integer triples have T₂ below any given bound. Completing squares,

    T₂ = 3(a + 8c/3)² + d₂(b + kc)² + d₃c²,    d₂ = 3λ₁² − 8,  k = 4λ₁/d₂,  d₃ = 716/(3d₂).

## 2. Division with a multiplier

**(2.1) Lemma.**
- (a) For every rational triple ξ in D = [−½, ½] × [−½, ½] × [0, ½] there are s ∈ {1, 2} and κ ∈ R with
  |N(sξ − κ)| < 1.
- (b) For every rational triple ξ with 2ξ ∈ R and ξ ∉ R there is κ ∈ R with 0 < |N(ξ − κ)| < 1.

*Proof of (a).* The table [dh_boxes.txt](scripts/dh_boxes.txt) lists 398 closed boxes in D. Each has ends that are
multiples of 1/64, a multiplier s ∈ {1, 2} and a translation κ ∈ R.

*The bound on a box.* For a box with centre ξ₀ and half-widths h, put v = sξ₀ − κ and w = s(ξ − ξ₀), so |wᵢ| ≤ shᵢ.
The exact Taylor expansion of the cubic N [05: 2.1],

    N(v + w) = N(v) + Σᵢ gᵢ(v)wᵢ + Σ_{i≤j} q_{ij}(v)wᵢwⱼ + N(w),

with, for v = (x, y, z),

    g₁ = 3x² + 16xz − 4y² + 18yz + 16z²,   g₂ = −8xy + 18xz − 18y² + 24z²,   g₃ = 8x² + 18xy + 32xz + 48yz + 108z²,
    q₁₁ = 3x + 8z,  q₁₂ = −8y + 18z,  q₁₃ = 16x + 18y + 32z,  q₂₂ = −4x − 18y,  q₂₃ = 18x + 48z,  q₃₃ = 16x + 24y + 108z,

gives |N(sξ − κ)| ≤ B on the box, where

    B = |N(v)| + Σᵢ |gᵢ(v)|shᵢ + Σ_{i≤j} |q_{ij}(v)|s²hᵢhⱼ + Σ_m |n_m|(sh)^m

and n_m runs over the coefficients of N (Σ|n_m| = 113). B is computed exactly, and B < 1 for every box; the largest
is 8191/8192 [05: 2.2].

*The boxes cover D.* Each of the 64·64·32 = 131072 cells of side 1/64 that tile D lies in exactly one box. So every
point of D lies in a closed cell, hence in a box [05: 2.3]. Two more checks confirm the cover:
- the boxes are exactly the leaves of a binary tree of midpoint splits of D, listed in depth-first order;
- their volumes add up to vol(D) = ½.

So the listed (s, κ) gives |N(sξ − κ)| ≤ B < 1 at every point of D. The argument works for real ξ as well; only
rational ξ are needed. ∎

*Proof of (b).* Modulo R there are seven such ξ, the triples x/2 with x ∈ {0, 1}³ \ {0}. With the listed κ, the
values of N(ξ − κ) are [05: 2.4]:

| ξ | (0,0,½) | (0,½,0) | (0,½,½) | (½,0,0) | (½,0,½) | (½,½,0) | (½,½,½) |
|---|---|---|---|---|---|---|---|
| κ | (−1,1,0) | (0,0,0) | (2,0,0) | (0,0,0) | (0,1,1) | (0,1,0) | (−1,2,0) |
| N(ξ − κ) | ½ | −¾ | −¼ | ⅛ | −⅞ | ⅜ | ⅜ |

All lie strictly between 0 and 1 in absolute value. ∎

**(2.2) Corollary.** For every rational triple ξ ∉ R there are s ∈ {1, 2} and κ ∈ R with 0 < |N(sξ − κ)| < 1.

*Proof.*
- If 2ξ ∈ R, apply (b) with s = 1.
- Otherwise, let r be the triple of nearest integers to the coordinates of ξ, and put ξ′ = ξ − r, or ξ′ = r − ξ if the
  third coordinate of ξ − r is negative. Then ξ′ ∈ D.
- Take s, κ′ for ξ′ from (a). If ξ′ = ξ − r, then sξ − (κ′ + sr) = sξ′ − κ′. If ξ′ = r − ξ, then
  sξ − (sr − κ′) = −(sξ′ − κ′). Either way we get κ with |N(sξ − κ)| = |N(sξ′ − κ′)| < 1.
- Since ξ ∉ R and 2ξ ∉ R, sξ ∉ R. So sξ − κ ≠ 0, and N(sξ − κ) ≠ 0 by (1.4). ∎

**(2.3) Division with a multiplier.** Let α, β ∈ R with β ≠ 0 and β ∤ α. Then there are s ∈ {1, 2} and κ ∈ R with

    0 < |N(sα − κβ)| < |N(β)|.

*Proof.* N(β) ≠ 0 by (1.4). The rational triple ξ = α·β*/N(β) satisfies ξβ = α, and ξ ∉ R because β ∤ α. Take s, κ
from (2.2) and put ρ = sα − κβ = (sξ − κ)β ∈ R. Then |N(ρ)| = |N(sξ − κ)|·|N(β)|, and 0 < |N(sξ − κ)| < 1. ∎

This is the Dedekind–Hasse criterion with the multipliers {1, 2}.

**(2.4) Remark: the table.**
- *Shape.* 263 boxes use s = 1 and 135 use s = 2. 375 boxes have all |κᵢ| ≤ 4; the largest |κᵢ| is 22.
  The deepest box is reached by 15 midpoint splits [05: 2.2].
- *How it was found.* A search ([cover_search_K.py](scripts/cover_search_K.py)) tried translates near ε^k·sξ₀ for
  k ∈ {−1, 0, 1}, and five split rules. It kept the smallest of the five tables. The proof uses only the table.
- *Sharp bound.* The Taylor bounds are not sharp. Expand N(sξ − κ) on each box in the 64 products of cubic Bernstein
  polynomials, as in mordell_1963 (2.3). The largest coefficient over all boxes is 4093/4096. It is a corner value,
  at ξ = (−5/16, −5/16, ¼), in the box with s = 1 and κ = (−1, 0, 0). So the largest value of |N(sξ − κ)| over the
  boxes, each on its own closed box, is exactly 4093/4096 [05: 2.6]. A referee checked that the point lies in five
  other boxes, where the values are 2829/4096 (three times), −3267/4096 and −3603/4096 (review.md, round 4).
- *Confirmation.* 3302 random and planted pairs (α, β), 3301 of them with β ∤ α, were divided by (2.3), with the
  expected inequality every time [05: 2.7].

**(2.5) Remark: R is not norm-Euclidean.** The multiplier cannot be dropped. Let

    ξ* = 79/109 + (247/327)λ + (142/327)λ².

Then εξ* − ξ* = (−6, −27, −10) ∈ R, so the coset C = ξ* + R satisfies εC = C. We have N(ξ* − 2) = −121/109.
Moreover |N(ξ* − κ)| ≥ 121/109 > 1 for every κ ∈ R [05: 2.5]. So 327ξ* = 237 + 247λ + 142λ² cannot be divided by
327 with a remainder of smaller norm.

*Proof of the bound.* Suppose some y ∈ C has |N(y)| < 121/109.
- The orbit {ε^k y : k ∈ Z} lies in C, since εC = C and ε⁻¹C = C. Its elements with T₂ ≤ T₂(y) are finitely many,
  since T₂ is positive definite on the lattice C. So the orbit has an element y₀ of least T₂. Also N(y₀) = N(y).
- Neither εy₀ nor ε⁻¹y₀ has smaller T₂, so T₂(y₀)³ ≤ K₃·N(y)², where K₃ = 4(E² + E + 1)³/(E²(E + 1)²) with E = Q(ε)
  as in (4.2), and C₀ = K₃^{1/3} ≤ 113.6238 [05: 2.5.2].
  - This is the scaled lemma (4.3′) with its exact threshold. The rounded 114 of (4.3′) would need a slightly larger
    list below.
  - The expression is unchanged under E ↦ 1/E. 05 writes it with σ₁(ε) = 1/E, and bounds K₃ above by its value at
    the end of the rational enclosure of σ₁(ε).
- So T₂(y₀)³ < K₃·(121/109)² ≤ T_e³, with T_e := 60909/500 = 121.818.
- [05: 2.5.3] uses T_low(y) = yᵀ(G′ − 3ρI)y. Here G′ is a rational matrix within ρ = 2⁻²⁰ of the Gram matrix G of T₂,
  entrywise; that follows from the rational enclosure of λ₁.
  - T_low ≤ T₂: T₂ − T_low = yᵀ(G − G′)y + 3ρΣyᵢ² ≥ −ρ(Σ|yᵢ|)² + 3ρΣyᵢ² ≥ 0.
  - T_low is positive definite, by an exact LDLᵀ.
  - 05 lists the 204 elements y ∈ C with T_low(y) ≤ T_e. The list contains every y ∈ C with T₂(y) ≤ T_e. The least |N|
    on it is 121/109, attained only at ξ* − 2.
- y₀ is on that list and has |N(y₀)| < 121/109, a contradiction. ∎

*The other fixed points.* There are |N(ε − 1)| = 654 classes of points fixed by ε modulo R. The same argument, with 1 in
place of 121/109, shows that only the classes of ±ξ* have no translate of norm below 1. The next largest minimum among
those found is 99/109 [05: 2.5.4]. The multiplier s = 2 is what the table uses near ±ξ*: 2ξ* is another fixed point,
with a translate of norm below 1.

So Z[λ] is not norm-Euclidean, and its Euclidean minimum is at least 121/109. A plain Euclidean division is
impossible, so Avanesov's class number 1 cannot be replaced as norm.triples did for Mordell's ring. Section 2 uses a
multiplier instead; Section 6 avoids division altogether.

## 3. Common divisors, primes and square allocation

**(3.1) Gcd and Bézout.** Let α, β ∈ R, not both 0. Then there is δ ∈ R with δ | α, δ | β and δ = xα + yβ for some
x, y ∈ R. Every common divisor of α and β divides δ.

*Proof.* The set L = {xα + yβ : x, y ∈ R} contains a nonzero element. Choose a nonzero δ ∈ L with |N(δ)| least;
|N| is a positive integer on nonzero elements of R, by (1.4).
- Suppose δ ∤ α. Then (2.3) gives s, κ with 0 < |N(sα − κδ)| < |N(δ)|.
- But sα − κδ ∈ L, which contradicts the choice of δ. So δ | α, and likewise δ | β.
- A common divisor of α and β divides xα + yβ = δ. ∎

The rest of this section follows mordell_1963 §3 from (3.1) on, with the same proofs and with 𝓡 replaced by R. Those
proofs use only (1.2) (N is multiplicative and integer-valued on R), (1.3a), (1.4) and (3.1); no Euclidean division
enters after (3.1). The statements and short proofs are repeated here, with these changes:
- (3.2′) is new;
- the list of primes in (3.3) is this ring's;
- (3.6) here is mordell_1963 (3.7), with this ring's characteristic polynomial.

*Gcd and coprime.* We call δ a gcd of α and β. We call α and β coprime if every common divisor is a unit. This holds
exactly when xα + yβ = 1 for some x, y ∈ R.

**(3.2) Euclid's lemma.** If γ | αβ and γ, α are coprime, then γ | β.

*Proof.* From xγ + yα = 1 we get β = xγβ + yαβ, and γ divides both terms. ∎

**(3.2′) Products.** If γ is coprime to α and to β, then γ is coprime to αβ. In particular, powers of two coprime
elements are coprime.

*Proof.* Multiply xγ + yα = 1 and x′γ + y′β = 1: the product is γ(…) + yy′αβ = 1. ∎

Section 5 needs (3.2′), because there several primes are separated at once; mordell_1963 has only one.

**(3.3) Prime elements.** Let ϖ ∈ R with N(ϖ) = ±q for a prime number q. Then:
- ϖ is not a unit;
- every divisor of ϖ is a unit or an associate of ϖ;
- if ϖ | αβ then ϖ | α or ϖ | β.

This applies to π₂ (q = 2), π₃ₐ, π₃ᵦ, π₃꜀ (q = 3) and π₁₇₉ (q = 179).

*Proof.* N(ϖ) ≠ ±1, so ϖ is not a unit (1.3a). If ϖ = δγ, then N(δ)N(γ) = ±q. So δ is a unit, or γ is a unit and δ is
an associate of ϖ. Suppose ϖ | αβ and ϖ ∤ α. A gcd of ϖ and α divides ϖ, and it is not an associate of ϖ. So it is a
unit, ϖ and α are coprime, and ϖ | β by (3.2). ∎

Two prime elements are either associates or coprime: a gcd divides both, so it is a unit or an associate of each.

**(3.4) Exact powers.** Let ϖ be as in (3.3) and α ∈ R nonzero. If ϖⁿ | α then |N(ϖ)|ⁿ divides |N(α)| ≥ 1. So there
is a largest n ≥ 0 with ϖⁿ | α, and then α = ϖⁿα′ with ϖ ∤ α′.

**(3.5) Lemma: coprime square allocation.** Let α, β, W ∈ R with W ≠ 0, and let u be a unit, with

    αβ = uW²    and    α, β coprime.

Then α = u′R₀² for a unit u′ and some R₀ ∈ R.

*Proof.* Induction on |N(W)|, exactly as mordell_1963 (3.5). Let g be a gcd of α and W.
- *g a unit.* Then xα + yW = 1. Squaring, and using W² = u⁻¹αβ, α divides 1. So α is a unit, and α = α·1².
- *g not a unit.* Write α = gα₁ and W = gW₁. Cancelling g in gα₁β = ug²W₁² gives α₁β = ugW₁².
  - g and β are coprime, since a common divisor divides α and β. So g | α₁ by (3.2); write α₁ = gα₂.
  - Then α₂β = uW₁², α₂ and β are coprime, and |N(W₁)| < |N(W)|.
  - By induction α₂ = u′R₁², so α = g²α₂ = u′(gR₁)².

When |N(W)| = 1, W is a unit, so every gcd of α and W is a unit and the first case applies. This starts the
induction. ∎

**(3.6) Remark: integral closedness (replaces N1 and N2 of the descent).** By (1.3b) and (1.4), R_Q is a field; it is
K = Q(λ). R is integrally closed in it.
- *R is integrally closed.* Suppose ξ ∈ R_Q satisfies ξⁿ + c_{n−1}ξ^{n−1} + ⋯ + c₀ = 0 with integers cᵢ. Then
  ξ ∈ R, by the argument of mordell_1963 (3.7).
  - Write ξ = α/m with α ∈ R and m ≥ 1 an integer. Take a gcd δ of α and m (3.1) and write ξ = α′/m′ with α′, m′
    coprime.
  - The equation gives m′ | α′ⁿ. Raising xα′ + ym′ = 1 to the n-th power shows m′ | 1.
- *Every element of R is integral over Z.* Each triple is a root of its characteristic polynomial,
  t³ − (3a + 8c)t² + (3a² + 16ac − 4b² + 18bc + 16c²)t − N(α), which is monic with integer coefficients [04: 1.2].

So the elements of K that are roots of monic integer polynomials are exactly R; that is, O_K = R. In particular, if
ξ ∈ K and ξ² = c ∈ R, then ξ ∈ R: ξ is a root of p_c(t²), where p_c is the monic characteristic polynomial of c.
Section 7 uses this in place of the facts N1 and N2 of the descent script.

## 4. The units

**(4.1) Theorem.** The units of R are exactly ±ε^k (k ∈ Z), and k is determined by the unit. So every unit u can be
written u = σε^l v² with σ ∈ {1, −1}, l ∈ {0, 1} and v a unit.

**(4.2) Values of ε** [06: 4.2].
- σ₁(ε) ∈ (1/606, 1/605).
- E := Q(ε) = σ₁(ε*) = σ₁(ε⁻¹) = 40λ₁² − 101λ₁ + 95 ∈ (605.08, 605.09).
- σ₁(ε)·E = N(ε) = 1. So σ₁(ε⁻¹) = E and Q(ε⁻¹) = 1/E.

**(4.3) Descent lemma.** For a nonzero rational triple y put x = σ₁(y)² and q = Q(y). Both are positive (1.6), with
xq² = N(y)² and T₂(y) = x + 2q. By (1.6) and (4.2) [06: 4.3]:

    T₂(yε) − T₂(y) = (E − 1)(2q − x(E + 1)/E²),      T₂(yε⁻¹) − T₂(y) = (E − 1)(x(E + 1) − 2q/E).

So T₂(yε) < T₂(y) exactly when x/q > ρ_hi := 2E²/(E + 1), and T₂(yε⁻¹) < T₂(y) exactly when x/q < ρ_lo := 2/(E(E + 1)).
ρ_hi increases and ρ_lo decreases with E, so E < 60509/100 gives ρ_hi < 1209 and ρ_lo > 5453·10⁻⁹ [06: 4.3].

*Lemma.* Let u be a unit with T₂(uε) ≥ T₂(u) and T₂(uε⁻¹) ≥ T₂(u). Then T₂(u) < 113.85.

*Proof.*
- Here xq² = 1, so x/q = q⁻³ ∈ [ρ_lo, ρ_hi]. Hence q³ ∈ (1/1209, 10⁹/5453).
- So A < q < B, with A = 469/5000 and B = 2841/50, since A³·1209 ≤ 1 ≤ B³·5453·10⁻⁹.
- T₂(u) = h(q) := q⁻² + 2q. By the identities h(q) − h(A) = (q − A)(2 − (q + A)/(q²A²)) and
  h(q) − h(B) = (q − B)(2 − (q + B)/(q²B²)), with (q + A)/A² ≥ (q + B)/B², we get h(q) ≤ max(h(A), h(B)) = h(A) < 113.85
  for A ≤ q ≤ B [06: 4.3]. (Equivalently, h is convex.) ∎

**(4.3′) Scaled form.** Let y be a nonzero rational triple with |N(y)| = n. If T₂(yε) ≥ T₂(y) and
T₂(yε⁻¹) ≥ T₂(y), then T₂(y) ≤ h(A)·n^{2/3} < 114·n^{2/3}.

*Proof.* Put w = q/n^{2/3}. Then x/q = w⁻³, and T₂(y) = n^{2/3}h(w). Apply the argument of (4.3) to w [06: 4.3]. ∎

**(4.4) The units of small height.** By the completed squares of (1.6), T₂ ≤ 114 forces |c| ≤ 2, then |b| ≤ 5 and
|a| ≤ 11. Of the 1265 triples in this box, exactly 187 have T₂ ≤ 114, each decided by exact interval arithmetic. The
only units among them are ±1 [06: 4.4].

**(4.5) Proof of (4.1).** Let u be a unit.
- The set {uε^k : T₂(uε^k) ≤ T₂(u)} is finite, since T₂ is positive definite. It is nonempty.
- Let v = uε^j be an element of it of least T₂. Then T₂(vε) ≥ T₂(v) and T₂(vε⁻¹) ≥ T₂(v), so T₂(v) < 114 by (4.3).
- So v = ±1 by (4.4), and u = ±ε^{−j}.
- *Uniqueness.* 0 < σ₁(ε) < 1, so σ₁(ε^k) = σ₁(ε)^k is not ±1 unless k = 0. So ε^k = ±1 only for k = 0.
- Finally write k = 2k₁ + l with l ∈ {0, 1}: then ±ε^k = ±ε^l(ε^{k₁})². ∎

**(4.6) Remark.** The four classes ±1 and ±ε are distinct modulo squares of units [06: 4.5]. A square has σ₁ > 0,
while σ₁(−1) < 0 < σ₁(ε). And a unit v with v² = ε would have T₂(v) = 1/E + 2√E < 114, so v = ±1 by (4.4), which is
false. Theorem H does not use this. Avanesov's unit ε₀ = ε⁻¹ is therefore fundamental, which confirms his p. 411
claim, made by Billevich's method.

## 5. The class of X − θ

**(5.1) Proposition (Theorem H).** Let (X, Y) be an integer solution of (E). Then X − θ = μ₀S² with S ∈ R and
μ₀ ∈ {1, επ₃ₐπ₃ᵦ, −π₃ₐπ₃꜀, −επ₃ᵦπ₃꜀}.

The proof takes (5.2)–(5.8).

**(5.2) Two factors.** Put α = X − θ = (X, −6, 0) and β = X² + θX + θ² − 144 = (X² − 144, 6X, 36). Then [07: 5.1]

    αβ = X³ − 144X + 1296 = Y²,    N(α) = X³ − 144X + 1296 = Y²,    β − (X + 2θ)α = 36(3λ² − 4).

Here X + 2θ = (X, 12, 0). Since α ≠ 0, N(α) ≠ 0 by (1.4), so Y ≠ 0 and β ≠ 0.

**(5.3) Common divisors.** By (1.5) [07: 5.2],

    36(3λ² − 4) = ε⁻² · π₂⁸ · (π₃ₐπ₃ᵦπ₃꜀)² · π₁₇₉.

Let 𝒮 = {π₂, π₃ₐ, π₃ᵦ, π₃꜀, π₁₇₉}. These are prime elements by (3.3), and they are pairwise non-associate [07: 5.2].
A common divisor δ of α and β divides β − (X + 2θ)α, hence the right-hand side.

**(5.4) Separating 𝒮.**
- By (3.4), applied to each ϖ ∈ 𝒮 in turn, write
  - α = Π_{ϖ∈𝒮} ϖ^{e_ϖ}·α′,
  - β = Π_{ϖ∈𝒮} ϖ^{f_ϖ}·β′,
  - Y = Π_{ϖ∈𝒮} ϖ^{w_ϖ}·Y′,

  where no element of 𝒮 divides α′, β′ or Y′. Each later quotient divides the earlier ones, so a ϖ that was divided
  out to its exact power cannot divide the final quotient.
- *Matching exponents.* Fix ϖ ∈ 𝒮 and write αβ = Y² as ϖ^{e+f}A = ϖ^{2w}B, where A = α′β′·Π_{ϖ′≠ϖ}ϖ′^{e′+f′} and
  B = Y′²·Π_{ϖ′≠ϖ}ϖ′^{2w′}. Neither A nor B is divisible by ϖ, by (3.3).
  - If e + f > 2w, cancelling ϖ^{2w} (1.4) gives ϖ | B.
  - If e + f < 2w, it gives ϖ | A.

  Both are impossible, so e_ϖ + f_ϖ = 2w_ϖ for every ϖ, and α′β′ = Y′².
- *Coprime.* Let δ be a common divisor of α′ and β′. Then δ divides α and β, hence ε⁻²π₂⁸(π₃ₐπ₃ᵦπ₃꜀)²π₁₇₉ (5.3).
  - For each ϖ ∈ 𝒮, a gcd of δ and ϖ is a unit or an associate of ϖ (3.3).
  - An associate would give ϖ | δ | α′, which is impossible.
  - So δ is coprime to every ϖ ∈ 𝒮, hence to ε⁻²π₂⁸(π₃ₐπ₃ᵦπ₃꜀)²π₁₇₉ by (3.2′). Since δ divides that product, δ | 1.

  So α′ and β′ are coprime.
- *Allocation.* (3.5) with u = 1 and W = Y′ gives α′ = u′R₀², with u′ a unit. Thus

      α = u′ · π₂^{e₂} π₃ₐ^{e_a} π₃ᵦ^{e_b} π₃꜀^{e_c} π₁₇₉^{e₁₇₉} · R₀².

**(5.5) The norm.** Take norms in the last display. With N(π₂) = −2, N(π₃ₐ) = N(π₃ᵦ) = −3, N(π₃꜀) = 3,
N(π₁₇₉) = 179 and N(u′) = ±1:

    Y² = N(u′)·(−2)^{e₂}·(−3)^{e_a+e_b}·3^{e_c}·179^{e₁₇₉}·N(R₀)².

Y ≠ 0 and N(R₀) ≠ 0 are integers.
- Compare the exponents of 2, 179 and 3 on both sides: e₂ and e₁₇₉ are even, and e_a + e_b + e_c is even.
- Compare signs: N(u′)·(−1)^{e_a+e_b} = 1.

**(5.6) Eight representatives.** By (4.1), u′ = σε^k with σ = ±1, and N(u′) = σ. Write k = 2k₁ + l,
e_a = 2m_a + ν_a, e_b = 2m_b + ν_b and e_c = 2m_c + ν_c, with l, ν_a, ν_b, ν_c ∈ {0, 1}. Then

    X − θ = μ₀S²,    μ₀ = σ ε^l π₃ₐ^{ν_a} π₃ᵦ^{ν_b} π₃꜀^{ν_c},    S = ε^{k₁}π₂^{e₂/2}π₃ₐ^{m_a}π₃ᵦ^{m_b}π₃꜀^{m_c}π₁₇₉^{e₁₇₉/2}R₀.

S ∈ R: ε^{k₁} ∈ R for every integer k₁ because ε⁻¹ ∈ R, and every other exponent is nonnegative. The conditions of
(5.5) become:
- ν_a + ν_b + ν_c is even, so (ν_a, ν_b, ν_c) ∈ {(0,0,0), (1,1,0), (1,0,1), (0,1,1)};
- σ = (−1)^{ν_a+ν_b}, which is +1, +1, −1, −1 respectively.

With l ∈ {0, 1} this leaves eight μ₀ [07: 5.4]:
- the four of (5.1), T := {1, επ₃ₐπ₃ᵦ, −π₃ₐπ₃꜀, −επ₃ᵦπ₃꜀};
- their products with ε, up to the square ε²: ε, π₃ₐπ₃ᵦ, −επ₃ₐπ₃꜀, −π₃ᵦπ₃꜀. This set is called εT below.

**(5.7) The step modulo 4.** S has integer coordinates, so the coordinates of μ₀S² modulo 4 depend only on S modulo 4
[07: 5.5]. X − θ = (X, −6, 0) requires its λ-coordinate to be ≡ 2 and its λ²-coordinate to be ≡ 0 (mod 4).

For the four μ₀ of the coset εT (ε, π₃ₐπ₃ᵦ, −επ₃ₐπ₃꜀, −π₃ᵦπ₃꜀) this never happens [07: 5.5]. Write S = (x, y, z).
- The λ-coordinate of μ₀S² is ≡ ux² (mod 2) with u odd: −5, 3, 1, 1 respectively. So x is even.
- With x even, the pair (λ-coordinate, λ²-coordinate) modulo 4 depends only on (y, z) modulo 2.
- The four possible values are (0, 0), (2, 2), (2, 1) and (0, 3), in some order. None is (2, 0).

For the four μ₀ of (5.1), residues S with the required coordinates exist, as they must, since integral points exist in
each class. Modulo 8 the counts are 0, 0, 0, 0 for εT and 32, 48, 48, 32 for the four μ₀ of (5.1) [07: 5.5].

This proves (5.1). ∎

**(5.8) Hypothesis (H).** By [07: 5.6], with τ₁ = −45 + 48λ − 19λ², τ₂ = 5 + 2λ and τ₃ = −8 + 8λ − 3λ²:
- 9 − θ = επ₃ₐπ₃ᵦ·τ₁²;
- −15 − θ = −π₃ₐπ₃꜀·τ₂²;
- −12 − θ = −επ₃ᵦπ₃꜀·τ₃².

So μ₀ = μ·τ⁻² for the matching μ ∈ {1, 9 − θ, −15 − θ, −12 − θ}, with τ = 1 when μ₀ = 1. Hence X − θ = μ(S/τ)². This is
hypothesis (H) of proof.md §1.2.

The four μ are pairwise distinct modulo K*².
- The quadratic characters at the roots 3 mod 7 and 10 mod 23 take the four different values (1, 1), (−1, 1), (1, −1)
  and (−1, −1) on them [07: 5.6].
- Such a character is constant on a class: if μ = μ′ξ², then η = μ′ξ has η² = μμ′ ∈ R, so η ∈ R by (3.6).

As a confirmation not used by the proof, each of the 28 integral points has an explicit representation
X − θ = μ₀S² with S ∈ R, in the class that the 2-descent assigns to it [07: 5.7].

**(5.9) Remark: what the step modulo 4 is.** It is the local condition at 2 of the 2-descent, here for integral
points with S integral. The norm condition of (5.5) is the condition "norm a square". The local conditions at 3 and 179
are not needed. In particular, no integral point has X ≡ 2 (mod 4): then π₂ divides X − θ exactly three times
[07: 5.3, remark].

## 6. A second proof of Section 3: the box principle

This section proves unique factorization in R again, without the division of Section 2, as an independent check
[08]. It uses:
- polynomial identities and congruences;
- Euclid's algorithm for polynomials over F_p;
- determinants of integer 3×3 matrices;
- a Bernstein bound for an integer cubic form on a box;
- the box principle (Dirichlet's pigeonhole).

It uses no real root of f. The only real number is p^{1/3} in (6.3), and it can be replaced by the rational
u = maxᵢ Bᵢ/tᵢ, which satisfies u³ ≤ p (08 checks Bᵢ³ ≤ tᵢ³p).

**(6.1) Prime lattices.** For β = (a, b, c) write β(t) = a + bt + ct². For a prime p and a monic irreducible factor g of
f modulo p, put

    P(p, g) = {β ∈ R : g divides β(t) in F_p[t]}.

- *A ring map.* u(t)v(t) − (uv)(t) is a multiple of f(t), and g | f mod p. So φ: β ↦ β(t) mod (p, g) respects sums and
  products. It maps R onto the field F_p[t]/(g), which has p^{deg g} elements [08: 6.1(b)].
- *Consequences.* P(p, g) = ker φ has index p^{deg g} in R, and RP(p, g) ⊆ P(p, g). If β ∈ P(p, g), then p | N(β),
  because (N(β), 0, 0) = ββ* ∈ P(p, g).
- *Index lemma.* If π ∈ P(p, g) and |N(π)| = p^{deg g}, then P(p, g) = πR [08: 6.1(c), (d)]. P(p, g) has a triangular
  basis with determinant p^{deg g}, so M(π) = BC with B that basis matrix and C an integer matrix of determinant ±1.

**(6.2) Lemma.** If p | N(β), then β ∈ P(p, g) for some monic irreducible factor g of f mod p.

*Proof.* Otherwise β(t) and f(t) have no common factor in F_p[t]. Euclid's algorithm gives u(t)β(t) + v(t)f(t) ≡ 1, so
U = u(λ) satisfies Uβ ≡ 1 (mod p), and N(U)N(β) ≡ 1 (mod p). ∎

f mod p has a repeated root only for p = 2 (f ≡ t³) and p = 179 (f ≡ (t − 47)²(t − 85)). This follows from
(36t + 81)f + (−12t² − 27t + 32)f′ = 358 = 2·179 [08: 6.2].

**(6.3) Lemma (box principle).** For every prime p and every root r of f mod p, P(p, t − r) contains some β ≠ 0 with
|N(β)| ≤ 22p.

*Proof.*
- Take the basis b = (1, 1 + λ, λ² + λ − 3) of R and the weights (t₀, t₁, t₂) = (2, 1, ½), with t₀t₁t₂ = 1.
- F(x) = N(x₀b₀ + x₁b₁ + x₂b₂) is an integer cubic form. On the box |xᵢ| ≤ tᵢ, cut into four sub-boxes, its Bernstein
  coefficients are at most C₀ = 535/24 < 23 in absolute value. So |F(x)| ≤ C₀ there, and |F(x)| ≤ C₀u³ on
  |xᵢ| ≤ tᵢu, by homogeneity [08: 6.3].
- Put Bᵢ = ⌊tᵢp^{1/3}⌋. Then Π(Bᵢ + 1) > t₀t₁t₂·p = p, so two of the points c with 0 ≤ cᵢ ≤ Bᵢ have the same image
  under φ (into F_p).
- Their difference x gives β = Σxᵢbᵢ ≠ 0 in P(p, t − r), with |xᵢ| ≤ tᵢp^{1/3}. So |N(β)| = |F(x)| ≤ C₀p. ∎

**(6.4) The table.** Generators for every prime p ≤ 22 and for p = 179 [08: 6.4], each checked by the index lemma:

| p | f mod p | generators |
|---|---|---|
| 2 | t³ | π₂ |
| 3 | t(t − 1)(t − 2) | π₃꜀, π₃ₐ, π₃ᵦ |
| 5, 13, 17, 19 | irreducible | p |
| 7 | (t − 3)(t² + 3t + 5) | π₇ = λ² + 2λ − 1, and 7/π₇ = λ² − 4λ + 5 |
| 11 | (t − 6)(t² + 6t + 10) | π₁₁ = −λ² − 3λ − 1, and 11/π₁₁ = 4λ² − 9λ + 7 |
| 179 | (t − 47)²(t − 85) | π₁₇₉ (root 47), and 8λ² + 11λ − 23 (root 85) |

**(6.5) The descent for p > 22.** Let r be a root of f mod p.
- (6.3) gives β ∈ P(p, t − r) with |N(β)| = kp, where k ≤ 22.
- While k > 1, take a prime q | k. Then q < p, and β ∈ P(q, g) = π_{q,g}R by (6.2) and the table, so β = π_{q,g}β′.
- φ_p(π_{q,g}) ≠ 0, since p ∤ N(π_{q,g}). F_p is a field, so β′ ∈ P(p, t − r), with a smaller k.
- When k = 1, the index lemma gives P(p, t − r) = βR.

The other lattices:
- *Quadratic factor.* If f ≡ (t − r)g, then P(p, g) = (p/π_r)R.
- *f irreducible mod p.* P(p, f) = pR.

Only the six degree-one lattices of norm at most 22 can occur in the descent: those above 2, 3 (three) and 7, 11. As a
test, [08: 6.5] runs the descent for every prime 23 ≤ p < 3000. It checks the generators of all 781 prime lattices,
428 of them built by the box-principle descent.

**(6.6) Unique factorization.**
- Every P(p, g) is πR, and such a π is a prime element: π | xy gives φ(x)φ(y) = 0 in a field.
- A nonzero non-unit α has a prime p | N(α). By (6.2), α lies in some P(p, g) = πR, so a prime element divides α.
- Induction on |N| gives a factorization into prime elements. The prime property gives uniqueness up to order and
  units.

With unique factorization, (3.2)–(3.6) and every step of Section 5 follow by comparing exponents of prime elements.
(Bézout itself does not follow this way, but Section 5 does not use it.) So Sections 2 and 6 are two independent
proofs of what Section 5 needs. Only Section 2 is free of the kernels P(p, g); see the Remark below. [08: 6.6] factors X − θ completely for the 14 integral X.

*Remark.* In modern language P(p, g) is a prime ideal of R, and Sections 6.3–6.5 are the classical proof of class
number 1 from a bound on norms. Here the bound comes from the box principle, 22.3, rather than from Minkowski's convex
body theorem, about 7.6. The lattices enter only as the kernels of φ. Section 2 avoids even these kernels, and that
is why it is the main proof.

## 7. What changes in proof.md

- *§1.2.* Hypothesis (H) is now Theorem H. The rank of E(Q), the Mordell–Weil theorem (N5) and the LMFDB value are no
  longer used.
- *§1.3.* The descent script works in Z-span(1, θ/6, θ²/36), which is R, and calls it O_K. It uses these facts about
  Ξ_μ = {ξ ∈ K : μξ² ∈ R}, for the four μ (each in R):
  - *R ⊆ Ξ_μ = −Ξ_μ,* since μ ∈ R.
  - *Closed under addition.* μ(ξ + η)² = μξ² + 2μξη + μη², and μξη ∈ R: indeed (μξη)² = (μξ²)(μη²) ∈ R, so μξη ∈ R
    by (3.6).
  - *Bounded denominators.* If N(μ) | k², then kξ ∈ R for every ξ ∈ Ξ_μ. Indeed, (kξ)² = (k²/N(μ))·μ*·(μξ²) ∈ R, so
    kξ ∈ R by (3.6).

  So R ⊆ Ξ_μ ⊆ (1/k)R, and Ξ_μ is a full lattice. The script's enumeration of (1/k)R/R finds it, and k = 1, 27, 9, 36
  for μ = 1, 9 − θ, −15 − θ, −12 − θ.
  - "O_K" in the script may be read as R. (3.6) replaces N2 and makes N1 unnecessary.
  - The script still prints N1, N2 and N5 as its reasons (for the maximal order, Ξ_μ and completeness). For 1611.a1
    these labels are superseded by (3.6), this section and Theorem H.
- *Theorem B.* It now holds without the rank, by Theorem A, §1 (with the parametrization lemma added there in round 4;
  see review.md) and Theorem H.

The 2-descent itself is unchanged: Theorem H only says in which class X − θ lies.

## 8. Provenance and novelty

Details and reading marks are in [literature_check.md](literature_check.md) §5.
- **The field.** LMFDB 3.1.716.1 lists class number 1 and the fundamental unit λ² − 5λ − 19 = ε. Cavallar and
  Lemmermeyer (1998, §7 table) and Lemmermeyer's survey (updated 2004, p. 28) list discriminant −716 as "N 121/109":
  class number 1, not norm-Euclidean, Euclidean minimum 121/109. Their value rests on a floating-point covering,
  with exact checks at the exceptional points; they do not present it as a proof. Remark 2.5 gives a written proof of
  the lower bound at ξ*; whether their exact check already amounts to one is not recorded. The upper bound
  M(K) ≤ 121/109 is not proved here, and it is not needed.
- **Division with a multiplier.** It is the Dedekind–Hasse criterion: Dedekind (Werke II 373–375, from the Nachlass;
  not accessed, cited by Lemmermeyer), and Hasse, J. reine angew. Math. 159 (1928).
  - Hurwitz (Math. Z. 3, 1919; not accessed, via Lemmermeyer) gives, for every field, an integer j < m with
    |N(jx − y)| < 1, where jx = y is allowed. That is not a Dedekind–Hasse division. Integer multipliers need not
    give one even when the class number is 1: in Z[(1+√−19)/2] no integer s works at ω/2, ω = (1+√−19)/2.
  - Whether a given finite set works is decided by a certificate, including the points ξ with sξ ∈ R. For Z[λ], (2.1)(b)
    shows that all seven half-lattice classes have translates with 0 < |N| < 1.
  - For Z[(1+√−19)/2], hand coverings use the multipliers 1 and 2, and θ or θ̄ at the 2-torsion points: Cámpoli, Amer.
    Math. Monthly 95 (1988), with multipliers 1, 2, 5, θ̄, 2θ; R. A. Wilson's note (2016).
  - Covering certificates for Euclidean minima are standard: Cavallar–Lemmermeyer; Lezowski, Math. Comp. 83 (2014).
    Two-stage divisions are known for other complex cubic fields: −199, −307, −335, −351.
  - No published multiplier certificate was found for −716, or of this kind for any cubic field.
- **The box principle (Section 6)** is Dedekind's own proof of the finiteness of the class number: Supplement XI to
  Dirichlet's Vorlesungen, 4th ed. 1894, §181, Satz I (the box lemma) and Satz II (finiteness). Here it is combined
  with an explicit table, and its constant is sharpened by another basis and a Bernstein bound.
  - For the power basis the constant is about 239 with unit weights (computed by hand), or 171 with optimized weights
    (08, 6.3(c)).
  - Here it is 535/24 ≈ 22.3; the true maximum on the box is ≈ 19.74 (verify_H_box).
- **The units.** The descent on T₂ is the one-unit version of the trace-height descent of mordell_1963 §4 (see the
  provenance on [norm.triples](../../methods/norm.triples.md)).
- **(H) itself.** No earlier elementary treatment of (H) for 1611.a1 was found, and none of the factorization step of
  Avanesov's proof. Avanesov took class number 1 from the Delone–Faddeev table. Stroeker and de Weger (Math. Comp. 68,
  1999, Prop. 1) proved rank 2 with Cremona's mwrank, which they call unconditional (literature_check.md §5.4).
  Uchiyama's
  sister problem (8604.a1) has a norm-Euclidean cubic field (discriminant −239, M = 8/9), but he too took its
  invariants from the table.
- **Novelty.** None is claimed, by the project's criterion. Theorem B's proof remains Avanesov's mechanism: cubic
  factorization, quartic Thue equations, Skolem's method. Here the cubic factorization is made table-free. The
  {1, 2} certificate and the one-unit T₂ descent are presented as certificates, not as new mechanisms.
