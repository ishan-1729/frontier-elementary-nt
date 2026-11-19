# Hypothesis (H) without the rank or Cassels: arithmetic of Z + Zθ + Zω, θ³ = −17

Written 2026-10-15 by Claude (Opus 5.5), revised 2026-10-26 after review round 3 ([review.md](review.md)). It proves
hypothesis (H) of [proof.md](proof.md) §1.2 for integral points,
with elementary means only. With it, Theorem B of proof.md depends neither on the rank of E(Q) nor on Cassels's
theorem. The method is that of the proof for 1611.a1
([triangular_tetrahedral/hypothesis_H.md](../triangular_tetrahedral/hypothesis_H.md),
[norm.triples.multiplier](../../methods/norm.triples.multiplier.md)): the ring is **not** norm-Euclidean (Remark 2.5),
and a division with a multiplier s ∈ {1, 2} (Dedekind–Hasse) replaces Euclid's. Scripts: [scripts/](scripts/), each with
its output beside it as `*.out.txt`. Every check in them is exact. Floating point only seeds some searches, and each
result is then verified in exact arithmetic.

**Theorem H.** Let (X, Y) be an integer solution of

    Y² = X³ + 17                                                       (E)

and let θ³ = −17 and ω = (1 + θ + θ²)/3. Then

    X − θ = μ₀ S²   with S ∈ Z + Zθ + Zω and μ₀ ∈ {1, ρσ, η, ηρσ},

where η = 18 + 7θ, ρ = 1 − θ + ω and σ = 1 + θ + ω. Moreover

    −2 − θ = ρσ,   −1 − θ = η·τ₁²,   4 − θ = ηρσ·τ₂²,

with τ₁ = −11 + 7θ − 6ω and τ₂ = η⁻¹σ = 30 − 19θ + 16ω. So X − θ ∈ μK*² with μ ∈ {1, −2 − θ, −1 − θ, 4 − θ}, which is
hypothesis (H).

**What the proof uses.**
- Integer and rational arithmetic, polynomial identities, congruences, induction (descent), and real inequalities.
- Real numbers enter only through rational interval enclosures of θ₁ = −17^{1/3}, the real root of t³ + 17, obtained
  from exact sign changes. They are used for the units (Section 4) and for Remark 2.5, which Theorem H does not use.
  The division of Section 2 needs only rational triples and a table of 795 boxes, checked in exact rational arithmetic.
- No ideal theory, class number, reciprocity law, unit theorem, Mordell–Weil theorem or rank is used.
- Cassels (1950, §§11–12, 22) proved (H) with ideal factorization in Q(∛17), its class number 1 and the unit η, the
  last two from his Table 2. Here these are replaced by proofs. Sections 2–3 give the divisibility theory that class
  number 1 provided (Bézout, prime elements, square allocation), and Section 4 gives the unit.

**Where (H) was used, and what replaces it.** Every item below is proved here.

| input of proof.md, Theorem B | before | now |
|---|---|---|
| (H), the list of four classes (§1.2) | rank E(Q) ≤ 2 (LMFDB, mwrank) with Mordell–Weil (N5); also Cassels 1950 §22 | Theorem H, Sections 1–5 |
| N2: integral elements form a ring; α² integral implies α integral (§1.3) | standard facts, not proved | (3.6): R is integrally closed |
| N1: O_K = Z + Zθ + Z(1 + θ + θ²)/3 (§1.1) | the index–discriminant formula and p-maximality | (3.6): R is the set of algebraic integers of K |

**Checks.** A reference such as [04: 1.5] means section 1.5 of `04_ring.py`, with the result in `04_ring.out.txt`.
The ring is defined in [hring.py](scripts/hring.py), and [hlib.py](scripts/hlib.py) holds helpers shared with the card
[triangular_square_pyramidal](../triangular_square_pyramidal/hypothesis_H.md).

| Script | Checks |
|---|---|
| [04_ring.py](scripts/04_ring.py) ([output](scripts/04_ring.out.txt)) | ring laws, norm, adjugate, nonvanishing of the norm, the two factors, special elements, real values (Section 1) |
| [05_cover.py](scripts/05_cover.py) ([output](scripts/05_cover.out.txt)) | the box table of the division with multiplier, its cover of D, the half-lattice points, the point ξ* that needs norm 1115/1028 (Section 2) |
| [06_units.py](scripts/06_units.py) ([output](scripts/06_units.out.txt)) | root interval, descent constants, the units of small height (Section 4) |
| [07_classes.py](scripts/07_classes.py) ([output](scripts/07_classes.out.txt)) | the two factors, the norm condition, the identification, the integral points (Section 5) |
| [cover_search.py](scripts/cover_search.py) ([output](scripts/cover_search.out.txt)) | how the box table [dh_boxes.txt](scripts/dh_boxes.txt) was found (not used by the proof) |
| [verify_H_units.py](scripts/verify_H_units.py), [verify_H_cover.py](scripts/verify_H_cover.py), [verify_H_classes.py](scripts/verify_H_classes.py) (with outputs) | independent re-implementations of Sections 1 and 4, 2, and 5 (review.md, round 3) |
| [lit_H_lemmermeyer.py](scripts/lit_H_lemmermeyer.py), [lit_H_lmfdb.py](scripts/lit_H_lmfdb.py) (with outputs) | the identifications with the tables of the literature check, §7 (not used by the proof) |

---

## 1. The ring of integer triples

**(1.1) Multiplication.** A triple (a, b, c) stands for a + bθ + cω, where θ³ = −17 and ω = (1 + θ + θ²)/3. Then
θ² = −1 − θ + 3ω, θω = ω − 6 and ω² = −4 − 2θ + ω, so

    (a,b,c)(d,e,g) = (ad − be − 6bg − 6ce − 4cg,  ae + bd − be − 2cg,  ag + 3be + bg + cd + ce + cg).

This is the product of the corresponding polynomials in θ, reduced modulo f(t) = t³ + 17 [04: 1.1]. The structure
constants are integers, so the integer triples are closed under the product; with coordinatewise addition they form a
commutative ring with unit (1, 0, 0) [04: 1.1]. Write R for the ring of integer triples and R_Q for the ring of
rational triples. θ = (0, 1, 0), and X − θ is the triple (X, −1, 0).

**(1.2) Norm.** Let M(α) be the matrix whose columns are α·1, α·θ, α·ω. For α = (a, b, c),

    M(α) = [[a, −b − 6c, −6b − 4c], [b, a − b, −2c], [c, 3b + c, a + b + c]],

and M(αβ) = M(α)M(β) [04: 1.2]. The norm is the cubic form

    N(α) = det M(α) = a³ + a²c + 17abc + 6ac² − 17b³ − 17b²c + 12c³.

The proof uses only N(αβ) = N(α)N(β), N(rα) = r³N(α) and N(−α) = −N(α) [04: 1.2], and N(X − θ) = X³ + 17 [04: 1.5].
The basis (1, θ, ω) has discriminant det(Tr(bᵢbⱼ)) = −867 = −3·17² [04: 1.2].

**(1.3) Adjugate triple.** For α = (a, b, c) put

    α* = (a² + ac − b² + 5bc + 2c²,  −ab − b² − bc − 2c²,  −ac + 3b² + 2bc).

Then α·α* = (N(α), 0, 0) [04: 1.3]. So α ∈ R is a unit exactly when N(α) = ±1 (1.3a), and a rational triple with
N(α) ≠ 0 has the inverse α*/N(α) (1.3b).

**(1.4) Lemma: the norm vanishes only at 0.** If α ∈ R_Q and N(α) = 0, then α = 0.

*Proof.* Scale a nonzero α to integers with gcd 1. Modulo 7, N(a, b, c) is prime to 7 for each of the 342 nonzero
residue triples [04: 1.4] (it is the norm of F₇[t]/(t³ + 17), a field because −17 ≡ 4 is not a cube mod 7). So
N(a, b, c) = 0 forces a, b, c ≡ 0 (mod 7), against gcd 1. ∎

**(1.5) Special elements** [04: 1.6]:

| element | triple | norm | identity |
|---|---|---|---|
| η = 18 + 7θ | (18, 7, 0) | 1 | 18³ − 17·7³ = 1; η⁻¹ = (275, −175, 147) |
| ρ = 1 − θ + ω | (1, −1, 1) | 3 | 3 = ρ²σ |
| σ = 1 + θ + ω | (1, 1, 1) | 3 | |
| θ | (0, 1, 0) | −17 | θ³ = −17; f′(θ) = 3θ² = ρ²σθ² |

ρ and σ are not associate (neither divides the other) [04: 1.6].

**(1.6) Real values.** f(−3) = −10 < 0 < 9 = f(−2), and f is strictly increasing, so θ₁ = −17^{1/3} ∈ (−3, −2) is its
only real root; bisection encloses it in a rational interval of width below 2⁻¹⁹⁹ [04: 1.7]. Then
f(t) = (t − θ₁)(t² + θ₁t + θ₁²), with complex roots z, z̄, z + z̄ = s := −θ₁ and zz̄ = p := θ₁².

For α = (a, b, c) put σ₁(α) = a + bθ₁ + c(1 + θ₁ + θ₁²)/3 and Q(α) = A² + sAB + pB² = A² − θ₁AB + θ₁²B², where A + Bt
is α, written in powers of θ, reduced modulo t² − st + p:

    A = a + c(1 − θ₁²)/3,    B = b + c(1 − θ₁)/3    [04: 1.7].

Q(α) = |α(z)|². Modulo f(θ₁) = 0 [04: 1.7]: σ₁(α)Q(α) = N(α), and σ₁ and Q are multiplicative.

*Signs.* Q = (A + sB/2)² + (3θ₁²/4)B². So Q ≥ 0 on real triples, where it is only semidefinite (it vanishes where
A = B = 0). Q(α) > 0 for every nonzero rational triple α, since σ₁(α)Q(α) = N(α) ≠ 0 by (1.4).

*Height.* T₂ = σ₁² + 2Q has, modulo f(θ₁), the Gram matrix [[3, 0, 1], [0, 3θ₁², θ₁²], [1, θ₁², (θ₁² − 17θ₁ + 1)/3]],
with leading minors 3, 9θ₁² ≈ 59.50 and 867 [04: 1.7]. So T₂ is positive definite, and only finitely many integer
triples have T₂ below any given bound. (867 = 3·17² = |disc(1, θ, ω)| by (1.2); −867 is the discriminant of K, since
R = O_K by (3.6).)

## 2. Division with a multiplier

**(2.1) Lemma.**
- (a) For every rational triple ξ in D = [−½, ½] × [−½, ½] × [0, ½] there are s ∈ {1, 2} and κ ∈ R with
  |N(sξ − κ)| < 1.
- (b) For every rational triple ξ with 2ξ ∈ R and ξ ∉ R there is κ ∈ R with 0 < |N(ξ − κ)| < 1.

*Proof of (a).* The table [dh_boxes.txt](scripts/dh_boxes.txt) lists 795 closed boxes in D. Each has ends that are
multiples of 1/128, a multiplier s ∈ {1, 2} and a translation κ ∈ R.

*The bound on a box.* For a box with centre ξ_c and half-widths h put v = sξ_c − κ and w = s(ξ − ξ_c), so |wᵢ| ≤ shᵢ.
The exact Taylor expansion N(v + w) = N(v) + Σᵢ gᵢ(v)wᵢ + Σ_{i≤j} q_{ij}(v)wᵢwⱼ + N(w) [05: 2.1], with, for
v = (x, y, z),

    g₁ = 3x² + 2xz + 17yz + 6z²,   g₂ = 17xz − 51y² − 34yz,   g₃ = x² + 17xy + 12xz − 17y² + 36z²,
    q₁₁ = 3x + z,  q₁₂ = 17z,  q₁₃ = 2x + 17y + 12z,  q₂₂ = −51y − 17z,  q₂₃ = 17x − 34y,  q₃₃ = 6x + 36z,

gives |N(sξ − κ)| ≤ B on the box, where

    B = |N(v)| + Σᵢ |gᵢ(v)|shᵢ + Σ_{i≤j} |q_{ij}(v)|s²hᵢhⱼ + Σ_m |n_m|(sh)^m.

Here n_m is the coefficient of N at the exponent triple m, (sh)^m = Πᵢ(shᵢ)^{mᵢ}, and Σ|n_m| = 71. B is computed
exactly, and B < 1 for every box; the largest is 131071/131072 [05: 2.2].

*The boxes cover D.* Each of the 128·128·64 = 1048576 cells of side 1/128 that tile D lies in exactly one box, so every
point of D lies in a closed cell, hence in a box [05: 2.3]. The boxes are also exactly the leaves of a binary tree of
midpoint splits of D, and their volumes add up to ½. So every point of D has the listed (s, κ) with |N(sξ − κ)| < 1. ∎

*Proof of (b).* Modulo R there are seven such ξ, the triples x/2 with x ∈ {0, 1}³ \ {0}. With the listed κ, the values
of N(ξ − κ) are [05: 2.4]:

| ξ | (0,0,½) | (0,½,0) | (0,½,½) | (½,0,0) | (½,0,½) | (½,½,0) | (½,½,½) |
|---|---|---|---|---|---|---|---|
| κ | (1,0,0) | (−3,0,1) | (−1,1,0) | (0,0,0) | (0,0,1) | (0,−1,−2) | (0,0,0) |
| N(ξ − κ) | −½ | ⅝ | ¼ | ⅛ | −¾ | ¼ | ⅜ |

All lie strictly between 0 and 1 in absolute value. Any other ξ with 2ξ ∈ R and ξ ∉ R is ξ = x/2 + r with one of these
x and some r ∈ R; take κ_x + r, so that N(ξ − κ_x − r) = N(x/2 − κ_x). ∎

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
from (2.2) and put ρ′ = sα − κβ = (sξ − κ)β ∈ R. Then |N(ρ′)| = |N(sξ − κ)|·|N(β)|, and 0 < |N(sξ − κ)| < 1. ∎

This is the Dedekind–Hasse criterion with the multipliers {1, 2}.

**(2.4) Remark: the table.**
- *Shape.* 560 boxes use s = 1 and 235 use s = 2; the deepest is reached by 19 midpoint splits; the largest |κᵢ| is 104
  [05: 2.2].
- *How it was found.* The search ([cover_search.py](scripts/cover_search.py)) tried, on each box with centre ξ_c, the
  translates κ = η^{−k}κ′ with κ′ an integer triple near η^k·sξ_c and k ∈ {−1, 0, 1}, with two split rules. The proof
  uses only the table.
- *Sharp bound.* The Taylor bounds are not sharp. The largest Bernstein coefficient over all boxes is 1046475/1048576.
  It is a corner value, at ξ = (−1/8, 49/128, 33/128), in the box YYZYxzyXzyzXyxzZy with s = 1 and κ = (2, 1, 0). So the
  largest value of |N(sξ − κ)| over the boxes, each with its own (s, κ) on its own closed box, is exactly
  1046475/1048576 [05: 2.6]. The corner lies in five other boxes, all with smaller values [05: 2.6].
- *Confirmation.* 3000 random pairs (α, β), none with β | α, were divided by (2.3), with the expected inequality every
  time [05: 2.7].

**(2.5) Remark: R is not norm-Euclidean.** Theorem H does not use this remark; it shows that the multiplier cannot be
dropped. Let ξ* = (94/257, 233/514, 19/1028). Then −ηξ* − ξ* = (−3, −8, −10) and −η⁻¹ξ* − ξ* = (211, −135, 113) are in
R, so −η and its inverse −η⁻¹ map the coset C = ξ* + R onto itself. (η maps C onto −C, where |N| takes the same values:
ηξ* + ξ* = (3, 8, 10) ∈ R.) Its element ξ* + (−1, −2, −2) = (−163/257, −795/514, −2037/1028) has |N| = 1115/1028, and
no element of C has smaller |N| [05: 2.5]. So no κ ∈ R has |N(ξ* − κ)| < 1. (The multiplier 2 helps: 2ξ* has a
translate with |N| = 3/514.)

*Proof of the bound.* As in [triangular_square_pyramidal/hypothesis_H.md](../triangular_square_pyramidal/hypothesis_H.md)
(2.5), with the unit −η in place of ε (the sign does not change T₂). Suppose some y ∈ C has |N(y)| < 1115/1028.
- The elements (−η)^k y with T₂ ≤ T₂(y) lie in C and are finitely many, so the orbit has an element y₀ of least T₂, and
  N(y₀) = ±N(y).
- Neither (−η)y₀ nor (−η)⁻¹y₀ has smaller T₂, so T₂(y₀)³ ≤ K₃·N(y)², with the exact threshold
  K₃ = 4(E² + E + 1)³/(E²(E + 1)²) of (4.3′), E = Q(η) as in (4.2); K₃ ≤ 3782604 [05: 2.5.2]. The rounded C of (4.3)
  would not fit here: C³·(1115/1028)² > T_e³ [05: 2.5.3].
- So T₂(y₀)³ < K₃·(1115/1028)² ≤ T_e³, with T_e := 82241/500.
- [05: 2.5.3] lists the 303 elements y ∈ C with T_low(y) ≤ T_e, where T_low = yᵀ(G′ − 3ρI)y ≤ T₂ is built as in the
  sister card (2.5), with ρ = 2⁻²⁰. The list contains every y ∈ C with T₂(y) ≤ T_e. Its least |N| is 1115/1028,
  attained only at ξ* + (−1, −2, −2).
- y₀ is on that list and has |N(y₀)| < 1115/1028, a contradiction. ∎

*The other fixed points.* There are |N(η − 1)| = 918 classes fixed by η modulo R, and |N(−η − 1)| = 1028 fixed by −η.
Only the classes of ±ξ* have no element of norm below 1; the next largest least norm among the fixed classes is
293/306 [05: 2.5.4]. So R has Euclidean minimum at least 1115/1028, the value that Lemmermeyer and Cavallar–Lemmermeyer
record (§7).

## 3. Common divisors, primes and square allocation

Exactly as in [triangular_square_pyramidal/hypothesis_H.md](../triangular_square_pyramidal/hypothesis_H.md) §3, with R
this ring and (2.3) giving s ∈ {1, 2}: in the proof of (3.1), δ ∤ α gives s, κ with 0 < |N(sα − κδ)| < |N(δ)|, and
sα − κδ ∈ L. The statements:
- **(3.1)** Gcd and Bézout: a gcd δ = xα + yβ exists.
- **(3.2)** Euclid's lemma; **(3.2′)** products of elements coprime to γ are coprime to γ.
- **(3.3)** An element of prime norm ±q is a prime element. This applies to ρ, σ (q = 3) and θ (q = 17).
- **(3.4)** Exact powers of a prime element exist.
- **(3.5)** Coprime square allocation: αβ = uW² with α, β coprime gives α = u′R₀², u′ a unit.
- **(3.6)** R is integrally closed, and every element of R is integral (characteristic polynomial [04: 1.2]); so
  O_K = R, and ξ ∈ K with ξ² ∈ R lies in R.

## 4. The units

**(4.1) Theorem.** The units of R are exactly ±η^k (k ∈ Z), with k determined by the unit. Every unit is σ′η^l v² with
σ′ ∈ {1, −1}, l ∈ {0, 1} and v a unit.

**(4.2) Values of η** [06: 4.2].
- e := σ₁(η) ∈ (0.0010288653, 0.0010288654), so 0 < e < 1.
- E := Q(η) = σ₁(η⁻¹) = 1/e, enclosed directly in (971.9444, 971.9445).
- So σ₁(η⁻¹) = E and Q(η⁻¹) = e.

**(4.3) Descent lemma.** As in [triangular_square_pyramidal/hypothesis_H.md](../triangular_square_pyramidal/hypothesis_H.md)
(4.3), with η in place of ε. For a nonzero rational triple y put x = σ₁(y)² > 0 and q = Q(y) > 0. Then T₂(yη) < T₂(y)
exactly when x/q > ρ_hi := 2E²/(E + 1), and T₂(yη⁻¹) < T₂(y) exactly when x/q < ρ_lo := 2/(E(E + 1)) [06: 4.3]. ρ_hi
increases and ρ_lo decreases with E, so E < 971.9445 gives ρ_hi < 1941.8911 and ρ_lo > 2.114951·10⁻⁶ [06: 4.3].

*Lemma.* Let u be a unit with T₂(uη) ≥ T₂(u) and T₂(uη⁻¹) ≥ T₂(u). Then T₂(u) ≤ C := max(h(A), h(B)), where
h(q) := q⁻² + 2q, A = 80153/10⁶ and B = 77905221/10⁶. Here C = h(A) = 155.81436… ≤ 155.8144 [06: 4.3].

*Proof.* Here xq² = 1, so q³ = (x/q)⁻¹ ∈ (1/1941.8911, 1/(2.114951·10⁻⁶)). Since A³·1941.8911 ≤ 1 ≤ B³·2.114951·10⁻⁶,
A < q < B [06: 4.3]. T₂(u) = h(q) ≤ max(h(A), h(B)), by the endpoint identities for h in the sister card (4.3)
[06: 4.3]. ∎

**(4.3′) Scaled form.** Let y be a nonzero rational triple with |N(y)| = n. If T₂(yη) ≥ T₂(y) and T₂(yη⁻¹) ≥ T₂(y),
then T₂(y) ≤ C·n^{2/3}, as in the sister card. Its exact threshold, used in Remark 2.5, is T₂(y)³ ≤ K₃·N(y)² with
K₃ = 4(E² + E + 1)³/(E²(E + 1)²) = 3782603.77… ≤ 3782604, the largest value of ψ(t) = (t + 2)³/t on [ρ_lo, ρ_hi]
[05: 2.5.2].

**(4.4) The units of small height.** 06 lists the integer triples with T_low ≤ C, where T_low = yᵀ(G′ − 3ρI)y with G′
a rational matrix within ρ = 2⁻³⁰ of the Gram matrix of T₂, entrywise. T_low ≤ T₂, by the inequality of the sister card
(2.5), so the list contains every triple with T₂ ≤ C. It has 275 triples, and the only units among them are ±1
[06: 4.4].

**(4.5) Proof of (4.1).** Let u be a unit.
- The set {uη^k : T₂(uη^k) ≤ T₂(u)} is finite, since T₂ is positive definite. It is nonempty.
- Let v = uη^j be an element of it of least T₂. Then T₂(vη) ≥ T₂(v) and T₂(vη⁻¹) ≥ T₂(v), so T₂(v) ≤ C by (4.3).
- So v = ±1 by (4.4), and u = ±η^{−j}.
- *Uniqueness.* 0 < σ₁(η) < 1, so σ₁(η^k) = σ₁(η)^k is not ±1 unless k = 0.
- Finally write k = 2k₁ + l with l ∈ {0, 1}: then ±η^k = ±η^l(η^{k₁})². ∎

**(4.6) Remark.** By (4.1), η is a fundamental unit. In addition, ±1 and ±η are distinct modulo squares of units: a
square has σ₁ > 0, and a square root v of η would have T₂(v) = e + 2√E ≤ 62.36 < C, so v = ±1 by (4.4), which is false
[06: 4.5]. Cassels's Table 2 (p. 270) lists, for D = 17, the class number 1 and the unit 18 − 7δ with δ = ∛17 = −θ,
which is η. The entry has no star, which by his note 2 (p. 271) marks a fundamental unit proved so by Dedekind or by
Nagell; his §11 (p. 249) needs only an odd power of the fundamental unit.

## 5. The class of X − θ

**(5.1) Proposition (Theorem H).** Let (X, Y) be an integer solution of (E). Then X − θ = μ₀S² with S ∈ R and
μ₀ ∈ {1, ρσ, η, ηρσ}.

**(5.2) Two factors.** Put α = X − θ = (X, −1, 0) and β = X² + θX + θ² = (X² − 1, X − 1, 3). Then [07: 5.1]

    αβ = X³ + 17 = Y²,    N(α) = X³ + 17 = Y²,    β − (X + 2θ)α = 3θ² = ρ²σθ².

Since α ≠ 0, N(α) ≠ 0 by (1.4), so Y ≠ 0 and β ≠ 0.

**(5.3) Common divisors.** Let 𝒮 = {ρ, σ, θ}. These are prime elements by (3.3), and they are pairwise non-associate:
ρ and σ do not divide each other [07: 5.2], and θ has another norm. A common divisor of α and β divides
β − (X + 2θ)α = ρ²σθ².

**(5.4) Separating 𝒮.** The argument of
[triangular_square_pyramidal/hypothesis_H.md](../triangular_square_pyramidal/hypothesis_H.md) (5.4), with 𝒮 = {ρ, σ, θ}
and ρ²σθ² in place of (π₃ₐπ₃ᵦπ₃꜀)²π₂₃₉:
- By (3.4), applied to each ϖ ∈ 𝒮 in turn, write α = Π_{ϖ∈𝒮} ϖ^{e_ϖ}·α′, β = Π_{ϖ∈𝒮} ϖ^{f_ϖ}·β′ and
  Y = Π_{ϖ∈𝒮} ϖ^{w_ϖ}·Y′, where no element of 𝒮 divides α′, β′ or Y′ (each later quotient divides the earlier ones).
- *Matching exponents.* For each ϖ ∈ 𝒮, αβ = Y² reads ϖ^{e+f}A = ϖ^{2w}B with A = α′β′·Π_{ϖ′≠ϖ}ϖ′^{e′+f′} and
  B = Y′²·Π_{ϖ′≠ϖ}ϖ′^{2w′}, neither divisible by ϖ (3.3), (5.3). Cancelling the smaller power of ϖ (1.4) shows
  e_ϖ + f_ϖ = 2w_ϖ. So α′β′ = Y′².
- *Coprime.* A common divisor δ of α′ and β′ divides ρ²σθ² (5.3). For each ϖ ∈ 𝒮 a gcd of δ and ϖ is a unit or an
  associate of ϖ (3.3), and an associate would give ϖ | α′. So δ is coprime to every ϖ ∈ 𝒮, hence to ρ²σθ² by (3.2′),
  and δ | 1.
- *Allocation.* (3.5) with u = 1 and W = Y′ gives α′ = u′R₀², with u′ a unit. Thus, with e_ρ, e_σ, e₁₇ the exponents
  of ρ, σ, θ,

      α = u′ · ρ^{e_ρ} σ^{e_σ} θ^{e₁₇} · R₀².

**(5.5) The norm.** Take norms. With N(ρ) = N(σ) = 3, N(θ) = −17 and N(u′) = ±1:

    Y² = N(u′)·3^{e_ρ+e_σ}·(−17)^{e₁₇}·N(R₀)².

Y ≠ 0 and N(R₀) ≠ 0 (since α′ ≠ 0) are integers. Comparing the exponents of 17 and 3, e₁₇ and e_ρ + e_σ are even;
then comparing signs, N(u′) = 1.

**(5.6) Four representatives.** By (4.1), u′ = σ′η^k with σ′ = ±1, and N(u′) = σ′; so σ′ = 1 by (5.5). Write
k = 2k₁ + l, e_ρ = 2m_ρ + ν and e_σ = 2m_σ + ν, with l, ν ∈ {0, 1} (one ν, since e_ρ and e_σ have the same parity).
Then

    X − θ = μ₀S²,    μ₀ = η^l(ρσ)^ν,    S = η^{k₁}ρ^{m_ρ}σ^{m_σ}θ^{e₁₇/2}R₀.

S ∈ R: η^{k₁} ∈ R for every integer k₁ because η⁻¹ = (275, −175, 147) ∈ R (1.5), and every other exponent is
nonnegative. So μ₀ ∈ {1, ρσ, η, ηρσ} [07: 5.3]. This proves (5.1). ∎

No further condition, modulo 4 or otherwise, can remove a candidate: all four occur, X = 8, −2, −1, 4 giving
μ₀ = 1, ρσ, η, ηρσ [07: 5.5].

**(5.7) Hypothesis (H).** By [07: 5.4], −2 − θ = ρσ, −1 − θ = η·τ₁² and 4 − θ = ηρσ·τ₂². So μ₀ = μ·τ⁻² for the matching
μ ∈ {1, −2 − θ, −1 − θ, 4 − θ} and τ ∈ {1, 1, τ₁, τ₂}, and X − θ = μ(S/τ)². This is hypothesis (H) of proof.md §1.2.

The four μ are pairwise distinct modulo K*².
- *Characters.* For a prime p ≠ 3 and an integer r with r³ ≡ −17 (mod p), a + bθ + cω ↦ a + br + c(1 + r + r²)/3
  (mod p) is a ring map R → F_p. It is the restriction to R of the ring map from the polynomials in θ with coefficients
  in Z[1/3] to F_p that sends θ to r; that map exists because 3 is invertible mod p and r³ ≡ −17. So
  χ_{p,r}(α) = ((a + br + c(1 + r + r²)/3)/p), a Legendre symbol, is multiplicative on R.
- The characters at the roots 2 mod 5 and 3 mod 11 take the values (1, 1), (1, −1), (−1, −1), (−1, 1) on 1, −2 − θ,
  −1 − θ, 4 − θ [07: 5.4]. All these values are nonzero, and the four pairs differ.
- If μ = μ′ξ² with ξ ∈ K, then ζ = μ′ξ has ζ² = μμ′ ∈ R, so ζ ∈ R by (3.6), and χ(μ)χ(μ′) = χ(ζ)². When χ(μ) and
  χ(μ′) are nonzero, the left side is nonzero, so it is 1 and χ(μ) = χ(μ′).

As a confirmation not used by the proof, each of the 16 integral points has X − θ = μ₀S² with S ∈ R for exactly one
candidate μ₀, in the class that the 2-descent assigns to it [07: 5.5].

**(5.8) Remark.** The norm condition of (5.5) is the condition "norm a square" of the 2-descent; here it alone leaves
the four classes, so no local condition at 2, 3 or 17 is needed. This is Cassels's §22 "a priori" list μ = 1, η, λ, λη,
here without ideals: his λ and λη are ρσ and ηρσ in some order, since his definition fixes λ only up to the unit η.

## 6. What changes in proof.md

- *§1.2.* Hypothesis (H) is now Theorem H. Neither the rank of E(Q) nor Cassels's theorem is used.
- *§1.3.* The descent works in O_K = Z + Zθ + Zω, which is R. For the four μ (each in R), Ξ_μ = {ξ ∈ K : μξ² ∈ R} is a
  lattice, closed under addition, with R ⊆ Ξ_μ ⊆ (1/k)R whenever N(μ) | k², by (3.6), exactly as in the sister card §6.
  Here N(1) = 1, N(−2 − θ) = 9, N(4 − θ) = 81 and N(−1 − θ) = 16, so k = 1, 3, 9, 4 respectively. (3.6) replaces N2 and
  makes N1 unnecessary.
- *Theorem B.* It now holds without the rank or Cassels, by Theorem A, proof.md §1 (with the parametrization lemma
  1.3′), proof.md §5 (the reconstruction of the points from the orbits) and Theorem H.

## 7. Provenance and novelty

Details and reading marks are in [literature_check.md](literature_check.md) §7.
- **The field.** K = Q(∛17) has discriminant −867 = −3·17² [04: 1.2]. Class number 1 and the unit η are classical:
  Cassels 1950, Table 2; Selmer's tables (1956); LMFDB 3.1.867.1, whose fundamental unit is −η.
- **Not norm-Euclidean.** Cioffari (1979) showed that Q(∛d) is norm-Euclidean only for d = 2, 3, 10 (abstract read).
  Lemmermeyer (1989) gives, for this field, the point −ξ* with M(K, −ξ*) = 1115/1028, and Cavallar and Lemmermeyer list
  −867 as "N 1115/1028" (arXiv:1202.6019v1, §7 table, p. 17, read in review round 3). Remark 2.5 re-proves the lower
  bound 1115/1028 at that point, with an exact finiteness argument. The upper bound M(K) ≤ 1115/1028 is a computation of
  theirs; it is not proved here, and it is not needed.
- **Division with a multiplier.** The Dedekind–Hasse criterion with multipliers {1, 2}, as for 1611.a1; see the
  provenance there ([norm.triples.multiplier](../../methods/norm.triples.multiplier.md)). No published division
  algorithm of any kind for this ring was found: no multiplier certificate, two-stage division or weighted norm.
- **The units.** The one-unit descent on T₂, as for 1611.a1 and 8604.a1.
- **(H) itself.** Cassels's §22 proves it, for all rational points, with ideals, class number 1 and the tabulated η;
  here it is table-free and ideal-free for integral points. No earlier ideal-free proof was found.
- **Novelty.** None is claimed, by the project's criterion. The descent is Cassels's, written for elements; the table
  and the unit descent are certificates for classical facts.
