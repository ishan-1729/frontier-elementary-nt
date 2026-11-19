# Mordell's equation y(y+1) = x(x+1)(x+2) without tables

**Theorem.** The integer solutions of

    y(y+1) = x(x+1)(x+2)                                            (1)

are exactly: x ∈ {0, −1, −2} with y ∈ {0, −1}; x = 1 with y ∈ {2, −3}; x = 5 with y ∈ {14, −15}.

This is L. J. Mordell, Pacific J. Math. 13 (1963) 1347–1351. Mordell works in the cubic field of θ³ − 4θ + 2 = 0
and takes three facts from the table of Delaunay and Faddeev: the integers of the field are a + bθ + cθ², the
class number is 1, and θ − 1, 2θ − 1 are fundamental units. Here those facts are replaced by proofs:
- a Euclidean division on integer triples, with 24 explicit rounding boxes checked in exact rational arithmetic
  (Section 2), which gives gcds, Bézout identities and square allocation (Section 3);
- a descent on the positive integer height T(u) = Tr(u²), which shows that every unit is ±(θ−1)^r(2θ−1)^s
  (Section 4).

The fact "the integers are a + bθ + cθ²" is not needed: all arithmetic stays inside the ring of integer
triples (Remark 3.7 proves it anyway). Mordell's case analysis of simultaneous quadratics is redone in full
in Section 6.

**What the proof uses:** integer and rational arithmetic, polynomial identities (including Girard's
identities for symmetric functions of three numbers), congruences, induction (descent), and real inequalities.
Section 3 builds the divisibility theory of the ring of triples (gcd, Bézout, primes, square allocation) with
every lemma proved. Real numbers are used only for the units (Section 4), through the three real roots of
t³ − 4t + 2 introduced in (1.6); they exist by the intermediate value theorem (bisection) applied to sign
changes. The Euclidean division of Section 2 needs only rational triples. No ideal theory, class-number theorem,
reciprocity law, unit theorem or table of invariants is used.

**Checks.** Every identity, bound and finite enumeration below is checked by a script in
[scripts/](scripts/), whose output is stored next to it. A reference such as [01: 1.5] means section 1.5 of
`01_triples.py`, with the result in `01_triples.out.txt`.

| Script | Checks |
|---|---|
| [01_triples.py](scripts/01_triples.py) ([output](scripts/01_triples.out.txt)) | ring laws, norm, adjugate, nonvanishing of the norm, special elements, real values (Section 1) |
| [02_euclid_cover.py](scripts/02_euclid_cover.py) ([output](scripts/02_euclid_cover.out.txt)) | the 24 boxes and their bounds, the cover of D (Section 2) |
| [03_units.py](scripts/03_units.py) ([output](scripts/03_units.out.txt)) | root intervals, descent constants, the 20 small units, independence (Section 4) |
| [04_reduction.py](scripts/04_reduction.py) ([output](scripts/04_reduction.out.txt)) | the substitution, parity, the factors A and B, their common divisors, the norm and mod-8 step (Sections 0 and 5) |
| [05_cases.py](scripts/05_cases.py) ([output](scripts/05_cases.out.txt)) | the four cases (Section 6) and the final list (Section 7) |
| [06_crosscheck.py](scripts/06_crosscheck.py) ([output](scripts/06_crosscheck.out.txt)) | confirmations not used by the proof: a direct search to \|X\| ≤ 2·10⁶, explicit representations, 20000 random divisions |
| [cover_search.py](scripts/cover_search.py) ([output](scripts/cover_search.out.txt)) | how the box table was found (not used by the proof) |

The shared definitions are in [triples.py](scripts/triples.py). Run each script from the `scripts/` folder with
Python 3 and sympy (checked with Python 3.13.5, sympy 1.14.0).

---

## 0. Reduction to 2Y² = X³ − 4X + 2

**(0.1)** Put Y = 2y + 1 and X = 2x + 2. Then

    2Y² − (X³ − 4X + 2) = 8[y(y+1) − x(x+1)(x+2)]                  [04: 5.1]

so (x, y) solves (1) exactly when (X, Y) solves

    2Y² = X³ − 4X + 2.                                              (E)

**(0.2)** Let (X, Y) be an integer solution of (E). If X were odd, X³ − 4X + 2 would be odd, but 2Y² is even.
So X = 2X₁ is even, and then Y² = 4X₁³ − 4X₁ + 1 [04: 5.2]. Since X₁³ − X₁ = (X₁−1)X₁(X₁+1) is even,
Y² ≡ 1 (mod 8); in particular Y is odd. Conversely X even and Y odd give integers x = X/2 − 1 and
y = (Y − 1)/2. So it suffices to find all integer solutions of (E).

## 1. The ring of integer triples

**(1.1) Multiplication.** A triple (a, b, c) stands for a + bθ + cθ², where θ³ = 4θ − 2 (so θ⁴ = 4θ² − 2θ).
Define

    (a,b,c)(d,e,f) = (ad − 2bf − 2ce,  ae + bd + 4bf + 4ce − 2cf,  af + be + cd + 4cf).

This is the product of a + bt + ct² and d + et + ft², reduced modulo t³ − 4t + 2 [01: 1.1]. With
coordinatewise addition the triples form a commutative ring with unit (1, 0, 0) [01: 1.2]. Write 𝓡 for the
ring of integer triples and 𝓡_Q for the ring of rational triples; an integer r is the triple (r, 0, 0). We
write θ = (0, 1, 0), and a triple (a, b, c) also as a + bθ + cθ².

**(1.2) Norm.** Let M(α) be the 3×3 matrix whose columns are the triples α·1, α·θ, α·θ². For α = (a,b,c),

    M(α) = [[a, −2c, −2b], [b, a+4c, 4b−2c], [c, b, a+4c]],

and M(αβ) = M(α)M(β) [01: 1.3]. Define the norm as the cubic form

    N(α) = a³ + 8a²c − 4ab² + 6abc + 16ac² − 2b³ + 8bc² + 4c³,

which is the determinant of M(α) [01: 1.3]. The determinant only explains where N comes from: what we use are
the polynomial identities N(αβ) = N(α)N(β), N(rα) = r³N(α) for an integer or rational r, and N(−α) = −N(α),
checked by expanding both sides [01: 1.4]. N(α) is an integer for α ∈ 𝓡. The norm is defined without
conjugates; the real values of (1.6), which are the conjugates, are used only in Section 4.

**(1.3) Adjugate triple.** For α = (a, b, c) put

    α* = (a² + 8ac − 4b² + 2bc + 16c²,  −ab − 2c²,  b² − ac − 4c²).

Then α·α* = (N(α), 0, 0), and α* is the first column of the adjugate matrix of M(α) [01: 1.5]. Consequences:
- **(1.3a)** α ∈ 𝓡 is a unit of 𝓡 exactly when N(α) = ±1. If αβ = 1 then N(α)N(β) = 1 with both integers; if
  N(α) = ±1 then α⁻¹ = N(α)·α* ∈ 𝓡.
- **(1.3b)** A rational triple α with N(α) ≠ 0 has the inverse α*/N(α) in 𝓡_Q.

**(1.4) Lemma: the norm vanishes only at 0.** If α ∈ 𝓡_Q and N(α) = 0, then α = 0.

*Proof.* Suppose α ≠ 0. Multiply α by a common denominator (N(rα) = r³N(α)) and divide by the gcd of the
coordinates. This gives integers (a, b, c) with gcd 1 and N(a, b, c) = 0. We show a, b, c are all even, a contradiction.
- N(a, b, c) − a³ has even coefficients [01: 1.6], so a³ is even: a = 2a′.
- N(2a′, b, c) = 2·P with P an integer polynomial and P + b³ having even coefficients [01: 1.6]. Since
  N = 0, P = 0, so b³ is even: b = 2b′.
- N(2a′, 2b′, c) = 4·P′ with P′ − c³ having even coefficients [01: 1.6]. So c³ is even and c is even. ∎

Consequently 𝓡 has cancellation: if αβ = αγ and α ≠ 0 then β = γ (multiply by α*, then divide by the nonzero
integer N(α)). For α, β ∈ 𝓡 we write β | α if α = βγ with γ ∈ 𝓡; then N(β) divides N(α). Two elements are
*associates* if each divides the other, that is, if they differ by a unit factor.

**(1.5) Special elements.** Let

    ε = θ − 1 = (−1, 1, 0),   η = 2θ − 1 = (−1, 2, 0),   π = 4θ − 3 = (−3, 4, 0).

Then [01: 1.7]:
- N(θ) = −2, N(ε) = 1, N(η) = −1, N(π) = 37;
- θ³ = 2η, so 2 = θ³·η⁻¹;
- ε⁻¹ = (−3, 1, 1) and η⁻¹ = (15, −2, −4).

Mordell's paper prints N(4θ − 3) = −37. The exact value is +37 = −64·f(3/4) with f(t) = t³ − 4t + 2; only
|N(π)| = 37 matters below.

**(1.6) Real values.** The polynomial f(t) = t³ − 4t + 2 changes sign on [−3, −2], [0, 1] and [1, 2]
(f(−3) = −13, f(−2) = 2, f(0) = 2, f(1) = −1, f(2) = 2) [01: 1.8]. By the intermediate value theorem (repeated
bisection) there is a root in each of the open intervals (−3, −2), (0, 1), (1, 2); since f(1) ≠ 0 they are three
distinct roots θ₁ < θ₂ < θ₃, with θ₁ ≈ −2.21432, θ₂ ≈ 0.53919, θ₃ ≈ 1.67513, each enclosed in a rational
interval of width below 2⁻⁵⁰ [03: 4.1]. By the factor theorem f(t) = (t − θ₁)(t − θ₂)(t − θ₃), so the
elementary symmetric functions are

    e₁ = θ₁ + θ₂ + θ₃ = 0,   e₂ = θ₁θ₂ + θ₁θ₃ + θ₂θ₃ = −4,   e₃ = θ₁θ₂θ₃ = −2.

For a triple α = (a, b, c) put αᵢ = a + bθᵢ + cθᵢ². Because θᵢ³ = 4θᵢ − 2, the map α ↦ αᵢ respects sums and
products: (αβ)ᵢ = αᵢβᵢ. We need two identities [01: 1.8]:

    α₁α₂α₃ = N(α),          α₁² + α₂² + α₃² = T(α) := 3a² + 16ac + 8b² − 12bc + 32c².

*Sum of squares.* The power sums pⱼ = θ₁ʲ + θ₂ʲ + θ₃ʲ are p₀ = 3, p₁ = e₁ = 0, p₂ = e₁² − 2e₂ = 8, and
pⱼ = 4pⱼ₋₂ − 2pⱼ₋₃ for j ≥ 3 (sum θᵢʲ = 4θᵢʲ⁻² − 2θᵢʲ⁻³ over i), so p₃ = −6 and p₄ = 32. Then
Σαᵢ² = a²p₀ + 2abp₁ + (b² + 2ac)p₂ + 2bcp₃ + c²p₄ = T(α).

*Product.* Expanding Π(a + bθᵢ + cθᵢ²), the coefficient of each monomial in a, b, c is a symmetric sum of the
θᵢ, evaluated by Girard's identities Σθᵢ² = e₁² − 2e₂, Σ_{i≠j} θᵢθⱼ² = e₁e₂ − 3e₃, Σ_{i<j} θᵢ²θⱼ² = e₂² − 2e₁e₃:

| monomial | a³ | a²b | a²c | ab² | abc | ac² | b³ | b²c | bc² | c³ |
|---|---|---|---|---|---|---|---|---|---|---|
| coefficient | 1 | e₁ | Σθᵢ² | e₂ | Σ_{i≠j}θᵢθⱼ² | Σ_{i<j}θᵢ²θⱼ² | e₃ | e₃e₁ | e₃e₂ | e₃² |
| value | 1 | 0 | 8 | −4 | 6 | 16 | −2 | 0 | 8 | 4 |

The values are the coefficients of N, so α₁α₂α₃ = N(α).

Writing Tr(β) = β₁ + β₂ + β₃ = ap₀ + bp₁ + cp₂ = 3a + 8c for β = (a, b, c) (also the trace of M(β)), we have
T(α) = Tr(α²), and
we call T the *trace height* [01: 1.8]. T is a positive definite quadratic form:
T = 3(a + 8c/3)² + 8(b − 3c/4)² + (37/6)c² [03: 4.3], with leading minors 3, 24, 148 [01: 1.8]. For α ∈ 𝓡,
T(α) is a nonnegative integer.

## 2. Euclidean division

**(2.1) Lemma.** For every rational triple ξ there is an integer triple κ with |N(ξ − κ)| < 1.

*Proof.* Let r be the integer triple of nearest integers to the coordinates of ξ (either choice at a tie). If
κ′ works for ξ − r, then r + κ′ works for ξ, since ξ − (r + κ′) = (ξ − r) − κ′. So we may assume that
ξ = (x, y, z) has x, y, z ∈ [−½, ½]. If z < 0, replace ξ by −ξ: since N(−ξ − κ) = −N(ξ + κ), a translation κ
that works for −ξ gives the translation −κ for ξ. So it suffices to treat ξ in

    D = [−½, ½] × [−½, ½] × [0, ½].

*The bound on a box.* Let a box be [x₀ − h₁, x₀ + h₁] × [y₀ − h₂, y₀ + h₂] × [z₀ − h₃, z₀ + h₃], let κ be an
integer triple, and put v = (x₀, y₀, z₀) − κ. For |uᵢ| ≤ hᵢ, the exact Taylor expansion of the cubic N is
[02: 2.1]

    N(v + u) = N(v) + g₁(v)u₁ + g₂(v)u₂ + g₃(v)u₃ + Σ_{i≤j} q_{ij}(v)uᵢuⱼ + N(u),

where, for v = (x, y, z),

    g₁ = 3x² + 16xz − 4y² + 6yz + 16z²,   g₂ = −8xy + 6xz − 6y² + 8z²,   g₃ = 8x² + 6xy + 32xz + 16yz + 12z²,
    q₁₁ = 3x + 8z,  q₁₂ = −8y + 6z,  q₁₃ = 16x + 6y + 32z,  q₂₂ = −4x − 6y,  q₂₃ = 6x + 16z,  q₃₃ = 16x + 8y + 12z.

Hence for every point of the box

    |N(ξ − κ)| ≤ B := |N(v)| + Σᵢ |gᵢ(v)| hᵢ + Σ_{i≤j} |q_{ij}(v)| hᵢhⱼ + (h₁³ + 8h₁²h₃ + 4h₁h₂² + 6h₁h₂h₃ + 16h₁h₃² + 2h₂³ + 8h₂h₃² + 4h₃³).

*The table.* The 24 boxes below lie in D; the ranges are in units of 1/16. For each, B is computed exactly
from the rational centre and half-widths, and B < 1 [02: 2.2].

| # | x | y | z | κ | B |
|---|---|---|---|---|---|
| 1 | [−8, −4] | [−8, −4] | [0, 8] | (0, −1, 0) | 31/32 |
| 2 | [−8, −4] | [−4, 0] | [0, 4] | (0, 0, 0) | 15/32 |
| 3 | [−8, −4] | [−4, 0] | [4, 8] | (0, 0, 0) | 247/256 |
| 4 | [−8, −6] | [0, 8] | [0, 8] | (0, 0, 0) | 243/256 |
| 5 | [−6, −4] | [0, 8] | [0, 8] | (0, 0, 0) | 1597/2048 |
| 6 | [−4, 0] | [−8, 8] | [0, 4] | (0, 0, 0) | 229/256 |
| 7 | [−4, −2] | [−8, 0] | [4, 8] | (0, 0, 0) | 2023/2048 |
| 8 | [−4, −2] | [0, 8] | [4, 8] | (0, 0, 0) | 479/512 |
| 9 | [−2, 0] | [−8, 0] | [4, 8] | (0, 0, 0) | 1921/2048 |
| 10 | [−2, 0] | [0, 8] | [4, 8] | (1, 0, 0) | 2017/2048 |
| 11 | [0, 4] | [−8, 0] | [0, 4] | (0, 0, 0) | 37/64 |
| 12 | [0, 2] | [−8, 0] | [4, 6] | (0, 0, 0) | 325/512 |
| 13 | [0, 2] | [−8, 0] | [6, 8] | (−2, 0, 1) | 59/64 |
| 14 | [2, 4] | [−8, −4] | [4, 8] | (1, −1, 0) | 1613/2048 |
| 15 | [2, 4] | [−4, 0] | [4, 8] | (1, −1, 0) | 1933/2048 |
| 16 | [4, 6] | [−8, 0] | [0, 4] | (0, 0, 0) | 443/512 |
| 17 | [6, 8] | [−8, 0] | [0, 4] | (1, −1, 0) | 221/256 |
| 18 | [4, 8] | [−8, 0] | [4, 8] | (1, −1, 0) | 245/256 |
| 19 | [0, 2] | [0, 8] | [0, 4] | (0, 0, 0) | 361/512 |
| 20 | [2, 4] | [0, 8] | [0, 4] | (1, 0, 0) | 503/512 |
| 21 | [4, 8] | [0, 8] | [0, 4] | (1, 0, 0) | 191/256 |
| 22 | [0, 8] | [0, 8] | [4, 6] | (1, 0, 0) | 365/512 |
| 23 | [0, 8] | [0, 8] | [6, 7] | (1, 0, 0) | 3813/4096 |
| 24 | [0, 8] | [0, 8] | [7, 8] | (2, 0, 0) | 805/1024 |

*The boxes cover D.* All box ends are multiples of 1/16, and each of the 16·16·8 = 2048 cells of side 1/16 that
tile D lies inside one of the boxes [02: 2.3]. (Their volumes add up to exactly 2048 cells, so the boxes have
disjoint interiors.) Every point of D lies in a closed cell and hence in a box, where the listed κ gives
|N(ξ − κ)| ≤ B < 1. The same argument works for real ξ; only rational ξ are needed. ∎

In words: the nearest rounding (κ = 0) works on 12 of the boxes; on the other 12 the translation is (1,0,0),
(2,0,0), (0,−1,0), (1,−1,0) or (−2,0,1). The table was found by the search in
[cover_search.py](scripts/cover_search.py); the proof uses only the table.

**(2.2) Division with remainder.** For α, β ∈ 𝓡 with β ≠ 0 there are κ, ρ ∈ 𝓡 with α = κβ + ρ and
|N(ρ)| < |N(β)|.

*Proof.* N(β) ≠ 0 by (1.4). The rational triple ξ = α·β*/N(β) satisfies ξβ = α. Take κ from (2.1) and put
ρ = α − κβ = (ξ − κ)β ∈ 𝓡. Then |N(ρ)| = |N(ξ − κ)|·|N(β)| < |N(β)|. ∎

**(2.3) Remark.** The constant 1 in (2.1) cannot be lowered below ½. For ξ = θ²/2 and any κ, the triple
2(ξ − κ) has even first two coordinates and odd third coordinate, and N(2a, 2b, c) ≡ 4c³ ≡ 4 (mod 8) for odd c
[02: 2.4]; so |N(ξ − κ)| ≥ 4/8 = ½. The Taylor bounds B of the table are not sharp. On each box, write
N(ξ − κ) in the variables s, t, w ∈ [0, 1] that run over the box, and expand it in the 64 products
B_I(s)B_J(t)B_K(w) of the cubic Bernstein polynomials B_i(s) = C(3, i)sⁱ(1 − s)³⁻ⁱ. These products are ≥ 0 and
sum to 1, so |N(ξ − κ)| is at most the largest absolute value of the 64 coefficients, and the eight corner
coefficients are the values of N(ξ − κ) at the corners. Computed exactly, the largest coefficient over all 24
boxes is 395/512 ≈ 0.77, and it is the value at the corner (3/8, 0, 1/4) of box 16, where κ = 0 [02: 2.5]. So
with the table |N(ξ − κ)| ≤ 395/512 on all of D, with equality there. So the division of (2.2), with κ chosen by the table as in the proof of (2.1),
always gives |N(ρ)| ≤ (395/512)|N(β)|.

## 3. Common divisors, primes and square allocation

**(3.1) Gcd and Bézout.** Let α, β ∈ 𝓡, not both 0. Then there is δ ∈ 𝓡 with δ | α, δ | β and δ = λα + μβ
for some λ, μ ∈ 𝓡. Every common divisor of α and β divides δ.

*Proof.* If β = 0, take δ = α. Otherwise put ρ₀ = β and divide: α = κ₁ρ₀ + ρ₁ with |N(ρ₁)| < |N(ρ₀)| (2.2),
then ρ₀ = κ₂ρ₁ + ρ₂, and so on. The values |N(ρⱼ)| are nonnegative integers that strictly decrease, so some
remainder is 0 (by (1.4), N(ρ) = 0 only for ρ = 0). The last nonzero term δ of the sequence ρ₀, ρ₁, ρ₂, …
divides all the earlier ones and hence α and β; back-substitution writes it as λα + μβ. A common divisor of α
and β divides λα + μβ = δ. ∎

We call δ a *gcd* of α and β. We call α and β *coprime* if every common divisor is a unit. This holds exactly
when λα + μβ = 1 for some λ, μ ∈ 𝓡: if they are coprime, the gcd δ is a unit and we multiply by δ⁻¹;
conversely a common divisor of α and β divides 1.

**(3.2) Euclid's lemma.** If γ | αβ and γ, α are coprime, then γ | β.

*Proof.* From λγ + μα = 1 we get β = λγβ + μαβ, and γ divides both terms. ∎

**(3.3) Prime elements.** Let ϖ ∈ 𝓡 with N(ϖ) = ±q for a prime number q. Then ϖ is not a unit, every divisor of
ϖ is a unit or an associate of ϖ, and if ϖ | αβ then ϖ | α or ϖ | β. This applies to θ (q = 2) and to
π = 4θ − 3 (q = 37).

*Proof.* N(ϖ) ≠ ±1, so ϖ is not a unit (1.3a). If ϖ = δγ, then N(δ)N(γ) = ±q, so |N(δ)| = 1 (δ is a unit) or
|N(γ)| = 1 (γ is a unit and δ is an associate of ϖ). Suppose ϖ | αβ and ϖ ∤ α. A gcd of ϖ and α divides ϖ; it
is not an associate of ϖ, since ϖ ∤ α; so it is a unit, ϖ and α are coprime, and ϖ | β by (3.2). ∎

**(3.4) Exact powers.** Let ϖ be as in (3.3) and α ∈ 𝓡 nonzero. Since |N(ϖ)|ⁿ divides |N(α)| ≥ 1 whenever
ϖⁿ | α, there is a largest n ≥ 0 with ϖⁿ | α, and α = ϖⁿα′ with ϖ ∤ α′.

**(3.5) Lemma: coprime square allocation.** Let α, β, W ∈ 𝓡 with W ≠ 0, and u a unit, with

    αβ = uW²    and    α, β coprime.

Then α = u′R² for a unit u′ and some R ∈ 𝓡.

*Proof.* Induction on the positive integer |N(W)|. Let g be a gcd of α and W (3.1).

*Case 1: g is a unit.* Then λα + μW = 1 for some λ, μ, and squaring gives
λ²α² + 2λμαW + μ²W² = 1. Since W² = u⁻¹αβ, α divides every term on the left, so α | 1 and α is a unit:
α = α·1².

*Case 2: g is not a unit.* Write α = gα₁ and W = gW₁. Then gα₁β = ug²W₁², and cancelling g (1.4) gives
α₁β = ugW₁². A common divisor of g and β divides α and β, so g and β are coprime, and g | α₁β gives g | α₁
(3.2). Write α₁ = gα₂. Then gα₂β = ugW₁², so α₂β = uW₁². Divisors of α₂ divide α, so α₂ and β are coprime. As g
is not a unit, |N(g)| ≥ 2 and |N(W₁)| = |N(W)|/|N(g)| < |N(W)|. By induction α₂ = u′R₁², and α = g²α₂ = u′(gR₁)².

When |N(W)| = 1, W is a unit, so any gcd of α and W divides W and is a unit: Case 1 applies. ∎

**(3.6) What this replaces.** Mordell takes the class number 1 from the table and concludes that unique
factorization holds (p. 1348). The proof below uses only (3.1)–(3.5), which follow from the Euclidean division
(2.2).

**(3.7) Remark: Mordell's fact (1).** By (1.3b) and (1.4) every nonzero rational triple has an inverse, so
𝓡_Q is the field Q(θ). 𝓡 is integrally closed in it: if a rational triple ξ satisfies
ξⁿ + c_{n−1}ξ^{n−1} + ⋯ + c₀ = 0 with integers cᵢ, then ξ ∈ 𝓡. *Proof.* Write ξ = α/m with α ∈ 𝓡 and m a
positive integer, take a gcd δ of α and m, and write α = δα′, m = δm′; then ξ = α′/m′ (cancellation in 𝓡_Q),
and α′, m′ are coprime: a common divisor γ of α′ and m′ gives the common divisor δγ of α and m, so δγ | δ and
γ is a unit. Multiplying the equation by m′ⁿ shows m′ | α′ⁿ. Raising λα′ + μm′ = 1 to the n-th power,
every term except λⁿα′ⁿ is a multiple of m′, and so is λⁿα′ⁿ; hence m′ | 1, m′ is a unit and ξ = α′m′⁻¹ ∈ 𝓡.
∎ Conversely, every integer triple α = (a, b, c) is a root of the monic cubic
t³ − (3a + 8c)t² + (3a² + 16ac − 4b² + 6bc + 16c²)t − N(α), which has integer coefficients [01: 1.8]. So the
elements of Q(θ) that are roots of monic integer polynomials are exactly the integer triples, which is
Mordell's fact (1). The proof of the theorem does not use this remark.

## 4. The units

**(4.1) Theorem.** The units of 𝓡 are exactly the elements ±ε^r η^s (r, s ∈ Z), with ε = θ − 1 and η = 2θ − 1,
and r, s are determined by the unit. Consequently every unit u can be written

    u = σ ε^l η^m v²,    σ ∈ {1, −1},  l, m ∈ {0, 1},  v a unit.

The proof is a descent on the trace height T (1.6). For a unit u, |u₁u₂u₃| = |N(u)| = 1.

**(4.2) Three multipliers.** Put

    μ₁ = ε = (−1, 1, 0),    μ₂ = εη⁻¹ = (−7, 1, 2),    μ₃ = ε⁻²η = (−1, 2, 1)

(products computed exactly [03: 4.2]; all three are units). Their real values, enclosed by exact interval
arithmetic on the root intervals of (1.6) [03: 4.2], are approximately

| k | (μ_k)₁ | (μ_k)₂ | (μ_k)₃ | c_k = max_{i≠k} (μ_k)ᵢ² ≤ | d_k = (μ_k)_k² ≤ | (d_k − c_k)/(1 − c_k) ≤ |
|---|---|---|---|---|---|---|
| 1 | −3.21432 | −0.46081 | 0.67513 | 0.455802 | 10.331852 | 18.148 |
| 2 | 0.59210 | −5.87936 | 0.28726 | 0.350588 | 34.566896 | 52.689 |
| 3 | −0.52543 | 0.36910 | 5.15633 | 0.276075 | 26.587690 | 36.346 |

So μ_k is small at the two places other than k: c_k < 1. The script checks with exact rational bounds that
c_k < 1, d_k > c_k and (1 − c_k)·53 > d_k − c_k for k = 1, 2, 3 [03: 4.2].

**(4.3) Descent lemma.** Let u be a unit with T(u) ≥ 53, and let k be an index with |u_k| ≤ |uᵢ| for all i.
Then T(u·μ_k) < T(u).

*Proof.* |u_k|³ ≤ |u₁u₂u₃| = 1, so u_k² ≤ 1. Using (uμ_k)ᵢ = uᵢ(μ_k)ᵢ,

    T(uμ_k) = Σᵢ uᵢ²(μ_k)ᵢ² ≤ c_k(T(u) − u_k²) + d_k u_k² = c_k T(u) + (d_k − c_k)u_k² ≤ c_k T(u) + (d_k − c_k),

and c_k T(u) + (d_k − c_k) < T(u) because (1 − c_k)T(u) ≥ 53(1 − c_k) > d_k − c_k. ∎

**(4.4) The units of small height.** By (1.6), T(a,b,c) ≤ 52 forces c² ≤ 52·6/37 < 9, then
8(b − 3c/4)² ≤ 52 and 3(a + 8c/3)² ≤ 52, so |c| ≤ 2, |b| ≤ 4, |a| ≤ 9 [03: 4.3]. Among the triples in this box,
exactly 20 have T ≤ 52 and N = ±1, and each is ±ε^r η^s with |r| ≤ 3, |s| ≤ 1, checked by exact multiplication
[03: 4.4]:

| T | units |
|---|---|
| 3 | ±1 |
| 7 | ±ε⁻¹ |
| 11 | ±ε, ±ε²η⁻¹ |
| 15 | ±ε⁻¹η |
| 27 | ±ε⁻², ±ε⁻²η |
| 35 | ±η, ±εη⁻¹ |
| 39 | ±ε³η⁻¹ |

**(4.5) Proof of (4.1).** Let u be a unit. While T(u) ≥ 53, replace u by uμ_k with k as in (4.3); T is a
positive integer and strictly decreases, so after finitely many steps we reach a unit u′ = u·μ_{k₁}⋯μ_{k_j}
with T(u′) ≤ 52. By (4.4), u′ = ±ε^r η^s; each μ_k is a product of powers of ε and η; so
u = u′(μ_{k₁}⋯μ_{k_j})⁻¹ is also of the form ±ε^r η^s.

Uniqueness: from the root intervals, |ε₂|, |η₂|, |ε₃| < 1 < |η₃| [03: 4.5]. It suffices to show that ε^r η^s = ±1
only for r = s = 0. If ε^r η^s = ±1, then |ε₂|^r |η₂|^s = 1 = |ε₃|^r |η₃|^s. The first equation rules out
r, s ≥ 0 and r, s ≤ 0 unless r = s = 0. The second rules out r > 0 > s (both factors < 1) and r < 0 < s (both
factors > 1). So r = s = 0.

Finally write r = 2r₁ + l and s = 2s₁ + m with l, m ∈ {0, 1}: then ±ε^r η^s = ±ε^l η^m (ε^{r₁}η^{s₁})². ∎

**(4.6) Remark.** The eight classes ±ε^l η^m (l, m ∈ {0, 1}) are distinct modulo squares of units, so none of
Mordell's four cases can be reduced to another by absorbing the square of a unit into S². The signs of the real values are (−,−,−) for −1, (−,−,+) for ε and (−,+,+)
for η; these three sign patterns are independent modulo 2, while a square has all values positive [03: 4.5].
This remark is not used.

## 5. The factorization of X − θ

**(5.1) Proposition.** Let (X, Y) be an integer solution of (E). Then there are σ ∈ {1, −1}, l, m ∈ {0, 1} and
S ∈ 𝓡 with

    X − θ = σ ε^l η^m θ S².                                          (11)

This is Mordell's equation (11). The proof takes (5.2)–(5.7).

**(5.2) Two factors.** By (0.2), X = 2X₁ and Y² = 4X₁³ − 4X₁ + 1 with Y odd. Put

    A = (4X₁ − 1, 0, −X₁),       B = (16X₁² + 4X₁ − 15,  2 − 2X₁²,  −4X₁² − X₁ + 4).

Then [04: 5.3, 5.4]

    θA = X − θ,     θ²B = X² + θX + θ² − 4,     (X − θ)(X² + θX + θ² − 4) = X³ − 4X + 2,     AB = η⁻¹Y².

The last identity is the factorization of 2Y² = X³ − 4X + 2 after dividing by θ³ = 2η. Since Y ≠ 0, A and B
are nonzero.

**(5.3) Common divisors of A and B divide π.** We have the identities [04: 5.5, 5.6]

    θB − (X + 2θ)A = (−8, 3, 2) = θπη⁻¹,        A + 1 = X₁θ²η⁻¹.

Let δ be a common divisor of A and B. The first identity gives δ | θπ. The second gives θ | A + 1, so θ ∤ A
(otherwise θ | 1, but θ is not a unit). A gcd of δ and θ divides θ, so by (3.3) it is a unit or an associate
of θ; an associate of θ would divide δ and so A, which is impossible. So δ and θ are coprime, and δ | π by
(3.2).

**(5.4) Separating π.** By (3.4) write A = π^e A′, B = π^f B′ and Y = π^w Y′ (Y as the triple (Y, 0, 0)) with
π dividing none of A′, B′, Y′. From AB = η⁻¹Y²,

    π^{e+f} A′B′ = η⁻¹ π^{2w} Y′².

If e + f > 2w, cancelling π^{2w} leaves π | η⁻¹Y′², so π | Y′² and π | Y′ (3.3), a contradiction. If
e + f < 2w, then π | A′B′, so π | A′ or π | B′, a contradiction. So e + f = 2w and A′B′ = η⁻¹Y′².

A common divisor of A′ and B′ divides A and B, hence π (5.3), hence is a unit or an associate of π (3.3); an
associate of π would give π | A′. So A′ and B′ are coprime, and (3.5), applied with the unit η⁻¹ and W = Y′,
gives A′ = νR² with ν a unit. Thus

    A = ν π^e R².

**(5.5) e is even.** N(A) = −Y² [04: 5.7]. Taking norms in A = νπ^eR² gives −Y² = N(ν)·37^e·N(R)², so
Y² = ±37^e N(R)². If e were odd, then Y² = ±37·Z² with Z = 37^{(e−1)/2}N(R) an integer. Z is odd because Y is
odd, so Z² ≡ 1 (mod 8), and Y² ≡ ±37 ≡ 5 or 3 (mod 8). But Y² ≡ 1 (mod 8) by (0.2) [04: 5.8]. So e is even, and
A = ν(π^{e/2}R)².

**(5.6) Units modulo squares.** By (4.1), ν = σε^l η^m v² with σ = ±1, l, m ∈ {0, 1} and v a unit. Hence
A = σε^l η^m S² with S = vπ^{e/2}R, and X − θ = θA gives (11). ∎

**(5.7) Comparison with Mordell.** Mordell shows that θ divides X − θ exactly once and that the only possible
common prime of the two factors besides θ is 4θ − 3, and he removes an odd power of 4θ − 3 by the same congruence
modulo 8. Here the exact division by θ and θ² is done by the explicit triples A and B, and unique
factorization is replaced by (3.3)–(3.5).

## 6. The four cases

Write S = (a, b, c) and S² = (A₀, A₁, A₂), where [05: 6.0]

    A₀ = a² − 4bc,    A₁ = 2ab + 8bc − 2c²,    A₂ = b² + 2ac + 4c².

Note A₀ + 4A₂ = (a + 4c)² + 4b² − 4bc ≡ a² (mod 4), and a square is 0 or 1 modulo 4. In each case we compare
the three coordinates of the two sides of (11), or of θ times (11).

**(6.1) Case (l, m) = (0, 0): X − θ = σθS².** Here θS² = (−2A₂, A₀ + 4A₂, A₁) [05: 6.1], so

    X = −2σA₂,     A₀ + 4A₂ = −σ,     A₁ = 0.

Since A₀ + 4A₂ ≡ a² ∈ {0, 1} (mod 4) and −σ = ±1, we get −σ = 1: σ = −1, A₀ + 4A₂ = 1 and X = 2A₂. Put
d = a + 4c. As A₁ = 2(bd − c²), the equations read

    bd = c²,       d² + 4b² − 4bc = 1.

- If b = 0, then c = 0 and d = a = ±1; X = 2A₂ = 0.
- If b ≠ 0, multiply the second equation by b² and use bd = c²: c⁴ − 4b³c + 4b⁴ = b². By the identity
  c⁴ − 4b³c + 3b⁴ = (c − b)²((c + b)² + 2b²) [05: 6.1] this says

      (c − b)²((c + b)² + 2b²) = b²(1 − b²).

  The left side is ≥ 0, so b² ≤ 1, hence b² = 1 and the left side is 0. As (c + b)² + 2b² > 0, c = b; then
  d = c²/b = b and a = d − 4c = −3b. X = 2A₂ = 2(b² + 2ac + 4c²) = 2(1 − 6 + 4) = −2.

So this case gives X ∈ {0, −2}, from (a, b, c) = (±1, 0, 0) and ±(−3, 1, 1) [05: 6.1].

**(6.2) Case (l, m) = (1, 0): X − θ = σεθS².** Here εθ = θ² − θ and
(θ² − θ)S² = (2A₂ − 2A₁, −A₀ + 4A₁ − 6A₂, A₀ − A₁ + 4A₂) [05: 6.2], so

    X = σ(2A₂ − 2A₁),     −A₀ + 4A₁ − 6A₂ = −σ,     A₀ − A₁ + 4A₂ = 0.

A₁ is even and A₀ ≡ a² (mod 2). The third equation gives a² ≡ 0, so a is even; the second gives a² ≡ 1 (mod 2),
so a is odd [05: 6.2]. No solutions.

**(6.3) Case (l, m) = (0, 1): X − θ = σηθS².** Multiply by θ and put S′ = θS: θ(X − θ) = σηS′². From now on
(a, b, c) denote the coordinates of S′, and A₀, A₁, A₂ those of S′² (same formulas). Forgetting that S′ is a
multiple of θ can only enlarge the set of candidate X, and every candidate is checked in Section 7. In fact
nothing is lost: θ·(p, q, r) = (−2r, p + 4r, q), so S′ is a multiple of θ exactly when a is even, which the first
equation below forces. Now θ(X − θ) = (0, X, −1) and ηS′² = (−A₀ − 4A₂, 2A₀ − A₁ + 8A₂, 2A₁ − A₂) [05: 6.3], so

    A₀ + 4A₂ = 0,     X = σ(2A₀ − A₁ + 8A₂),     A₂ − 2A₁ = σ.

The first equation gives a² ≡ 0 (mod 4), so a is even; write a = 2(k − 2c), that is k = a/2 + 2c. Then
[05: 6.3]

    A₀ + 4A₂ = 4(k² − b(c − b)),        A₂ − 2A₁ = 4k(c − 2b) + b².

So b² ≡ σ (mod 4), which forces σ = 1 and b odd, and the equations become

    k² = b(c − b),        4k(c − 2b) = 1 − b².

Squaring the second and substituting the first gives (1 − b²)² = 16b(c − b)(c − 2b)² [05: 6.3]. So b divides
(1 − b²)² = 1 − b(2b − b³), hence b | 1 and b = ±1. Then 1 − b² = 0, so (c − b)(c − 2b)² = 0:
- c = b: then k² = 0, k = 0 and a = −4b;
- c = 2b: then k² = b² = 1, k = ±1 and a = 2(k − 4b).

With b = 1 this gives (a, b, c) = (−4, 1, 1), (−6, 1, 2), (−10, 1, 2), and b = −1 gives their negatives. All six
satisfy the equations, and X = 2A₀ − A₁ + 8A₂ equals 2, 4, 12 respectively [05: 6.3].

Mordell splits this case by c = 2b or c ≠ 2b; the divisibility b | (1 − b²)² handles both at once.

**(6.4) Case (l, m) = (1, 1): X − θ = σεηθS².** Multiply by θ and again put S′ = θS = (a, b, c):
θ(X − θ) = σεηS′² with εη = 2θ² − 3θ + 1 = (1, −3, 2). Here
εηS′² = (A₀ − 4A₁ + 6A₂, −3A₀ + 9A₁ − 16A₂, 2A₀ − 3A₁ + 9A₂) [05: 6.4], so

    A₀ − 4A₁ + 6A₂ = 0,     X = σ(−3A₀ + 9A₁ − 16A₂),     2A₀ − 3A₁ + 9A₂ = −σ.

The first equation is ≡ a² (mod 2), so a is even; with a even it is ≡ 2b² (mod 4), so b is even. With a and
b even, 2A₀ − 3A₁ + 9A₂ is even [05: 6.4], so it cannot be ±1. No solutions.

## 7. Conclusion

By (5.1) every integer solution of (E) falls into one of the four cases, so X ∈ {−2, 0, 2, 4, 12}. Then
2Y² = X³ − 4X + 2 gives Y = ±1, ±1, ±1, ±5, ±29, and x = X/2 − 1, y = (Y − 1)/2 give [05: 6.5]

| X | Y | x | y | y(y+1) = x(x+1)(x+2) |
|---|---|---|---|---|
| −2 | ±1 | −2 | 0, −1 | 0 |
| 0 | ±1 | −1 | 0, −1 | 0 |
| 2 | ±1 | 0 | 0, −1 | 0 |
| 4 | ±5 | 1 | 2, −3 | 6 |
| 12 | ±29 | 5 | 14, −15 | 210 |

All of these are solutions of (1). This proves the theorem. ∎

As a confirmation that the proof does not use: a direct search finds no other even X with |X| ≤ 2·10⁶ for which
(X³ − 4X + 2)/2 is a square [06: 7.1], and each of the five values has an explicit representation (11), with
(l, m) = (0, 0) for X = −2, 0 and (0, 1) for X = 2, 4, 12 [06: 7.2].

## 8. Where Mordell's inputs went

| Mordell's input | Source in his paper | Here |
|---|---|---|
| (1) the integers of Q(θ) are a + bθ + cθ² | Delaunay–Faddeev table | Not needed: all arithmetic is in the ring 𝓡 of integer triples. Proved anyway in (3.7). |
| (2) class number 1, hence unique factorization | Delaunay–Faddeev table | Replaced by the Euclidean division (2.1)–(2.2), and its consequences (3.1)–(3.5). |
| (3) fundamental units θ − 1, 2θ − 1 | Delaunay–Faddeev table | Theorem (4.1), by descent on T(u) = Tr(u²). |
| 2 = θ³ × unit; θ prime | derived by Mordell (pp. 1348–1349) | The identity θ³ = 2η (1.5) and (3.3) with N(θ) = −2. |
| 4θ − 3 prime, of norm −37 | derived by Mordell (p. 1349) | (3.3) with N(4θ − 3) = +37 (1.5). |
