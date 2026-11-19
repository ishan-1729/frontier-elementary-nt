# x² + x + 1 = y³ by rational ghosts and 2-adic lifting

**Theorem.** The integer solutions of

    x² + x + 1 = y³                                                    (1)

are exactly (x, y) = (0, 1), (−1, 1), (18, 7), (−19, 7).

With (X, Y) = (y, x), equation (1) is the curve

    E :  Y² + Y + 1 = X³.                                              (E)

That is, Y² + Y = X³ − 1, the curve 243.a1 of LMFDB (written there as y² + y = x³ − 1). Nagell (1921) reduced
x² + x + 1 = yⁿ to n = 3, and that case to the Thue equation a³ − 3ab² + b³ = 1. The theorem itself is Ljunggren's
(Acta Math. 75, 1942). He solved the Thue equation with units of a sextic field and Skolem's method, at the prime 2.
Later proofs use quartic Thue equations, Baker's method or elliptic logarithms. The sources, the closest uses of
this proof's mechanism and the novelty verdict are in [literature_check.md](literature_check.md).

The proof below follows the k = 7 proof of x² + 7 = y³ ([theorem card](../mordell_k7/README.md),
[papers/mordell_k7.tex](../../../shared/prior/elementary_diophantine_paper/papers/mordell_k7.tex)) step by step:
1. **Ghosts.** The ghosts are the rational points of E, written in reduced coordinates (A, B, C). The marker is
   B = 1 and the height is the numerator A (Section 0).
2. **Parents.** A 3-isogeny descent through the auxiliary curve E′ : y² + y = x³ + 20. The parents exist by
   Eisenstein pair arithmetic, with ω² = −1 − ω (Section 1), and by integer cube allocation. The second ring is
   the integers, because E′ factors over Z as (y − 4)(y + 5) = x³ (Section 2).
3. **One index.** Chords through G = (1, 0) and the intertwining identity TS = S³T, checked as an identity of
   rational functions, put every rational point at exactly one index n (Sections 3 and 4).
4. **Marker.** The reduced denominators are |h_n|, with h_{n+2}h_{n−2} = h_{n+1}h_{n−1} + 6h_n² (Section 5). The
   multiples of 3 have 3 | B (4.8). A factored stride identity gives two Somos-4 sequences (Section 6).
5. **Exclusion.** Exact lifting at the prime 2 excludes |h_n| = 1 for every n ≢ 0 (mod 3) except n = ±1, ±2
   (Sections 7 and 8).

**What the proof uses:** integer and rational arithmetic, polynomial identities, congruences, induction and
descent, and elementary inequalities between real numbers. It uses Euclid's algorithm and unique factorization in
Z, and Euclidean division in the Eisenstein pairs, which is proved in Section 1. It uses no reciprocity law, no
ideal or class-number theory, no unit theorem, no group law on the cubic, no p-adic numbers or completions, and no
linear forms in logarithms. Every map identity is a finite polynomial identity, checked in exact arithmetic and
applied at points through Lemma (0.7).

**Checks.** Every identity and finite enumeration below is checked by a script in [scripts/](scripts/), whose
output is stored next to it. Inequalities are derived in the text; the scripts confirm them on examples. A
reference such as [02: 2.3] means section 2.3 of `02_maps.py`, with the result in `02_maps.out.txt`. The section
numbers of the scripts 00–07 are those of the proof; the sections of 08 are named after what they confirm. A line marked NOTE in an output is a step of reasoning that the
text proves and the script does not compute. A line marked "confirmation" checks examples that the proof does not
use.

| Script | Checks |
|---|---|
| [00_reduction.py](scripts/00_reduction.py) ([output](scripts/00_reduction.out.txt)) | reduced coordinates, 3 ∤ A, the three classes, the reflection, the arithmetic of the function field (Section 0) |
| [01_pairs.py](scripts/01_pairs.py) ([output](scripts/01_pairs.out.txt)) | Eisenstein pairs: ring laws, norm, division, Bézout, units, the cube rule, cube formulas, residues (Section 1) |
| [02_maps.py](scripts/02_maps.py) ([output](scripts/02_maps.out.txt)) | the curve E′, the maps Φ and Ψ, their reduced forms, the parent identities (Section 2) |
| [03_chords.py](scripts/03_chords.py) ([output](scripts/03_chords.out.txt)) | chord maps, the reduced chord image, the class shift, the height bound, tangents (Section 3) |
| [04_index.py](scripts/04_index.py) ([output](scripts/04_index.out.txt)) | S³ = the chord through Q₃, TS = S³T, TN = NT, exceptional inputs, the orbit (Section 4) |
| [05_marker.py](scripts/05_marker.py) ([output](scripts/05_marker.out.txt)) | the gcd product, the signed sequence, the chord product, the Somos-4 recurrence of h_n (Section 5) |
| [06_stride.py](scripts/06_stride.py) ([output](scripts/06_stride.out.txt)) | the products 𝒫₃ and 𝒫₆, the stride recurrence, the two sequences s_k (Section 6) |
| [07_lifting.py](scripts/07_lifting.py) ([output](scripts/07_lifting.out.txt)) | the scaled map modulo 4, L⁸ ≡ I, the table modulo 64 (period 8 modulo 16), the drifts (Section 7) |
| [08_crosscheck.py](scripts/08_crosscheck.py) ([output](scripts/08_crosscheck.out.txt)) | confirmations not used by the proof: a direct search, the descent run on 88 points, the prime 3, division-polynomial values |

The shared definitions are in [curve.py](scripts/curve.py). Run each script from the `scripts/` folder with Python 3
and sympy (checked with Python 3.13.5 and sympy 1.14.0; 08 also uses mpmath).

---

## 0. Reduced coordinates and three classes

**(0.1)** Integers x, y satisfy (1) exactly when (X, Y) = (y, x) is an integer point of E. A rational point of E is
a pair (X, Y) of rationals satisfying (E). Write O for an extra symbol, the point at infinity, and Ê for the set of
rational points together with O.

**(0.2) Reduced coordinates.** Every rational point (X, Y) of E can be written uniquely as

    X = A/B²,  Y = C/B³,  with B > 0, gcd(A, B) = 1, and C² + CB³ + B⁶ = A³,        (0.2)

and every such triple (A, B, C) gives a point.

*Proof.* Write Y = c/e and X = a/d in lowest terms with d, e > 0. Then X³ = (c² + ce + e²)/e² [00: 0.2]. Since
c² + ce + e² = c² + e(c + e), its gcd with e equals gcd(c², e) = 1, so this fraction is in lowest terms. Comparing
it with the reduced fraction a³/d³ gives d³ = e² and a³ = c² + ce + e². For each prime p, 3v_p(d) = 2v_p(e), so
v_p(d) is even and v_p(e) is a multiple of 3. Thus d = B², e = B³ with B > 0, and (A, C) = (a, c).
- *Uniqueness.* gcd(A, B) = 1 gives gcd(A, B²) = 1, so A/B² is the lowest-terms form of X; this fixes A and
  B > 0, and then C = YB³.
- *Converse.* (C/B³)² + C/B³ + 1 = (C² + CB³ + B⁶)/B⁶ = A³/B⁶ = (A/B²)³ [00: 0.2]. ∎

The integer points are exactly the triples with B = 1. In (0.2), gcd(B, C) = 1, since a prime dividing B and C
divides A³.

**(0.3) Bounds.** 4A³ = (2C + B³)² + 3B⁶ [00: 0.3]. Hence:
- A ≥ 1, since the right side is positive;
- 3B⁶ ≤ 4A³, so B⁶ ≤ (4/3)A³ < (2A)³ and B² < 2A;
- (2C + B³)² ≤ 4A³, so |2C + B³| ≤ 2A^{3/2}.

**(0.4) Lemma.** 3 ∤ A and A ≡ 1 (mod 3).

*Proof.* If 3 | A then 3 ∤ B, because gcd(A, B) = 1. For 3 ∤ B, the value C² + CB³ + B⁶ modulo 9 lies in
{1, 3, 4, 7}, never 0 [00: 0.4]. This is an enumeration over all residues of B and C modulo 9, on which the value
modulo 9 depends. But 9 | A³. So 3 ∤ A. Modulo 3, C² + CB³ + B⁶ takes only the values 0 and 1 [00: 0.4], and
A³ ≡ A (mod 3). Hence A ≡ 1. ∎

**(0.5) Definition and lemma (classes).** A reduced point is in
- **class 0** if 3 | B;
- **class 1** if 3 ∤ B and 9 | C;
- **class 2** if 3 ∤ B and 9 | C + B³.

Every reduced point lies in exactly one class.

*Proof.* Suppose 3 ∤ B. Since A ≡ 1 (mod 3), A³ ≡ 1 (mod 9), and B⁶ ≡ 1 (mod 9) [00: 0.5]. So
C(C + B³) ≡ 0 (mod 9). The factors C and C + B³ differ by B³, which is a unit modulo 3, so they are not both
divisible by 3. Hence exactly one of them is divisible by 9 [00: 0.5]. ∎

**(0.6) Reflection.** N(X, Y) = (X, −1 − Y) maps E to E, since (−1 − Y)² + (−1 − Y) = Y² + Y. On triples
N(A, B, C) = (A, B, −C − B³) [00: 0.6]. N preserves A, B and class 0, and it exchanges classes 1 and 2 [00: 0.6].
Put N(O) = O. Four integer points are

    G = (1, 1, 0) (class 1),  −G := N(G) = (1, 1, −1) (class 2),
    2G := (7, 1, −19) (class 2),  −2G := (7, 1, 18) (class 1)          [00: 0.6].

The names are labels for now; Section 4 shows that these are the points at indices ±1, ±2.

**(0.7) Lemma (evaluating identities on E).** Let f_E = Y² + Y + 1 − X³ and Q(E) = Q(X)[Y]/(f_E).
- *The field Q(E).* f_E has no root in Q(X), because its discriminant 4X³ − 3 has odd degree and so is not a
  square. Every element of Q(E) has exactly one form a(X) + b(X)Y with a, b ∈ Q(X), and products reduce Y² to
  X³ − 1 − Y. Moreover (a + bY)((a − b) − bY) = a² − ab − b²(X³ − 1), and this norm is a² if b = 0 and
  b²·f_E(−a/b) if b ≠ 0 [00: 0.7]. So it is nonzero for (a, b) ≠ (0, 0), every nonzero element has the inverse
  ((a − b) − bY)/(a² − ab − b²(X³ − 1)), and Q(E) is a field.
- *The lemma.* Let e₁ and e₂ be rational expressions in X and Y built from +, −, × and ÷, for example two
  composites of the maps below. Suppose their images in Q(E) are equal, which is what the scripts check. Let Q be
  a rational point of E at which every division met in evaluating e₁ and e₂ step by step is by a nonzero number.
  Then e₁(Q) = e₂(Q).

*Proof.* f_E is monic in Y [00: 0.7]. So if f_E divides a polynomial g ∈ Q[X, Y] in Q(X)[Y], then division by f_E
in Q[X][Y] has the same remainder, 0, and g = f_E·g₁ with g₁ ∈ Q[X, Y]; in particular g vanishes at every point of E.
1. *Each expression is a fraction that can be evaluated at Q.* By induction on the expression, e = p/q with
   p, q ∈ Q[X, Y], q(Q) ≠ 0 and e(Q) = p(Q)/q(Q). For a quotient e′/e″ with e′ = p′/q′ and e″ = p″/q″, take
   p = p′q″ and q = q′p″. Here p″(Q) ≠ 0 because e″(Q) ≠ 0. Since q(Q) ≠ 0, q is not a multiple of f_E, so q is
   invertible in Q(E) and e = p/q holds there.
2. *The difference is a multiple of f_E.* Π = p₁q₂ − p₂q₁ maps to 0 in Q(E), so Π = f_E·π with π ∈ Q[X, Y].
3. *Evaluation.* At Q, f_E(Q) = 0, so Π(Q) = 0. Since q₁(Q)q₂(Q) ≠ 0, e₁(Q) = e₂(Q). ∎

The lemma is about the expressions, not about their normal forms a + bY: the coefficients a and b can have poles at
such a point. For example, the coefficients of S₊(S₊(X, Y)) (Section 3) have the denominators (X − 7)² in the
X-coordinate and (X − 7)³ in the Y-coordinate, while S₊(S₊(Q₂)) = Q₄ is an ordinary rational point [00: 0.7].
Every identity of rational functions below is applied at points through this lemma.

## 1. Eisenstein pairs

**(1.1) Pairs.** An integer pair (a, b) stands for a + bω, where ω² = −1 − ω. The product is

    (a, b)(c, d) = (ac − bd, ad + bc − bd).

This is the product of a + bω and c + dω reduced modulo ω² + ω + 1. With coordinatewise addition the pairs form a
commutative ring with unit (1, 0), and ω³ = 1 [01: 1.1]. Rational pairs are multiplied by the same rule.

**(1.2) Conjugate and norm.** Put conj(a, b) = (a − b, −b), which sends ω to ω² = −1 − ω, and
Nm(a, b) = a² − ab + b². Then [01: 1.2]:
- μ·conj(μ) = (Nm(μ), 0);
- conj(μν) = conj(μ)conj(ν) and Nm(μν) = Nm(μ)Nm(ν);
- 4Nm(a, b) = (2a − b)² + 3b², so Nm(μ) ≥ 0, with equality only at μ = 0.

Hence a nonzero pair cancels: if μν = μρ with μ ≠ 0, multiplying by conj(μ) gives Nm(μ)ν = Nm(μ)ρ, and an integer
acts coordinatewise [01: 1.2]. Write μ | ν when ν = μκ for an integer pair κ; then Nm(μ) divides Nm(ν). A **unit** is
a pair that divides (1, 0). A pair is a unit exactly when its norm is 1: if μκ = (1, 0) then Nm(μ)Nm(κ) = 1, and if
Nm(μ) = 1 then μ·conj(μ) = (1, 0).

**(1.3) Lemma (division with remainder).** For integer pairs μ, ν with ν ≠ 0 there are integer pairs κ, ρ with
μ = κν + ρ and Nm(ρ) ≤ (7/16)Nm(ν).

*Proof.* Let (x, y) be the rational pair μ·conj(ν)/Nm(ν), so that μ = (x, y)ν. Choose an integer n₂ with
y = n₂ + δ and |δ| ≤ 1/2. Then choose an integer n₁ with x = n₁ + t and |t − δ/2| ≤ 1/2. Put κ = (n₁, n₂) and
ρ = μ − κν = (t, δ)ν. Then Nm(ρ) = Nm(t, δ)Nm(ν), and

    Nm(t, δ) = (t − δ/2)² + (3/4)δ² ≤ 1/4 + 3/16 = 7/16                [01: 1.3].  ∎

[01: 1.3] confirms the bound on 20000 random divisions.

**(1.4) Gcd and Bézout.** Let μ, ν be integer pairs, not both 0. Put r₀ = μ, r₁ = ν, and divide repeatedly,
r_{i−1} = κ_i r_i + r_{i+1}. The norms of the remainders strictly decrease, so some r_{m+1} = 0. Put g = r_m.
- *Forward induction.* Each r_i is a combination ξ_iμ + η_iν with integer pairs ξ_i, η_i. This holds for i = 0, 1,
  and r_{i+1} = r_{i−1} − κ_i r_i. So g = ξμ + ην.
- *Backward induction.* g divides r_m and r_{m+1} = 0, and r_{i−1} = κ_i r_i + r_{i+1}. So g divides every r_i, in
  particular μ and ν.
- *Every common divisor divides g,* since it divides ξμ + ην.

We call μ and ν **coprime** when every common divisor is a unit. Then g is a unit, and multiplying by its inverse
gives ξ′μ + η′ν = (1, 0). [01: 1.4] confirms the construction on 3000 random pairs.

**(1.5) Units.** Nm(a, b) = 1 means (2a − b)² + 3b² = 4, so |b| ≤ 1. The six solutions are ±(1, 0), ±(0, 1) and
±(1, 1), that is ±1, ±ω and ±ω² = ∓(1 + ω) [01: 1.5]. Each unit is ±ω^j with j ∈ {0, 1, 2}, and −1 = (−1)³.

**(1.6) Lemma (cube rule).** If U, V are coprime integer pairs, W is an integer pair and UV = W³ ≠ 0, then
U = ω^j θ³ for an integer pair θ and some j ∈ {0, 1, 2}.

*Proof.* We use induction on Nm(W). Let g be a gcd of U and W; g ≠ 0, since g | U and U ≠ 0.
- *g is a unit.* Write (1, 0) = ξU + ηW. Cubing gives (1, 0) = Uλ + η³W³ with
  λ = ξ³U² + 3ξ²UηW + 3ξη²W² [01: 1.6]. Since W³ = UV, (1, 0) = U(λ + η³V). So U is a unit, U = ±ω^j = ω^j(±1)³.
- *g is not a unit.* Write U = gU₁ and W = gW₁. Cancelling g in gU₁V = g³W₁³ gives U₁V = g²W₁³.
  - g and V are coprime, since a common divisor divides U and V. So (1, 0) = ξg + ηV. Squaring gives
    (1, 0) = ξ′g² + η′V with ξ′ = ξ² and η′ = 2ξηg + η²V [01: 1.6].
  - Then U₁ = ξ′g²U₁ + η′(VU₁) = g²(ξ′U₁ + η′W₁³). So g² | U₁. Write U₁ = g²U₂.
  - Cancelling g² gives U₂V = W₁³. Here U₂ and V are coprime, since U₂ divides U, and
    Nm(W₁) = Nm(W)/Nm(g) < Nm(W).
  - By induction U₂ = ω^j θ₁³, so U = g³U₂ = ω^j(gθ₁)³. ∎

[01: 1.6] runs this proof literally, with the Bézout coefficients, on 318 random examples and recovers
U = unit·θ³ each time.

**(1.7) Cubes.** For integers u, v [01: 1.7]:

    (u + vω)³ = (u³ − 3uv² + v³) + 3uv(u − v)·ω,
    ω·(x, y) = (−y, x − y),   ω²·(x, y) = (y − x, −x).

In particular the three pairs (u, v), ω(u, v), ω²(u, v) have the same cube and the second coordinates v, u − v
and −u.

**(1.8) Residues modulo 3** [01: 1.8]: u³ − 3uv² + v³ ≡ u + v and u² − uv + v² ≡ (u + v)². If 3 ∤ u + v, then
exactly one of u, v, u − v is divisible by 3.

## 2. The auxiliary curve and the two maps

**(2.1) The auxiliary curve.** Let

    E′ :  y² + y = x³ + 20,   equivalently (y − 4)(y + 5) = x³,          [02: 2.1]

and let N′(x, y) = (x, −1 − y). Every rational point of E′ can be written uniquely as x = D/E², y = F/E³ with
E > 0, gcd(D, E) = 1 and

    F² + FE³ = D³ + 20E⁶,   that is (F − 4E³)(F + 5E³) = D³           [02: 2.1].

*Proof.* As in (0.2), write y = c/e and x = a/d in lowest terms with d, e > 0. Then
x³ = (c² + ce − 20e²)/e², and c² + ce − 20e² = c² + e(c − 20e) [02: 2.1], so this fraction is in lowest terms; its
numerator may be negative or zero.
- If x ≠ 0, comparing with a³/d³ gives d³ = e², hence d = E², e = E³, and (D, F) = (a, c).
- If x = 0, then y² + y = 20, so y = 4 or y = −5 [02: 2.1], and (D, E, F) = (0, 1, 4) or (0, 1, −5).

Uniqueness follows as in (0.2). Finally gcd(F, E) = 1, since a prime dividing F and E divides D³. ∎

Inside a triple (D, E, F), E denotes a coordinate; the curves are always named E and E′ with their equations.
A triple (D, E, F) and the triple (D, −E, −F) give the same point; so do (A, B, C) and (A, −B, −C) on E. Below,
a triple is sometimes written with a negative middle entry and normalized at the end.

*Why the second ring is Z.* In the k = 7 proof the auxiliary curve F² = D³ + 189E⁶ needed pair arithmetic for the
real quadratic norm a² + ab − 5b², of discriminant 21. Here the quadratic y² + y − 20 has discriminant
1 + 80 = 81 = 9² [02: 2.1], so it splits over Z. The "pair ring" of E′ is Z × Z, and its cube extraction is ordinary
cube allocation in Z. The two rational points (0, 4) and (0, −5) are the poles of Ψ below.

**(2.2) The maps.** Define

    Φ(X, Y) = ( (X³ − 3)/X²,  (X³Y + 6Y + 3)/X³ ),
    Ψ(x, y) = ( (x³ + 81)/(9x²),  (x³y − 13x³ − 162y − 81)/(27x³) ).

They satisfy

    X⁶(Φ_y² + Φ_y − Φ_x³ − 20) = −(X³ + 6)²(X³ − Y² − Y − 1),
    729x⁶(Ψ_Y² + Ψ_Y + 1 − Ψ_X³) = −(x³ − 162)²(x³ − y² − y + 20)       [02: 2.2].

So Φ maps rational points of E to E′, and Ψ maps rational points of E′ with x ≠ 0 to E. Φ is defined at every
rational point of E, since X = A/B² > 0 there.

**(2.3) Lemma (Φ on reduced triples).** For a reduced point (a, b, c) of E,

    Φ(a, b, c) = (D, E, F) = (a³ − 3b⁶,  ab,  c(a³ + 6b⁶) + 3b⁹)            [02: 2.3]

is the reduced triple of the image. It has 3 ∤ D and

    F − 4E³ = (c − b³)³,   F + 5E³ = (c + 2b³)³,   D = (c − b³)(c + 2b³)      [02: 2.3].

*Proof.* The first display is (2.2) after substituting X = a/b² and Y = c/b³ [02: 2.3]. Here E = ab > 0.
Suppose a prime p divides ab and a³ − 3b⁶.
- If p | b, then p | a³, contradicting gcd(a, b) = 1.
- If p | a, then p | 3b⁶ with p ∤ b, so p = 3. But 3 ∤ a by (0.4).

So the triple is reduced. Since D − a³ = −3b⁶, 3 ∤ D. The three identities follow from
a³ = c² + cb³ + b⁶ [02: 2.3]. ∎

**(2.4) Lemma (Ψ on reduced triples).** Let (D, E, F) be a reduced point of E′ with 3 ∤ D. The reduced triple of
Ψ(D, E, F) is

    A = D³ + 81E⁶,   B = 3|D|E,   C = sgn(D)·(D³F − 13D³E³ − 162FE⁶ − 81E⁹)     [02: 2.4].

In particular 3 | B, and

    A = F² + FE³ + 61E⁶ = (F + E³/2)² + (243/4)E⁶ ≥ (243/4)E⁶          [02: 2.4].

*Proof.* Substituting x = D/E² and y = F/E³ in (2.2) gives X = (D³ + 81E⁶)/(3DE)² and
Y = (D³F − 13D³E³ − 162FE⁶ − 81E⁹)/(3DE)³ [02: 2.4]. It remains to show gcd(A, 3DE) = 1.
- Modulo 3, A ≡ D³ ≢ 0.
- A prime p | D has A ≡ 81E⁶ ≢ 0 (mod p), because p ∤ 3E.
- A prime p | E has A ≡ D³ ≢ 0 (mod p).

When D < 0, the triple (A, 3DE, C₀) with C₀ = D³F − ⋯ − 81E⁹ has negative middle entry; replacing it by
(A, −3DE, −C₀) changes neither X nor Y, and gives the stated B and C. [02: 2.4] compares the formula with the
reduced triple of the image point on examples with D of both signs. The height identity uses
D³ = F² + FE³ − 20E⁶. ∎

**(2.5) The map T.** For a rational point P of E put T(P) = Ψ(Φ(P)), and put T(O) = O. By (2.3) and (2.4), for
every rational point P = (a, b, c):
- T(P) is a rational point, not O;
- T(P) is in class 0;
- A(T(P)) ≥ (243/4)(ab)⁶ ≥ (243/4)a⁶.

Moreover T(N(P)) = N(T(P)), because Φ∘N = N′∘Φ and Ψ∘N′ = N∘Ψ [02: 2.5], [04: 4.3]. For example
T(G) = (73, 6, 487) and T(−G) = (73, 6, −703) [02: 2.5].

**(2.6) Lemma (Φ-parents).** Every reduced point (D, E, F) of E′ with 3 ∤ D is Φ(a, b, c) for a reduced point
(a, b, c) of E.

*Proof.* KEY STEPS:
1. *The two factors are coprime.* F + 5E³ − (F − 4E³) = 9E³. A common prime p of the two factors divides 9E³ and
   their product D³. If p | E then p | D, contradicting gcd(D, E) = 1. If p = 3 then 3 | D, contradicting the
   hypothesis. So the factors are coprime.
2. *Both are cubes.* Coprime integers whose product is a cube are cubes, signs included, since −1 = (−1)³. Write
   F − 4E³ = σ³ and F + 5E³ = τ³ with σ, τ ∈ Z. Then (στ)³ = D³ gives στ = D, so 3 ∤ στ, and
   τ³ − σ³ = 9E³.
3. *τ ≡ σ (mod 3).* For 3 ∤ z, z³ ≡ ±1 (mod 9) according as z ≡ ±1 (mod 3) [02: 2.6]. So 9 | τ³ − σ³ forces
   τ ≡ σ (mod 3) [02: 2.6].
4. *Split the difference.* τ³ − σ³ = (τ − σ)(τ² + τσ + σ²) and τ² + τσ + σ² = (τ − σ)² + 3τσ [02: 2.6]. As
   9 | (τ − σ)² and 3τσ is divisible by 3 exactly once, the second factor is divisible by 3 exactly once. Hence

       ((τ − σ)/3) · ((τ² + τσ + σ²)/3) = E³,

   a product of two integers.
5. *They are coprime.* The second factor is prime to 3. A prime p ≠ 3 dividing both factors divides
   (τ − σ)² + 3τσ − (τ − σ)² = 3τσ, hence τσ, hence both τ and σ. That contradicts step 1.
6. *Allocate the cubes.* The second factor is positive: 4(τ² + τσ + σ²) = (2τ + σ)² + 3σ² [02: 2.6], and σ ≠ 0
   because στ = D ≠ 0. So (τ² + τσ + σ²)/3 = a³ with a ≥ 1, and (τ − σ)/3 = E³/a³ = b³ with b ≥ 1,
   gcd(a, b) = 1 and ab = E.
7. *The parent.* Put c = σ + b³, so that c − b³ = σ and c + 2b³ = σ + 3b³ = τ. The identity
   c² + cb³ + b⁶ = (τ² + τσ + σ²)/3 [02: 2.6] gives a³ = c² + cb³ + b⁶. Hence (a, b, c) is a reduced point of E.
8. *Φ returns the given point.* By (2.3), Φ(a, b, c) = ((c − b³)(c + 2b³), ab, 4(ab)³ + (c − b³)³)
   = (στ, E, 4E³ + σ³) = (D, E, F). ∎

**(2.7) Lemma (Ψ-parents).** Every reduced point (A, B, C) of E in class 0 equals Ψ(D, E, F) for a reduced point
(D, E, F) of E′ with 3 ∤ D.

*Proof.* KEY STEPS:
1. *Eisenstein factorization.* Put γ = (C, −B³), that is C − B³ω. Then Nm(γ) = C² + CB³ + B⁶ = A³,
   conj(γ) = (C + B³, B³) and γ − conj(γ) = −B³·(1, 2), with Nm(1, 2) = 3 [02: 2.7].
2. *γ and conj(γ) are coprime.* A common divisor δ divides B³(1, 2), so Nm(δ) divides 3B⁶. It also divides
   Nm(γ) = A³. But gcd(A³, 3B⁶) = 1 by (0.2) and (0.4). So Nm(δ) = 1 and δ is a unit.
3. *γ is a cube times a unit.* γ·conj(γ) = (A, 0)³. The cube rule (1.6) gives γ = ω^j θ³ for an integer pair
   θ = (u, v) and some j ∈ {0, 1, 2}. Then A³ = Nm(γ) = Nm(θ)³, so A = Nm(θ) = u² − uv + v².
4. *j = 0.* Write θ³ = (p, 3q) with p = u³ − 3uv² + v³ and q = uv(u − v) (1.7). Then
   ω(p, 3q) = (−3q, p − 3q) and ω²(p, 3q) = (3q − p, −p) [01: 1.7]. By (1.8), 3 ∤ A ≡ (u + v)² gives
   p ≡ u + v ≢ 0 (mod 3). So for j = 1, 2 the second coordinate of γ would be prime to 3. It is −B³, and 3 | B.
   Hence j = 0, and

       C = u³ − 3uv² + v³,   −B³ = 3uv(u − v),   A = u² − uv + v².

5. *Choose the cube root.* θ, ωθ and ω²θ have the same cube. Their second coordinates are v, u − v and −u (1.7),
   and exactly one of u, v, u − v is divisible by 3 (1.8). Replace θ by the root whose second coordinate is
   divisible by 3, and rename it (u, v). The displays of step 4 still hold, now with 3 | v, 3 ∤ u and 3 ∤ u − v.
6. *Coprimality.* A prime dividing u and v divides A and B, so gcd(u, v) = 1. Hence u, v and u − v are pairwise
   coprime. None of them is 0, because B ≠ 0.
7. *Three integer cubes.* B³ = −3uv(u − v). Since u is coprime to 3v(u − v), u is a cube, say u = s³. In the
   same way u − v = t³, with s, t ∈ Z. Then v/9 = (−B/(3st))³ is the cube of a rational number. Also
   v₃(v) = 3v₃(B) − 1 ≥ 2, so v/9 is an integer. A rational number whose cube is an integer is an integer, so
   w := −B/(3st) ∈ Z, with v = 9w³. Then

       s³ − t³ = 9w³,   B = −3stw.

8. *The parent on E′.* Put (D, E, F) = (st, w, 4w³ − s³). Then F − 4E³ = −s³ and F + 5E³ = −t³, so the triple
   lies on E′ [02: 2.7].
   - gcd(D, E) = 1: a prime dividing s and w divides u and v, and a prime dividing t and w divides u − v and v.
   - 3 ∤ D = st.
   - If w < 0, use the same point in the form (st, −w, s³ − 4w³); this changes neither the point nor its
     Ψ-image [02: 2.7].
9. *Ψ returns the given point.* D³ + 81E⁶ = s³t³ + 81w⁶ = u(u − v) + v² = A, 3DE = 3stw = −B, and the third
   Ψ-coordinate equals −(u³ − 3uv² + v³) = −C [02: 2.7]. So Ψ(D, E, F) has X = A/B² and Y = (−C)/(−B)³ = C/B³. ∎

**(2.8) Corollary (T-parents).** A rational point P of E is T(P₀) for a rational point P₀ exactly when P is in
class 0. In that case (243/4)A(P₀)⁶ ≤ A(P), so A(P₀) < A(P).

*Proof.* If P = T(P₀), then P is in class 0 by (2.5). Conversely, let P be in class 0. Lemma (2.7) gives a reduced
point P′ of E′ with 3 ∤ D and Ψ(P′) = P. Lemma (2.6) gives P₀ with Φ(P₀) = P′. So T(P₀) = P. The height bound is
(2.5). ∎

This replaces the k = 7 statements "every auxiliary point has a Φ-parent" and "a point has a Ψ-parent exactly when
A is odd". Here the image of Φ is "3 ∤ D", and the image of Ψ restricted to it is class 0.

## 3. Chords and the descent

**(3.1) Chord maps.** Let R₀ = (t, v) be a rational point of E, and let Q = (X, Y) be a rational point with X ≠ t.
The line through Q and R₀ has slope m = (Y − v)/(X − t). It meets E in a third point with abscissa
X′ = m² − X − t, because

    z³ − (v + m(z − t))² − (v + m(z − t)) − 1 − (z − t)(z − X)(z − X′)

has degree at most 1 in z and vanishes at z = t and z = X [03: 3.1]. Reflecting the third point by N defines

    R_{R₀}(X, Y) = ( m² − X − t,  −1 − v − m(m² − X − 2t) ).

Put S₊ = R_G and S₋ = R_{−G}, the chord maps through G = (1, 0) and −G = (1, −1). In coordinates, S₊ has
m = Y/(X − 1), and S₋ has m = (Y + 1)/(X − 1). As identities of rational functions on E, used through (0.7),

    S₋(S₊(X, Y)) = (X, Y),   S₊(S₋(X, Y)) = (X, Y),   S₋ = N S₊ N           [03: 3.1].

S₊ and S₋ are defined at every rational point except ±G, since X = 1 forces A = B² = 1 (see (3.6)).

**(3.2) Lemma (reduced chord image).** Let P = (A, B, C) be a reduced point other than ±G. Put D₀ = A − B² ≠ 0,
h = gcd(|D₀|, |C|), d = D₀/h, n = C/h, and

    K₊ = A² + AB² − 2B⁴ − BC.

Then S₊(P) has the reduced triple

    A′ = K₊/h²,   B′ = |d|,   C′ = sgn(d)(n(d² − A′)/B − d³).

Likewise S₋(P) has A′ = K₋/h₋² and B′ = |D₀|/h₋, where h₋ = gcd(|D₀|, |C + B³|) and K₋ = A² + AB² − B⁴ + BC.

*Proof.* On E,

    C² − (A + B²)D₀² = B²K₊,   D₀² − K₊ = B(C + 3B³ − 3AB)                [03: 3.2].

- *A′ is an integer.* gcd(h, B) divides gcd(C, B) = 1. Dividing the first identity by h² gives
  n² − (A + B²)d² = B²K₊/h², so h² | K₊ and B²A′ = n² − (A + B²)d².
- *C′ is an integer.* The second identity gives h²(d² − A′) = B(C + 3B³ − 3AB), so B | d² − A′ (again because
  gcd(h, B) = 1).
- *Reduced.* A prime dividing A′ and d would divide n, by the first relation, but gcd(n, d) = 1. So
  gcd(A′, B′) = 1.
- *The image.* Substituting m = C/(B D₀) = n/(Bd) into (3.1) gives X(S₊P) = A′/d² and
  Y(S₊P) = (n(d² − A′)/B − d³)/d³ = C′/B′³ [03: 3.2].
- *S₋.* Use S₋ = N S₊ N: N(P) = (A, B, −C − B³), and K₊ with C replaced by −C − B³ is K₋ [03: 3.2].

[03: 3.2] compares these triples with the chord images of Q_n for 2 ≤ |n| ≤ 40, including cases with d < 0. ∎

**(3.3) Lemma (class shift).**
- If P is in class 1 and P ≠ G, then S₋(P) is in class 0.
- If P is in class 2 and P ≠ −G, then S₊(P) is in class 0.

*Proof.* Let P be in class 1. Then 3 ∤ B and A ≡ 1 (mod 3) give 3 | D₀ = A − B² [03: 3.3]. Also
C + B³ ≡ B³ ≢ 0 (mod 3) [03: 3.3], so 3 ∤ h₋. Hence 3 | |D₀|/h₋ = B(S₋P). Class 2 is the same with C in place of
C + B³: C ≡ −B³ ≢ 0 (mod 3). ∎

**(3.4) Lemma (height of a chord image).** In (3.2), A(S₊P) < (29/2)A² and A(S₋P) < (29/2)A².

*Proof.* By (0.3), B² < 2A and |2C + B³| ≤ 2A^{3/2}. So |B(2C + B³)| < 2√2·A². Since
2BC = B(2C + B³) − B⁴ [03: 3.4], and B⁴ < 4A²,

    |BC| < (√2 + 2)A² < (7/2)A²,

using √2 < 3/2. Then |K₊| ≤ A² + AB² + 2B⁴ + |BC| < (1 + 2 + 8 + 7/2)A² = (29/2)A², and
|K₋| ≤ A² + AB² + B⁴ + |BC| < (21/2)A² [03: 3.4]. Finally A′ = K₊/h² ≤ |K₊| for S₊, and likewise for S₋.
[03: 3.4] confirms the inequalities on the orbit in exact integers. ∎

**(3.5) Tangents and completed values.** At G the line Y = m(X − 1) gives, with X = 1 + z,
m²z² + mz + 1 − (1 + z)³ = z((m − 3) + (m² − 3)z − z²). A double root at z = 0 forces m = 3. The remaining root
z = 6 gives the point (7, 18), and reflection gives (7, −19) [03: 3.5]. Complete the maps by

    S₊(O) = G,  S₊(−G) = O,  S₊(G) = (7, −19) = 2G;
    S₋(O) = −G, S₋(G) = O,   S₋(−G) = (7, 18) = −2G.

Then S₊ and S₋ are mutually inverse bijections of Ê, and N S₊ N = S₋.
- *Generic points.* If Q ∉ {±G} and S₊(Q) ∉ {±G}, then S₋(S₊(Q)) = Q by (3.1) and (0.7). The same holds with S₊
  and S₋ exchanged.
- *When a chord image is ±G.* By (3.1), X(S₊Q) − 1 = (Y² − (X + 2)(X − 1)²)/(X − 1)². Subtracting this numerator
  from the curve equation gives Y = 3(X − 1), and then the numerator is −(X − 7)(X − 1)². For S₋ the same holds
  with Y + 1 = −3(X − 1) [03: 3.5]. So X(S₊Q) = 1 or X(S₋Q) = 1 forces X(Q) ∈ {1, 7}, and for Q ∉ {±G} this leaves
  Q = ±2G. There S₊(−2G) = −G and S₋(2G) = G, while S₊(2G) and S₋(−2G) are not ±G [03: 3.5]. By the completed
  values, S₋(−G) = −2G and S₊(G) = 2G.
- *The exceptional points.* At O, ±G and ±2G all ten compositions S₋S₊ and S₊S₋ are the identity, and
  N S₊ N = S₋ holds at these five points [03: 3.5]. Elsewhere N S₊ N = S₋ is the identity of (3.1).

**(3.6) Seeds.** The rational points with A = 1 are exactly G and −G. Indeed 3B⁶ ≤ 4 gives B = 1, and then
C² + C = 0 [03: 3.6].

**(3.7) Theorem (descent).** Let 𝒯 be the smallest set of rational points of E that contains G and −G and
contains T(P), S₊(T(P)) and S₋(T(P)) whenever it contains P. (Since T(P) is in class 0 by (2.5), T(P) ≠ ±G, so
S₊(T(P)) and S₋(T(P)) are rational points given by (3.1).) Then every rational point of E lies in 𝒯.

*Proof.* Strong induction on A. If A = 1 then P = ±G by (3.6). Let A > 1.
- *Class 0.* By (2.8), P = T(P₀) with A(P₀) < A. By induction P₀ ∈ 𝒯, so P ∈ 𝒯.
- *Class 1.* P ≠ G because A > 1. By (3.3), R = S₋(P) is a rational point of class 0, and A(R) < (29/2)A² by
  (3.4). By (2.8), R = T(P₀) with

      (243/4)A(P₀)⁶ ≤ A(R) < (29/2)A²,  so  A(P₀)⁶ < (58/243)A² < A²      [03: 3.4].

  Thus A(P₀)³ < A, and A(P₀) < A. By induction P₀ ∈ 𝒯. Since P ∉ {±G} and R ∉ {±G} (R is in class 0), (3.5)
  gives P = S₊(S₋(P)) = S₊(T(P₀)) ∈ 𝒯.
- *Class 2.* The same argument with S₊ and S₋ exchanged gives P = S₋(T(P₀)). ∎

The chord may increase the height, but the T-parent that follows decreases it strictly. [03: 3.7] checks the
class shift and the drop in A along Q_n for 2 ≤ |n| ≤ 40. [08: 8.2] constructs the parents of (2.7) and (2.6) for
the 88 points Q_n with 2 ≤ |n| ≤ 45: the Eisenstein cube root is found numerically and verified exactly, the
integer cube roots are exact, and the result is Q_m with n = 3m or n = 3m ± 1.

## 4. One signed index

**(4.1)** Write S = S₊, so S⁻¹ = S₋, and N S N = S⁻¹ on Ê (3.5). Define

    Q_n = Sⁿ(O)   (n ∈ Z).

Then Q₀ = O, Q₁ = G, Q₂ = (7, −19), Q₃ = (73/36, 487/216), Q₋₁ = −G and Q₋₂ = (7, 18). From N S N = S⁻¹,
N(Q_n) = Q₋ₙ. The triples of Q_n for |n| ≤ 6 are listed in [04: 4.1]; for example

    Q₄ = (2443, 37, −137826),  Q₅ = (296221, 179, 158277601),  Q₆ = (48833569, 7140, −312708123503).

**(4.2) Three chords are one chord.** Let R₃ be the chord map through Q₃ = (t, v) = (73/36, 487/216), and R₋₃ the
chord map through Q₋₃ = N(Q₃) = (73/36, −703/216). As identities of rational functions,

    S(S(S(X, Y))) = R₃(X, Y),   S⁻¹(S⁻¹(S⁻¹(X, Y))) = R₋₃(X, Y)           [04: 4.2].

Complete R₃ by R₃(O) = Q₃ and R₃(Q₋₃) = O (the line is vertical), and at Q₃ by the tangent there, reflected. Then
S³ = R₃ on all of Ê:
- *Generic points.* The composite S³ is evaluated by (3.1) when Q, SQ and S²Q are not in {O, G, −G}. By the
  completed values (3.5), this excludes Q ∈ {O, G, −G} ∪ S⁻¹{O, G, −G} ∪ S⁻²{O, G, −G} = {O, G, −G, −2G, Q₋₃}. R₃
  is given by (3.1) when Q ∉ {O, Q₃, Q₋₃}. At every other rational point, S³(Q) = R₃(Q) by (0.7).
- *The six exceptional points* O, G, −G, −2G, Q₋₃, Q₃ are evaluated directly [04: 4.2].

The same argument gives S⁻³(Q) = R₋₃(Q) at every rational point Q outside O, ±G, 2G, Q₃ and Q₋₃, which is all
that (6.4) uses.

**(4.3) The intertwining identity.** As identities of rational functions on E [04: 4.3],

    T(S(X, Y)) = R₃(T(X, Y)),   T(N(X, Y)) = N(T(X, Y)).

The check works in the field Q(E) of (0.7). For each coordinate, the coefficients of 1 and of Y are compared by
cross-multiplying their numerators and denominators. The cross products have degrees 34, 34, 52, 52 in X (35, 32,
53, 50 terms), and every difference is the zero polynomial [04: 4.3].

**(4.4) Lemma.** T S = S³ T on Ê.

*Proof.* Let Q be a rational point. The divisions in the two composites are:
- in S at Q, by X − 1, which vanishes only at Q = ±G;
- in T, at Q and at S(Q), by X and by X³ − 3, which never vanish at a rational point; for Q ∉ {O, ±G}, S(Q) is a
  rational point;
- in R₃ at T(Q), by X(T(Q)) − 73/36.

If X(T(Q)) = 73/36 then T(Q) = Q₃ or Q₋₃, so A(T(Q)) = 73 ≥ (243/4)A(Q)⁶ by (2.5). That forces A(Q) = 1, so
Q = ±G [04: 4.4]. For Q ∉ {O, ±G}, then, (4.3) and (0.7) give T(S(Q)) = R₃(T(Q)), and R₃ = S³ by (4.2). The
remaining inputs are O, G and −G, and there [04: 4.4]:
- T(S(O)) = T(G) = (73/36, 487/216) = S³(O);
- T(S(G)) = T(Q₂) and S³(T(G)) = S³(Q₃) are the same point, (48833569/7140², −312708123503/7140³);
- T(S(−G)) = T(O) = O and S³(T(−G)) = S³(Q₋₃) = O. ∎

**(4.5) Corollary.** T(Q_n) = Q_{3n} for every integer n.

*Proof.* T(Q₀) = O = Q₀. For n ≥ 0, T(Q_{n+1}) = T(S(Q_n)) = S³(T(Q_n)) = S³(Q_{3n}) = Q_{3n+3}. For negative
indices, T(Q₋ₙ) = T(N(Q_n)) = N(T(Q_n)) = N(Q_{3n}) = Q₋₃ₙ. [04: 4.5] confirms T(Q_m) = Q_{3m} for |m| ≤ 6. ∎

**(4.6) Lemma (distinct indices).** Q_i ≠ Q_j for i ≠ j.

*Proof.* If Q_i = Q_j with i < j, applying S⁻ⁱ gives Q_{j−i} = O. Let d ≥ 1 be least with Q_d = O; d ≠ 1, since
Q₁ = G. Then Q_{n+d} = Q_n for all n.
- *3 | d.* T(Q_{d/3}) = Q_d = O, so Q_{d/3} = O by (2.5), against the choice of d.
- *d ≡ 2 (mod 3).* Put m = (d + 1)/3. Then T(Q_m) = Q_{d+1} = Q₁ = G.
- *d ≡ 1 (mod 3).* Put m = (2d + 1)/3. Then T(Q_m) = Q_{2d+1} = Q₁ = G.

Both values are impossible: T(O) = O, and every other T-image is in class 0 by (2.4) and (2.5), while G has
B = 1. ∎

**(4.7) Theorem (one index).** Every rational point of E is Q_n for exactly one integer n ≠ 0.

*Proof.* G = Q₁ and −G = Q₋₁. If P = Q_m with m ≠ 0, then T(P) = Q_{3m}, S₊(T(P)) = Q_{3m+1} and
S₋(T(P)) = Q_{3m−1}, by (4.5); these indices are not 0. So the set {Q_n : n ≠ 0} contains ±G and is closed under
the moves of (3.7), and by the minimality of 𝒯 every rational point is some Q_n with n ≠ 0. Uniqueness is (4.6). ∎

In particular n ↦ Q_n is a bijection from Z onto Ê. This is the elementary content of "Mordell–Weil group Z" for
243.a1. No group law is used: S and T act on indices by n ↦ n + 1 and n ↦ 3n.

**(4.8) Multiples of 3.** For m ≠ 0, Q_m ≠ O by (4.6), so Q_{3m} = T(Q_m) is in class 0 by (2.5), and
3 | B(Q_{3m}). In particular Q_{3m} is never an integer point.

## 5. The denominator marker is exact

Write (A_n, B_n, C_n) for the reduced triple of Q_n (n ≠ 0).

**(5.1) Lemma (gcd product).** For a reduced point P ≠ ±G, put e₊ = gcd(|A − B²|, |C|) and
e₋ = gcd(|A − B²|, |C + B³|). Then e₊e₋ = |A − B²|.

*Proof.* On E, C(C + B³) = (A − B²)(A² + AB² + B⁴), and A² + AB² + B⁴ = (A − B²)² + 3AB² [05: 5.1]. The two
factors C and C + B³ are coprime, because gcd(C, B) = 1.
- *A prime p ≠ 3 with p | A − B².* Such a p does not divide A² + AB² + B⁴: otherwise it divides 3AB², hence A or B,
  and together with p | A − B² this contradicts gcd(A, B) = 1. So the full power of p in A − B² lies in exactly
  one of C and C + B³, and it contributes v_p(A − B²) to e₊e₋.
- *The prime 3, when 3 | A − B².* By (0.4), 3 ∤ A, hence 3 ∤ B, since 3 | A − B². So v₃(3AB²) = 1, while
  v₃((A − B²)²) ≥ 2; hence v₃(A² + AB² + B⁴) = 1. So exactly one of C and C + B³ is divisible by 3, and it has
  3-adic valuation v₃(A − B²) + 1. Its gcd with A − B² carries exactly v₃(A − B²).
- *Primes not dividing A − B².* They do not divide e₊ or e₋. ∎

[05: 5.1] confirms the lemma on Q_n for 2 ≤ |n| ≤ 40.

**(5.2) Corollary.** B(S₊P)·B(S₋P) = |A − B²|²/(e₊e₋) = |A − B²|, by (3.2) and (5.1). For |n| ≥ 2 this reads

    B_{n+1}B_{n−1} = |A_n − B_n²|.

**(5.3) The signed sequence.** Put h₀ = 0, h₁ = h₂ = 1 and

    h_{n+1} = −(A_n − B_n²)/h_{n−1}   (n ≥ 2),   h₋ₙ = −h_n.

So h₃ = −(7 − 1) = −6 and h₄ = −(73 − 36) = −37 [05: 5.3]. By induction from (5.2), h_n is a nonzero integer with
|h_n| = B_n for every n ≠ 0. The step is |h_{n+1}| = B_{n+1}B_{n−1}/B_{n−1}, and B₋ₙ = B_n because N preserves B.
The relation

    h_{n+1}h_{n−1} = −(A_n − B_n²)                                       (5.3)

holds for all n ≠ 0:
- for n ≥ 2 it is the definition;
- for n = ±1 both sides are 0;
- for n ≤ −2 it follows from the case −n.

Hence X(Q_n) − 1 = (A_n − B_n²)/B_n² = −h_{n−1}h_{n+1}/h_n². The marker is the reduced denominator itself, not a
multiple of it.

**(5.4) The chord product.** As an identity of rational functions on E [05: 5.4],

    (X(SQ) − 1)(X(S⁻¹Q) − 1) = (7 − X(Q))/(X(Q) − 1)².

Equivalently, with u_j = X(SʲQ) − 1, u₋₁u₀²u₁ = 6 − u₀. By (0.7) it holds at every rational Q ∉ {±G}.

**(5.5) The recurrence.** Substitute u_j = −h_{n+j−1}h_{n+j+1}/h_{n+j}² (j = −1, 0, 1) into (5.4) at Q = Q_n,
|n| ≥ 2. The numerator of the difference is h_{n+2}h_{n−2} − h_{n+1}h_{n−1} − 6h_n² [05: 5.5]. So

    h_{n+2}h_{n−2} = h_{n+1}h_{n−1} + 6h_n²   (n ≥ 2).                     (5.5)

At n = 2 this reads 0 = −6 + 6. For n ≥ 3, h_{n−2} ≠ 0 and (5.5) determines h_{n+2}. From h₁, …, h₄ = 1, 1, −6, −37
of (5.3) it gives [05: 5.5]

    (h₁, …, h₈) = (1, 1, −6, −37, 179, 7140, 11989, −8324963).

Unlike the k = 7 case there is no parity factor g_n. [05: 5.3] confirms |h_n| = B_n for n ≤ 40.

## 6. A factored stride identity

By (4.8), every B_{3m} with m ≠ 0 is divisible by 3. The classes n ≡ 1, 2 (mod 3) need a recurrence of stride 3.

**(6.1) Two chord identities.** With u_j = X(SʲQ) − 1, as identities of rational functions on E [06: 6.1]:

    u_{j−1}u_j²u_{j+1} = 6 − u_j,
    u₁ + u₋₁ = (6u₀ + 1)/u₀².

The first is (5.4) at SʲQ. The second follows from X(SQ) + X(S⁻¹Q) = (Y² + (Y + 1)²)/(X − 1)² − 2X − 2 and
Y² + (Y + 1)² = 2X³ − 1 on E [06: 6.1].

**(6.2) The stride-3 product.** Put 𝒫_m(Q) = ∏_{|j| ≤ m−1} u_j^{m−|j|}. By the first identity at j = ±1,
𝒫₃ = (6 − u₋₁)u₀(6 − u₁). The first identity at j = 0 gives u₁u₋₁ = (6 − u₀)/u₀², and the second gives the sum
u₁ + u₋₁. Hence (6 − u₋₁)(6 − u₁) = 36 − 37/u₀, and

    𝒫₃(Q) = 36u₀ − 37 = 36X(Q) − 73                                       [06: 6.2].

[06: 6.2] also checks this by multiplying the five factors directly. It vanishes at X = 73/36 = X(Q₃).

**(6.3) The stride-6 product.** Comparing exponents (6 − |j| on both sides) gives

    𝒫₆(Q) = 𝒫₃(Q)² 𝒫₃(S³Q) 𝒫₃(S⁻³Q)                                      [06: 6.3].

By (6.2) at S^{±3}Q and by (4.2), 𝒫₃(S^{±3}Q) = 36(X(R_{±3}Q) − t) with t = 73/36. For the chord through
Q₃ = (t, v),

    (X(R₃Q) − t)(X(R₋₃Q) − t) = (595/108)²(X₆ − X)/(X − t)²,   X₆ = 48833569/7140² = X(Q₆),

where 595/108 = 2v + 1 [06: 6.3]. Since 𝒫₃(Q)² = 36²(X − t)² and 36⁴(595/108)² = 7140², this yields

    𝒫₆(Q) = 48833569 − 7140²·X(Q)                                         [06: 6.3].

[06: 6.3] checks this both by the factored route and as the direct product of the eleven factors u_j^{6−|j|}.

**(6.4) The stride recurrence.** Let n ≥ 7. Then at Q = Q_n:
- For |i| ≤ 4, n + i ≥ 3, so Q_{n+i} ∉ {O, ±G} by (4.6). The composites used are S^{±1} and S^{±2} at Q_n and
  at Q_{n±3}, and S^{±3} at Q_n. Every chord map in them is applied at such a point, so (0.7) gives (6.1) and
  (6.2) at Q_n and Q_{n±3}.
- X(Q_n) ≠ t, since n ≠ ±3, and Q_n is not an exceptional point of S^{±3} = R_{±3} in (4.2). So
  R_{±3}(Q_n) = Q_{n±3}, and (0.7) gives the chord product for Q₃, hence (6.3), at Q_n.
- h_{n+j} ≠ 0 for |j| ≤ 6.

By (5.3) the products telescope [06: 6.4]:

    𝒫₃(Q_n) = −h_{n+3}h_{n−3}/h_n²,   𝒫₆(Q_n) = h_{n+6}h_{n−6}/h_n².

Eliminating X(Q_n) between (6.2) and (6.3) gives [06: 6.4]

    h_{n+6}h_{n−6} = α h_{n+3}h_{n−3} + β h_n²   (n ≥ 7),                   (6.4)

where α = 7140²/36 = 1190² = 1416100 and β = 48833569 − 73·1190² = −54541731. As in the k = 7 proof, the large
stride identity comes from two small chord identities, not from five forward and five backward substitutions.
[06: 6.4] confirms (6.4) for 7 ≤ n ≤ 400.

**(6.5) Two Somos-4 sequences.** For r ∈ {1, 2} and k ≥ −1 put s_k = s_k^{(r)} = (−1)^k h_{3k+r}. For k ≥ 2,
(6.4) at n = 3k + r ≥ 7 gives

    s_{k+2}s_{k−2} = α s_{k+1}s_{k−1} + β s_k².                              (6.5)

The sign (−1)^k cancels in each product. The entry k = 1 is checked directly from the values of h [06: 6.5], so
(6.5) holds for every k ≥ 1. The initial values, also checked in [06: 6.5], are

| k | −1 | 0 | 1 | 2 | 3 |
|---|---|---|---|---|---|
| s_k^{(1)} | 1 | 1 | 37 | 11989 | −57690006839 |
| s_k^{(2)} | 1 | 1 | −179 | −8324963 | −13536551707271 |

Only this forward range is used.

## 7. Finite 2-adic lifting

Fix r ∈ {1, 2} and write s = s^{(r)}.

**(7.1) Every term is 1 modulo 4.** α ≡ 0 and β ≡ 1 (mod 4), and s₋₁, s₀, s₁, s₂ ≡ 1 (mod 4) [07: 7.1]. By
(6.5), if four consecutive terms are ≡ 1 (mod 4), so is the next: s_{k+2} ≡ s_{k+2}s_{k−2} ≡ α + β ≡ 1. So every
s_k ≡ 1 (mod 4). In particular |s_k| = 1 forces s_k = 1.

**(7.2) The ring 𝓡.** Let 𝓡 be the set of rational functions U/V in variables z = (z₀, z₁, z₂, z₃), with
U, V ∈ Z[z] and V ∈ 1 + 2Z[z], that is, V ≡ 1 (mod 2) coefficientwise.
- *A ring.* 1 + 2Z[z] is closed under products, so 𝓡 is closed under sums and products.
- *Congruences.* For g, g′ ∈ 𝓡 write g ≡ g′ (mod 2^m) when g − g′ = 2^m·ρ with ρ ∈ 𝓡. Use the same notation
  coordinatewise for maps in 𝓡⁴.
- *Substitution.* If g ∈ 𝓡 and Ξ ∈ 𝓡⁴, then g(Ξ) ∈ 𝓡. A denominator V(Ξ) = 1 + 2V₁(Ξ), with V₁(Ξ) = U₀/V₀ ∈ 𝓡,
  equals (V₀ + 2U₀)/V₀ with V₀ + 2U₀ ∈ 1 + 2Z[z]; so its inverse is in 𝓡.
- *Congruences survive substitution.* If Ξ ≡ Ξ′ (mod 2^m), then g(Ξ) ≡ g(Ξ′) (mod 2^m). For polynomials this is the
  factorization of U(Ξ) − U(Ξ′); for a denominator use 1/V(Ξ) − 1/V(Ξ′) = (V(Ξ′) − V(Ξ))/(V(Ξ)V(Ξ′)). Also, if
  g ≡ g′ (mod 2^m) and Ξ ∈ 𝓡⁴, then g(Ξ) ≡ g′(Ξ) (mod 2^m), since (g − g′)(Ξ) = 2^m·ρ(Ξ) with ρ(Ξ) ∈ 𝓡.
- *Evaluation.* At an integer point ζ every denominator V(ζ) is odd. If U/V = U′/V′ then UV′ = U′V, so the value at
  ζ does not depend on the representation. A congruence modulo 2^m in 𝓡 becomes a congruence modulo 2^m between
  rationals with odd denominators. If g ∈ 𝓡 and Ξ ∈ 𝓡⁴, and Ξ(ζ) is again an integer point, then
  (g∘Ξ)(ζ) = g(Ξ(ζ)), because evaluating a polynomial is a ring homomorphism and all denominators met are odd.

All of this is finite algebra; no series and no completion appear.

**(7.3) Lemma (lifting at 2).** Let H ∈ 𝓡⁴ with H(z) = z + 4Δ(z), Δ ∈ 𝓡⁴. If q = 2^e·r′ with e ≥ 0 and r′ odd,
then

    H^q(z) ≡ z + 2^{e+2}r′Δ(z)   (mod 2^{e+3}).

*Proof.* KEY STEPS:
1. *Doubling.* Let J(z) = z + 2^aΛ(z) with a ≥ 2 and Λ ∈ 𝓡⁴. Then

       J(J(z)) = z + 2^aΛ(z) + 2^aΛ(z + 2^aΛ(z)).

   Since z + 2^aΛ(z) ≡ z (mod 2^a), (7.2) gives Λ(z + 2^aΛ(z)) = Λ(z) + 2^aΘ(z) with Θ ∈ 𝓡⁴. So

       J²(z) = z + 2^{a+1}(Λ(z) + 2^{a−1}Θ(z)),   with Λ + 2^{a−1}Θ ≡ Λ (mod 2).

   The condition a ≥ 2 makes 2a ≥ a + 2. This is where the prime 2 needs one more factor than the prime 5 of the
   k = 7 proof.
2. *Powers of 2.* Starting from H (a = 2) and doubling e times gives H^{2^e}(z) = z + 2^{e+2}Δ_e(z) with
   Δ_e ≡ Δ (mod 2).
3. *Odd multiples.* Put J = H^{2^e}. By induction on i, J^{i+1}(z) = J^i(z) + 2^{e+2}Δ_e(J^i(z)) ≡ J^i(z)
   (mod 2^{e+2}), since Δ_e(J^i(z)) ∈ 𝓡⁴ by (7.2); so J^i(z) ≡ z (mod 2^{e+2}). Then (7.2) gives
   Δ_e(J^i(z)) ≡ Δ_e(z) (mod 2^{e+2}), and by induction J^{r′}(z) ≡ z + r′2^{e+2}Δ_e(z) (mod 2^{2e+4}). Since
   2e + 4 ≥ e + 3 and Δ_e ≡ Δ (mod 2), the claim follows. ∎

**(7.4) The scaled map.** The four-term update of (6.5) is

    Ω(w₀, w₁, w₂, w₃) = (w₁, w₂, w₃, (α w₁w₃ + β w₂²)/w₀).

It sends the state (s_k, …, s_{k+3}) to (s_{k+1}, …, s_{k+4}) for every k ≥ −1. Put w = 1 + 4z and
χ(z) = (Ω(1 + 4z) − 1)/4. Its last coordinate is

    𝒩(z)/(1 + 4z₀),   𝒩(z) = (α + β − 1)/4 + α(z₁ + z₃) + 2βz₂ − z₀ + 4(αz₁z₃ + βz₂²)     [07: 7.4],

with α + β − 1 = −53125632 = −2⁹·3⁵·7·61. So 𝒩 ∈ Z[z] and χ ∈ 𝓡⁴. Every coefficient of
𝒩(z) − (2z₂ − z₀)(1 + 4z₀) is an integer divisible by 4 [07: 7.4]. Hence

    χ(z) ≡ L(z) := (z₁, z₂, z₃, 2z₂ − z₀)   (mod 4).

The linear recurrence c_{j+4} = 2c_{j+2} − c_j of L splits into two arithmetic progressions,
c_{j+2i} = c_j + i(c_{j+2} − c_j). So c_{j+8} = c_j + 4(c_{j+2} − c_j), that is L⁸ = I + 4M with M integral
[07: 7.4]. By (7.2) and induction, χ^j ≡ L^j (mod 4): χ(χ^j) ≡ χ(L^j) ≡ L(L^j), using both substitution rules.
Therefore

    Γ := χ⁸ = z + 4Δ(z),   Δ ∈ 𝓡⁴.

[07: 7.4] confirms χ⁸(z) ≡ z (mod 4) at 3000 random points; the coefficient statement is the argument just given.
A polynomial can vanish modulo 4 at every point without its coefficients being divisible by 4, so the
coefficientwise statement is the one used.

**(7.5) Return classes.** Put y_k = (s_k − 1)/4 ∈ Z and ζ_k = (y_k, …, y_{k+3}). Then χ(ζ_k) = ζ_{k+1}, and by the
evaluation rule of (7.2) Γ(ζ_k) = ζ_{k+8}: every intermediate point is an integer state, with odd denominator
1 + 4y = s. Since Γ ≡ z (mod 4), s_{k+8} ≡ s_k (mod 16). The table [07: 7.5], given modulo 64,

| k | −1 | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|---|
| s_k^{(1)} mod 64 | 1 | 1 | 37 | 21 | 9 | 9 | 45 | 29 | 17 | 17 |
| s_k^{(2)} mod 64 | 1 | 1 | 13 | 29 | 57 | 57 | 5 | 21 | 49 | 49 |

shows, with the period 8 modulo 16, that s_k ≡ 1 (mod 16) only for k ≡ 0 or −1 (mod 8). So s_k = 1 can occur only
in these two classes.

**(7.6) Exact lifting.** Let c ∈ {0, −1} and q ≥ 1, with q = 2^e·r′, r′ odd. Then Γ^q(ζ_c) = ζ_{c+8q}, by the
evaluation rule of (7.2). By (7.3), the first coordinate gives y_{c+8q} − y_c ≡ 2^{e+2}r′Δ₀(ζ_c) (mod 2^{e+3}).
Here Δ₀(ζ_c) = (y_{c+8} − y_c)/4 = (s_{c+8} − s_c)/16, and these numbers are odd for both sequences [07: 7.6]:

    (s₈ − s₀)/16 ≡ 9, 23   and   (s₇ − s₋₁)/16 ≡ 41, 23   (mod 64), for r = 1, 2.

Hence v₂(y_{c+8q} − y_c) = e + 2 and

    v₂(s_{c+8q} − s_c) = 4 + v₂(q)   (c = 0, −1; q ≥ 1).                     (7.6)

The right side is finite, so s_{c+8q} ≠ s_c = 1. [07: 7.6] confirms (7.6) for q ≤ 3000.

**(7.7) Conclusion.** For k ≥ 0 and r = 1, 2,

    |s_k^{(r)}| = 1   ⟺   k = 0.

*Proof.* k = 0 gives s₀ = h_r = 1. Let k ≥ 1. If k ≢ 0, −1 (mod 8) then s_k ≢ 1 (mod 16) by (7.5). Otherwise
k = c + 8q with c ∈ {0, −1} and q ≥ 1, and (7.6) gives s_k ≠ 1. By (7.1), |s_k| = 1 would force s_k = 1. ∎

The finite table located the two possible return classes; the nonzero drift excludes their entire infinite
tails. There is no search bound.

## 8. Recovery of the solutions

Let n ≥ 1.
- If 3 | n then 3 | B_n by (4.8).
- Otherwise n = 3k + r with k ≥ 0 and r ∈ {1, 2}. Then B_n = |h_n| by (5.3), and |h_n| = |s_k^{(r)}| by the
  definition in (6.5). This equals 1 only for k = 0, by (7.7).

So B_n = 1 for n ≥ 1 exactly when n ∈ {1, 2}. Since B₋ₙ = B_n, the integer points are Q_n with n ∈ {±1, ±2}, by
(4.7):

    Q₁ = (1, 0),  Q₂ = (7, −19),  Q₋₁ = (1, −1),  Q₋₂ = (7, 18).

With (x, y) = (Y, X) these are (0, 1), (−19, 7), (−1, 1), (18, 7). This proves the theorem. [08: 8.1] confirms by
direct search that there are no other solutions with y ≤ 3·10⁶.

---

## 9. Remarks

**Comparison with the k = 7 proof.** The structure is the same: ghosts, a 3-isogeny descent, chords and TS = S³T,
an exact Somos marker, a factored stride, and finite lifting. The details differ in five places.
1. *The second ring.* E′ factors over Z, so Φ-parents come from integer cube allocation. The k = 7 proof needed
   the pairs of discriminant 21, with the unit (2, 1) and a descent for the units.
2. *The obstruction.* The obstruction to a Ψ-parent is the unit ω modulo cubes, not the prime 2. It is read off
   from B mod 3 and C mod 9 (the three classes). There is no parity factor g_n. The values h₂ = 1, h₃ = −6 and
   h₄ = −37 agree with the division-polynomial values ψ₂(G), ψ₃(G) and ψ₄(G) [08: 8.4]; this is not used.
3. *The Ψ-parent.* Lemma (2.7) chooses one of the three cube roots of γ so that 3 | v. This selects the one
   Ψ-preimage, among three, that lies in the image of Φ. The other two differ from it by the rational 3-torsion
   points (0, 4) and (0, −5) of E′.
4. *The prime.* The lifting prime is 2, with period 8, instead of 5 with period 5. The doubling step (7.3) needs
   the start H = z + 4Δ rather than z + 2Δ.
5. *Why not 3.* The prime 3 divides h₃ = −6, as 5 divided h₃ = −15 for k = 7. Every s_k ≡ 1 (mod 9): the four
   initial values are ≡ 1 and α + β ≡ 1 (mod 9), so the induction of (7.1) works modulo 9 [08: 8.3]. But for every
   k < 3600 checked, the class k ≡ 5 (mod 9) of s^{(1)} (and k ≡ 3 of s^{(2)}) has s_k ≡ 1 (mod 3⁵) [08: 8.3].
   In that range, v₃(s_{5+9q} − 1) = 5 + v₃(q − 286) for q ≠ 286, and the valuation is at least 12 at q = 286
   [08: 8.3]. That is what a 3-adic zero of the interpolated sequence near q ≡ 286 (mod 3⁷) looks like, at a point
   that is not one of the known indices. If such a zero exists, no congruence modulo a power of 3 excludes the
   class. With the prime 2, the only return classes are those of the known solutions.

**Scope.** The argument proves the stated theorem. It also shows elementarily that the rational points of 243.a1
are exactly the Q_n.

**Provenance.**
- *The theorem.* It is Ljunggren's (1942), after Nagell's reduction (1921) to the Thue equation
  a³ − 3ab² + b³ = 1. That form appears here too, as the first coordinate of (u + vω)³ (1.7), but no Thue equation
  is solved.
- *The method.* It is that of the k = 7 proof.
- *Closest uses of the mechanism, all for other theorems:*
  - M. Ayad's division-polynomial scheme, with integrality as ψ_m = ±1 and congruences modulo prime powers:
    J. Number Theory 38 (1991), at p = 2; Manuscripta Math. 76 (1992); Ann. Inst. Fourier 43 (1993).
    His 1991 lemmas need ψ₃ odd. Here 2 | h₃ = −6, which is the case he left open, and the stride-3 subsequences
    handle it.
  - É. Lucas's partial descent on the model x³ + y³ = 9z³ (1873).
  - Skolem's 2-adic principle, in Ljunggren's own proof.
- *Novelty.* In the sources searched, no earlier proof of this theorem, of an equivalent form, or of a theorem
  containing it uses this combination of steps; see [literature_check.md](literature_check.md).
