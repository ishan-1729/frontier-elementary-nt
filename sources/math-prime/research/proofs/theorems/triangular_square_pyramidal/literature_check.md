# Literature check: triangular = square-pyramidal, and R1 on 8604.a1

Written 2026-10-15 by a Claude subagent for the R1 session, and edited by the author. Section 1 is from a full
reading of Uchiyama's paper, with every statement that can be
tested checked by exact arithmetic in [lit_uchiyama.py](scripts/lit_uchiyama.py) →
[lit_uchiyama.out.txt](scripts/lit_uchiyama.out.txt). Section 2 is about Avanesov, section 3 is the literature
sweep, section 4 gives the verdict, and section 5 lists what was not accessed. Section 6 (added 2026-10-26) covers
the arithmetic of Z[λ] behind [hypothesis_H.md](hypothesis_H.md). Each source is marked with how it
was read: [full], [pages], [review], [abstract], [snippet], [record] (bibliographic record only) or [not accessed].

## 1. Uchiyama (1984), read in full

**Source.** S. Uchiyama, *Solution of a Diophantine problem*, Tsukuba J. Math. 8 (1984), no. 1, 131–157, received
29 July 1983, DOI 10.21099/tkbjm/1496159951, Zbl 0544.10011 (review by K. Szymiczek), MR 86i:11010 (per Guy) [full].
- *Pages: 131–157.* The printed header of p. 1 reads 131–157, the scan has 27 pages, the last page is numbered 157,
  and zbMATH and Avanesov–Gusev (2012, ref. [6]) agree. Guy's *UPINT* (2nd and 3rd eds., D3) prints 131–137, a
  misprint; [rank2_signatures_2026-10/literature.md](../../targets/rank2_signatures_2026-10/literature.md) copied it.
- *Access.* Project Euclid showed a robot check on 2026-10-15, which was not passed. The scan with OCR text layer in
  the University of Tsukuba repository (https://tsukuba.repo.nii.ac.jp/records/15881) was fetched with WebFetch,
  and its 27 page images were rendered in the session scratchpad and read page by page. Nothing is quoted beyond
  equations.
- *Context.* He cites Guy's first edition (1981), where the problem is listed as open, and Avanesov (Acta Arith.
  1966/67) for the tetrahedral problem. He says his reasoning is alike in principle to Avanesov's but differs much in
  detail. He cites neither Avanesov 1971 nor Finkelstein–London 1972 (section 2, 3).

**1.1 Statement and model.**
- *Theorem.* The integer solutions of (4) 3Y² = X³ − X + 3 are X ∈ {−1, 0, 1, 3, 8, 11, 13, 171, 1704}. With
  X = 2x+1, Y = 2y+1, (4) is t_y = P_x. The odd X give x ∈ {−1, 0, 1, 5, 6, 85}, and the common values 1, 55, 91,
  208335.
- *Our curve.* 27(3Y² − X³ + X − 3) = (9Y)² − ((3X)³ − 9(3X) + 81). So his points are the points (3X, 9Y) of
  E: y² = x³ − 9x + 81, and they are exactly the points of E with 3 | x. His nine X are the nine x ∈ {−3, 0, 3, 9, 24,
  33, 39, 513, 5112} of our list divided by 3. The other three x (−5, 7, 1099) form our trivial class and are not on
  his model.

**1.2 The cubic field (§1.1).**
- K = Q(λ), λ³ − λ + 3 = 0. Our θ is 3λ. He states:
  - D_K = −239 and the integral basis {1, λ, λ²} (true: disc(g) = −239 is squarefree);
  - the fundamental unit ε = λ² + λ − 1, taken from the Delone–Faddeev table [7, p. 113] (N(ε) = 1; fundamentality not
    re-checked here; update 2026-10-26: [hypothesis_H.md](hypothesis_H.md) (4.1) proves it);
  - the class number h_K = 1, stated without proof or reference;
  - the prime ideal factorizations (6) (3) = (λ)(λ−1)(λ+1) and (7) (239) = (3λ²−1)²(3λ²−4) (true, as element
    identities up to units).
- (8): (X − λ)(X² + λX + λ² − 1) = ζλ(λ−1)(λ+1)Y² with ζ a unit (it is −1).
- *Coprimality.* The two factors have the same common divisors as X − λ and 3λ² − 1, a prime of norm 239. This
  common divisor is nontrivial iff X ≡ 124 (mod 239). Then 239 divides g(X) exactly once, since g(124) = 3·239·2659 and
  g′(124) = 193·239, so 3Y² = g(X) is impossible. All checked.
- *Factorization.* With h_K = 1, X − λ = ±π ε^k (a + bλ + cλ²)² with π = λ, λ − 1 or λ + 1 as X ≡ 0, 1 or −1 (mod 3),
  and k ∈ {0, 1}. This step uses ideals and the class number.

**1.3 The reduction (§§1.2–1.5).**
- (9) gives the square (a + bλ + cλ²)² = A + Bλ + Cλ². Equating coefficients gives one ternary system for each case.
- For k = 1, and for the wrong sign when k = 0, the system is impossible modulo 3, 4 or 8. Checked by exhaustive
  residues.
- The case X ≡ 0, 1, −1 (mod 3) with k = 0 is exactly our class (1,1), (1,0), (0,1). Indeed x − θ = 3(X − λ),
  −3λ = −θ, −3(λ−1) = 3 − θ and −3(λ+1) = −3 − θ.
- **So his §1 proves hypothesis (H) for the points with 3 | x**, from h_K = 1 and the unit of K. That is all that
  T(n) = P(m) needs, since there x = 6m + 3. It uses ideal theory and a class number, which this project excludes.
- The systems are (10) 2b(a+c) = 3c², (a+c)² + b² − 6bc = 1; (16) (2a−b+c)(b−c) = 3c², a² − 6bc = 1; and (25)
  (2a+b+c)(b+c) = 3c², a² − 6bc = 1.
- Coprime factorization of the first equation of each system, split by the power of 2 in c, gives sixteen
  parametrizations and twelve distinct equations. All sixteen were checked symbolically.
- The sub-cases bc = 0 give X = 0, 3 (class (1,1)), 1, 13 (class (1,0)) and −1, 8, 11 (class (0,1)).

**1.4 His twelve sub-case equations are our twelve equations F_cl(m,n) = c.** Each relation F_cl(T(U,V)) = s·L(U,V)
was found from the roots and verified exactly (§5 of the output). Under his parity conditions the gcd of T(U,V) is 1,
which fixes c.

| his sub-case | his equation L = r | our class | c | (m, n) = T(U,V) | his solutions → (m,n) → X |
|---|---|---|---|---|---|
| 1.2(i) | 9U⁴ − 24UV³ + 4V⁴ = 1, U odd | (1,1) | 9 | (V, 3U), s = 9 | none: excluded mod 16 |
| 1.2(ii) | (11) U⁴ − 72UV³ + 36V⁴ = 1 | (1,1) | 1 | (V, U) | (1,2) → (2,1), X = 1704; (1,0) → (0,1), X = 0 |
| 1.2(iii) | (12) 36U⁴ − 12UV³ + V⁴ = 1 | (1,1) | 36 | (V, 6U), s = 36 | (0,1) → (1,0), X = 3 |
| 1.2(iv) | (13) 4U⁴ − 36UV³ + 9V⁴ = 1 | (1,1) | 4 | (V, 2U), s = 4 | (2,1) → (1,4), X = 171 |
| 1.3(i) | 9U⁴ − 18U²V² + 24UV³ + V⁴ = 4, U, V odd | (1,0) | 36 | (V, −3U), s = 9 | none: excluded mod 16 |
| 1.3(ii) | U⁴ − 18U²V² + 72UV³ + 9V⁴ = 4, U, V odd | (1,0) | 4 | (V, −U) | none: excluded mod 16 |
| 1.3(iii),(v) | (17) 9U⁴ − 18U²V² + 24UV³ + V⁴ = 1 | (1,0) | 9 | (V, −3U), s = 9 | (0,1) → (1,0), X = 13 |
| 1.3(iv),(vi) | (18) U⁴ − 18U²V² + 72UV³ + 9V⁴ = 1 | (1,0) | 1 | (V, −U) | (1,0) → (0,−1), X = 1 |
| 1.4(i) | (26) 9u⁴ + 18u²v² − 24uv³ + v⁴ = 4, u, v odd | (0,1) | 36 | (v, 3u), s = 9 | (1,1) → (1,3), X = 8 |
| 1.4(ii) | u⁴ + 18u²v² − 72uv³ + 9v⁴ = 4, u, v odd | (0,1) | 4 | (v, u) | none: excluded mod 128 |
| 1.4(iii),(v) | (27) 9U⁴ + 18U²V² − 24UV³ + V⁴ = 1 | (0,1) | 9 | (V, 3U), s = 9 | (0,1) → (1,0), X = 11 |
| 1.4(iv),(vi) | (28) U⁴ + 18U²V² − 72UV³ + 9V⁴ = 1 | (0,1) | 1 | (V, U) | (1,0) → (0,1), X = −1 |

- The four congruence exclusions are our (class, c) = ((1,1), 9), ((1,0), 36), ((1,0), 4) and ((0,1), 4). Our
  02_exclusion.py also removes the first, second and fourth by a congruence alone. The third, ((1,0), 4), has the
  field-form solution G(1,0) = 4, which gives no class-form solution.
- Remaining are his eight Thue equations (11), (12), (13), (17), (18), (26), (27), (28). In §2 he replaces (13) by
  (34) u⁴ − 92u³v + 96u²v² − 32uv³ + 4v⁴ = 1 (U = 2u − v, V = u), and (26) by (39)
  U⁴ − 144U²V² + 936UV³ − 1692V⁴ = 1.

**1.5 The quartic fields and units (§2.1 and Propositions 1–8).**
- *The fields.* Every defining polynomial has discriminant −2¹²·3⁶·239 and two real roots. The resolvent cubic has no
  rational root, so the Galois group is S4. All checked.
- *Integral bases* come from Berwick's algorithm [1]. *Fundamental units* come from Berwick's criterion [2]: three
  units with the sign and size conditions (a), (b), (c), any two of which are fundamental. Tzanakis (Math. Comp. 46,
  1986) later filled a gap in this theorem (section 3). His data for Proposition 1 satisfy the sign and size
  conditions, with the smaller real root; the "least" clause was not checked.

In the last column ε₁, ε₂ are core.UNITS of the same class.

| his field | Props. | his polynomials | our field | his D | true D_L | his units in core.UNITS (mod ±1) |
|---|---|---|---|---|---|---|
| L_A | 1, 2, 3 | x⁴ − 72x + 36 (= our g of (1,1)), x⁴ − 12x³ + 36, x⁴ − 92x³ + 96x² − 32x + 4 | class (1,1) | −2⁸·3²·239 | −2⁶·3²·239 | Prop. 1: α = −ε₁, β = ε₂ |
| L_B | 4, 5 | x⁴ − 18x² + 72x + 9, x⁴ + 24x³ − 18x² + 9 | class (1,0) (θ = 1 − 2β) | −3²·239 | −3²·239 | Prop. 4: α = ε₁ε₂⁻¹, β = −ε₂ |
| L_C | 6, 7, 8 | x⁴ + 18x² − 72x + 9 (= our g of (0,1)), x⁴ − 24x³ + 18x² + 9, x⁴ − 144x² + 936x − 1692 | class (0,1) | −2⁸·3²·239 | −2⁶·3²·239 | Prop. 6: α = −ε₁ε₂, β = ε₁⁻¹ |

- *The unit groups agree.* In each field his fundamental pair and core.UNITS generate the same group modulo ±1.
  01_coverage.py proves core.UNITS fundamental as a by-product, so his three unit groups are right. This holds for
  the printed units of Props. 1, 2, 4, 5 and 6, and for corrected units of Props. 3, 7 and 8 (section 1.7).
- *The argument.* In each field he writes ±(U − Vθ) = αᵃβᵇ. He then excludes all but the known exponent pairs (a, b)
  by congruences modulo powers of 2 on the coordinates ω₂, ω₃ (the "markers"), split by the 2-adic valuations of a, b
  or of combinations such as a − 6b. There are no p-adic series and no completions; everything is a finite congruence.
  - **Prop. 1 (11)** is printed in full. α = ω₂ = ε₁⁻¹, β = ε₁ε₃, Ω = 1 − 2ω₁ = −α². The steps are:
    1. a ≡ b (mod 2), from a congruence mod 4;
    2. a, b odd is excluded by ω₃ mod 2⁵;
    3. for a = 2m, b = 2n: ω₂ ≢ 0 (mod 2⁵) unless m + n ≡ 0, 1 (mod 4), and then n ≡ 0 (mod 4);
    4. n = 2ʰn′ (h ≥ 2) is excluded by ω₃ mod 2^{h+5};
    5. b = 0 leaves m = 0, 1, by ω₂ mod 2^{h+4} and mod 2^{h+3}.

    The solutions are (a, b) = (0, 0) and (2, 0), that is (U, V) = (1, 0) and (1, 2).
  - **Prop. 2 (12)** is printed. It uses the single marker ω₃: α^{2ʰm} ≡ ±(1 + 2ʰmω₃) (mod 2^{h+1}), the same for β,
    and a formula for odd exponents.
  - **Prop. 4 (18)** is printed in part. α has period 7 modulo 2, which gives a + b ≡ 0 (mod 7). Then
    a − 6b = 2ʰ·7e and b = 2ᵏf, with α⁷ ≡ −1 + 4ω₂ and α⁶β ≡ 1 − 2(ω₁ − ω₂) (mod 8), give "a − 6b = 0 or b = 0".
    "Similarly" 6a − b = 0 or a = 0, and hence a = b = 0.
  - **Props. 3, 5, 6, 7, 8** print only the field, the integral basis, the matrix M between bases (for 3, 5, 7, 8),
    the units and a few congruences mod 8, then refer to "the same lines as" Prop. 2 or Prop. 4.
- *§2.10.* No two of the eight forms are GL₂(Z)-equivalent. True: no relation with |det T| = 1 and s = 1 exists.

**1.6 What is proved in print, and what is only asserted.**
- *Proved in print, and correct:* the reduction of §1, given h_K = 1 and the fundamental unit of K. Every printed
  identity, parametrization and congruence exclusion of §1 checks. There are two harmless misprints: p. 137 has
  bc = uv(u ± v) for uv²(u ± v), and p. 138 has 3v³ for 3v².
- *Taken from tables or asserted:* h_K = 1 (no source given); the unit ε of K (Delone–Faddeev table); the integral
  bases (Berwick's algorithm, no computation shown); the fundamental units (Berwick's criterion; data printed only for
  Prop. 1, and the criterion had a gap later filled by Tzanakis 1986); and Props. 3, 5, 6, 7 and 8 beyond their
  first congruences.
- *Printed but not correct as printed* (section 1.7): steps of Props. 1, 2 and 4.

**1.7 Errors found** (exact; [lit_uchiyama.out.txt](scripts/lit_uchiyama.out.txt), with an independent sympy
re-computation in its §9). In all, 117 of the 155 tested statements are true.
- **Field data.**
  - *Integral bases.* For L_A and L_C (Props. 1, 2, 3, 6, 7, 8) his "integral bases" span orders of index 2 in the
    maximal order. The field discriminant is −2⁶·3²·239 = −137664, not −2⁸·3²·239. For example, (1/12)θ³ + (1/2)θ is
    integral in the field of x⁴ − 72x + 36. 00_fields.py finds the same discriminants with sympy's round_two. L_B is
    right.
  - *Prop. 3.* The relation between the roots is (4 − ξ*)ξ = 2, not (4 − ξ)ξ* = 2. In M, row 1 begins 176, not
    186. The printed β = −1462 − 1073ω₁ + 2253ω₂ − 24ω₃ was computed with the wrong M and has norm 1680961, so it is
    not a unit. The unit is β = −1422 − 1073ω₁ + 2253ω₂ − 24ω₃.
  - *Prop. 7.* Its M is right, but the printed α, β have norms −357444935 and −73968435431, so neither is a unit. The
    units of Prop. 6, carried over by M, are α = −257 − 61ω₁ + 483ω₂ − 20ω₃ and β = −101 + 17ω₁ + 229ω₂ − 10ω₃.
  - *Prop. 8.* Rows 2 and 3 of M are wrong: they should be (56, −20, 5, 1) and (18, −13, 9, 1). The printed α, β were
    computed with the wrong M and are not units. The correct ones are α = −428 + 176ω₁ − 63ω₂ − 10ω₃ and
    β = 6 − 2ω₁ + ω₂.
  - *(39).* The printed substitution U = x − 6y should be U = x + 6y.
  - *False congruences mod 8.* These are, even up to sign:
    - Prop. 5: α⁷ and α⁶β. They repeat the printed β⁷ and αβ⁶; the true values are α⁷ ≡ −1 + 4ω₂ and
      α⁶β ≡ −1 + 2ω₁ + 2ω₂;
    - Prop. 6: α and β², where the sign of the ω₂-coefficient is wrong;
    - Prop. 7: α², β² and αβ;
    - Prop. 8: α² and αβ.
- **Arguments** (the statements of all eight propositions are true: our 02_exclusion.py certifies the twelve
  equations).
  - *Prop. 1, step 4 is not proved as printed.*
    - The congruence for ω₃ modulo 2^{h+5} is false: it fails for 2328 of 3120 tested triples (h, m, n).
    - Its conclusion, ω₃ ≢ 0 (mod 2^{h+5}), also fails. It fails for every tested m ≡ 1 (mod 2^{h+2}), and for half
      of the tested m with v₂(m) = h − 1.
    - These are the 2-adic neighbourhoods of the solutions (a, b) = (2, 0) and (0, 0). For example, ω₃(α²β⁸) =
      96184930944 has 2-adic valuation exactly 7 = h + 5.
    - No solution is lost, because ω₃ ≠ 0 there. But the printed congruence does not show it.
    - Among the failures, the 2-adic valuation of ω₃ reaches h + 13 in the tested range. So working modulo a fixed
      higher power of 2 is not shown to suffice. The natural repair uses both markers at the known zeros, as a
      Jacobian criterion does (our Proposition 17).
    - Two side congruences are misstated: the 2ʰ-congruence for Ω is false for h = 1 (it is used only for h ≥ 2), and
      the sign (−1)^{2ʰ−1} should be (−1)^{2^{h−1}n}. The coefficients that step 2 and step 3 use are right.
  - *Prop. 2: the odd case and the case v(a) = v(b) are not proved.* For odd a, b, the claim
    αᵃβᵇ ≡ ±(1 + 2c(a,b)ω₂ + 2ω₃) (mod 4) is false whenever a + b ≡ 0 (mod 4); then αᵃβᵇ ≡ −1 (mod 4). For example,
    α⁻¹β = 3 − 4ω₁ + 8ω₂ exactly: ω₃ = 0, and both markers are ≡ 0 (mod 8). Only ω₂ = 8 ≠ 0 (mod 16) excludes it.
  - *Prop. 4.*
    - The lemma (α⁶β)^{2ᵏf} ≡ 1 − 2^{k+1}f(ω₁ − ω₂) (mod 2^{k+2}) is false for k ≥ 1, because (ω₁ − ω₂)² ≢ 0
      (mod 2).
    - The inference "f ≡ 0 (mod 2) if k = 1" therefore has no basis. In fact ω₂ ≡ 4 (mod 8) in that case, which
      excludes it for another reason.
    - The case where a and b are both even is never treated; "if k > 0 then h = 0" is asserted.
    - The half "6a − b = 0 or a = 0" is not printed.
  - *Props. 3, 5, 6, 7, 8* inherit these gaps, since they say "as in Proposition 2" or "as in Proposition 4". On the
    tested box, the Prop. 2 valuation pattern holds for min(ω₂, ω₃) in Props. 6 and 7 and fails in Props. 3 and 8.
  - *Prop. 8 cannot be completed by congruences modulo powers of 2 alone* (found in round 2 of the review). Its
    equation (39) is a transform of (26), our class (0,1) with c = 36 (§1.4). The markers of that orbit have a common
    2-adic zero that is not an integer zero ([proof.md](proof.md) §4.5; Hensel certificate in
    [check_2adic_zero.out.txt](scripts/check_2adic_zero.out.txt)). Whether ε^{(a,b)}w lies in Q + Qβ depends neither on
    the basis nor on the generators chosen for the unit group, so the same 2-adic point meets Uchiyama's conditions,
    with his units corrected: no argument by congruences modulo powers of 2 can exclude it, and odd moduli are needed
    (43 and 13 in 02). Prop. 3, by contrast, is repaired 2-adically: its orbit is our (1,1), c = 4, which Proposition 17
    at q = 16 and a sieve by powers of 2 finish ([verify_exclusion.out.txt](scripts/verify_exclusion.out.txt), other
    modulus sets).
- **Summary.** The theorem is right: Theorem A of this card with Uchiyama's §1 (given h_K = 1 and the unit of K)
  proves it; within the card it rests on (H). The reduction of §1 is a sound proof, given h_K = 1 and the unit of K.
  The solution of the eight Thue equations is not complete in print, and for Prop. 8 it cannot be completed by
  congruences modulo powers of 2.
  - The zbMATH review (Szymiczek; read through the zbMATH API) describes the proofs as complete and fairly detailed.
    That is not what the page-by-page check finds.
  - The confirmations in lit_uchiyama.out.txt (marker zero sets for |a|, |b| ≤ 20, and box searches with |U|, |V| ≤
    1000) agree with his solution lists.

## 2. Avanesov

- **Guy, UPINT, 3rd ed. (2004), D3** [full; local text claude/sources/guy/guy_all.txt, lines 10180–10430].
  - Martin Gardner (Sci. Amer., July 1974) noted that "triangle = square pyramid" was the one open case among the
    pairings of four figurate families.
  - Guy says Schinzel sent Avanesov's affirmative answer, which Uchiyama later found again.
  - Guy's references: È. T. Avanesov, *The Diophantine equation 3y(y+1) = x(x+1)(2x+1)* (Russian), Volž. Mat. Sb.
    Vyp. 8 (1971) 3–6, MR 46 #8967; and Uchiyama with the wrong pages 131–137.
- **The 1971 paper** [not accessed].
  - It is not in zbMATH; the author search through the zbMATH API was checked on 2026-10-15.
  - Avanesov–Gusev (2012) cite it as «Диофантово уравнение 3y(y+1) = x(x+1)(2x+1)», Волжский математический сборник,
    1971, № 8, 3–6, and give it, with Uchiyama, as the source of the result.
  - Its method is unknown here. Four pages is short for an argument of Uchiyama's kind, unless it rests on earlier
    work such as the 1968 note.
- **The 1968 note** [record].
  - È. T. Avanesov, «К вопросу об одной теореме Сколема» (A note on a theorem of Skolem), Izv. Akad. Nauk Armyan.
    SSR, Mat. 3 (1968), no. 2, 160–165, Zbl 0165.06002 (no review).
  - The Armenian NAS digital library has the full text (arar.sci.am, record 110726), but access is restricted.
  - ProofWiki ("Square Pyramidal and Triangular Numbers") cites it for the result [full: the page has no proof and
    flags one as wanted].
  - Finkelstein–London (1972) list it among their four references (section 3), which suggests a Skolem-type tool. So
    the project's "Avanesov 1968/1971" is uncertain on the 1968 side: whether that note contains the result is not
    known.
- **Avanesov–Gusev, Chebyshevskii Sb. 13:2 (2012) 5–7** [full; mathnet.ru cheb28].
  - This is a survey without proofs.
  - For this problem it gives the change of variables to 3N² = M³ − M + 3 and u² = v³ − 9v + 81, which is 8604.a1.
  - It then says that Dem'yanenko's estimate bounds |x| below 940697295604 and that a computer calculation gives the
    result.
  - So a third route is in print: a height bound from Mordell–Weil generators, then a search. It is not said
    whether this was Avanesov's 1971 method. The survey uses the same route for the cannonball problem and cites its
    2007 paper for it.
- **Avanesov–Gusev, Eurasian Math. J. 2007, no. 3, 9–12** [not accessed].

## 3. Other proofs and related literature
A sweep on 2026-10-15 used WebSearch, WebFetch, the zbMATH API and the built-in browser for abstracts.

**Complete solutions of the theorem in print.**
- *Avanesov (1971)* [not accessed]; section 2.
- *R. Finkelstein, H. London, On triangular numbers which are sums of consecutive squares*, J. Number Theory 4
  (1972) 455–462, DOI 10.1016/0022-314X(72)90036-4, Zbl 0244.10015 [abstract; references via Crossref].
  - The abstract says they solve H. E. Thomas Jr.'s problem (Amer. Math. Monthly 75 (1968) 1018, Problem 5634): all
    integer solutions of Σ_{i≤n} i = Σ_{i≤r} i², that is, T(n) = P(r). The solution list is not legible in the online
    abstract.
  - Their references are Avanesov 1968, Billevich 1956 (units of cubic and quartic fields), Delone 1964 (the
    Delone–Faddeev monograph) and Thomas 1968. This toolkit matches Uchiyama's cubic-field descent, quartic units and
    Skolem. Not confirmed: the full text is behind a robot check, which was not passed.
  - **No project file mentions this paper before now.**
- *Uchiyama (1984)*, section 1: complete reduction, incomplete Thue step as printed.
- *Avanesov–Gusev (2012)*: Dem'yanenko bound plus computer search, described without details (section 2).
- *Machine computations of the integral points of 8604.a1.* LMFDB lists the 24 points
  ([lmfdb_data.md](../../targets/rank2_signatures_2026-10/lmfdb_data.md)). Such lists come from integral-point
  routines based on elliptic logarithms and a Mordell–Weil basis. LMFDB's reliability page, quoted in lmfdb_data.md,
  says that its integral points were computed rigorously by two independent implementations (Magma and SageMath)
  whose results were compared; they rest on LMFDB's Mordell–Weil basis. No paper was found that treats this curve by
  elliptic logarithms explicitly.

**Sources that state the result without a proof.**
- OEIS A039596 [full] lists 0, 1, 55, 91, 208335 and links Finkelstein–London and Thomas's problem.
- ProofWiki [full] has no proof.
- MathWorld, Square Pyramidal Number [full] cites Guy only.
- Guy D3 [full].
- Deza–Deza, *Figurate numbers* (2012), p. 97; Roberts, *Lure of the Integers*, p. 245; Wells' dictionaries
  [not accessed; cited by OEIS].
- Krawczyk, *On some Diophantine equations connected with triangular and pyramidal numbers*, Pr. Nauk. Inst. Mat.
  Politech. Wrocław 19 (1980) 49–53, Zbl 0453.10012 [record]. Its content is unknown.

**Wider statements that do not determine the solutions.**
- B. Brindza, Á. Pintér, S. Turjányi, *On equal values of pyramidal and polygonal numbers*, Indag. Math. 9 (1998)
  183–185, Zbl 0920.11015 [review]. It proves effective finiteness for pyramidal = polygonal, for all orders except
  effectively determinable pairs. That covers this case as a finiteness statement only.
- M. Kaneko, K. Tachibana, *When is a polygonal pyramid number again polygonal?*, Rocky Mountain J. Math. 32 (2002)
  149–165, Zbl 1039.11012 [review]. It finds families of solutions over function fields and does not treat this
  case.

**Elementary proof.** None was found. Searches covered the equation in English and Russian, figurate-number
coincidences, OEIS, ProofWiki and MathWorld.

**Related methods.**
- N. Tzanakis, *A remark on a theorem of W. E. H. Berwick*, Math. Comp. 46 (1986) 623–625, Zbl 0588.12005
  [abstract]. It indicates and fills a gap in Berwick's theorem on fundamental units of semireal biquadratic fields:
  the criterion Uchiyama used in 1984.
- R. J. Stroeker, N. Tzanakis, *On the application of Skolem's p-adic method to the solution of Thue equations*,
  J. Number Theory 29 (1988) 166–195 [abstract]. It treats totally real normal quartic fields; not our S4 fields.
- R. J. Stroeker, *On quartic Thue equations with trivial solutions*, Math. Comp. 52 (1989) 175–187 [abstract].
  Quartic fields of negative discriminant with a real quadratic subfield; our fields are S4 and have none.
- L. Hajdu, Á. Pintér, Sz. Tengely, N. Varga, *Equal values of figurate numbers*, J. Number Theory 137 (2014)
  [full, revised preprint]. It does not treat this pair.
- For the two-index Skolem tradition (Skolem 1933, Ljunggren, Avanesov 1966/67, Bremner–Tzanakis) and for the
  provenance of the coverage lemma (Koecher's reduction theory), see
  [triangular_tetrahedral/literature_check.md](../triangular_tetrahedral/literature_check.md), §2. It applies here
  unchanged.

## 4. Verdict
- **Theorem.** Classical.
  - Avanesov (1971) and Finkelstein–London (1972) solved it; Finkelstein–London answer Thomas's 1968 Monthly problem
    (Problem 5634). Neither was read here.
  - Uchiyama (1984) re-solved it; his reduction is complete and his Thue step is incomplete as printed.
  - Avanesov–Gusev (2012) report a computer route.
  - Modern integral-point software reproduces the 24 points.
- **Mechanism.** It is Uchiyama's: a descent in the cubic field Q(θ), then quartic Thue equations in the same three S4
  fields of signature (2,1), then units, then a two-index exclusion by finite congruences.
  - His twelve sub-case equations are exactly our twelve equations F_cl = c.
  - His three fields are ours, and his unit groups are ours.
  - Finkelstein–London very probably used the same mechanism, judging by their references; not confirmed.
  - By the project's criterion (novel unless the same proof, or a rewording of it, has proved this theorem), **no
    novelty is claimed for the mechanism.**
- **What differs.**
  - *Reduction.* His uses ideals, h_K = 1 and a tabulated unit, and is unconditional for the points with 3 | x.
    Ours uses no ideals or class numbers. It first rested on hypothesis (H), which [hypothesis_H.md](hypothesis_H.md)
    now proves for all integral points (update 2026-10-26).
  - *Units.* He takes them from Berwick's criterion. We prove coverage with Lemma 3.2's cone certificate, without the
    unit theorem; fundamentality is a by-product.
  - *Exclusion.*
    - His is a case analysis on 2-adic valuations of the exponents, mostly with one marker at a time. As printed it
      fails near the known solutions and in Prop. 2's odd case, and five propositions are only sketched. For Prop. 8
      no argument modulo powers of 2 can succeed (§1.7).
    - Ours applies Proposition 17 (a finite Jacobian criterion, Skolem's lemma in finite form) at each known zero,
      with a congruence sieve over the moduli 4, 7, 8, 13, 16, 23, 43, 49.
  - *Completeness.* All twelve equations have exact, re-runnable certificates.
  - The finite, completion-free shape itself is Uchiyama's, as it was for 1611.a1.
- **What this card adds is a complete and checked proof.** It also records the errors in Uchiyama's printed
  argument (section 1.7).
- **A pointer for the open question on (H).** Uchiyama's §1 shows that, for the points needed here (x ≡ 3 mod 6),
  (H) follows from the factorization in Z[λ] with h_K = 1 and the unit ε.
  - The field here (discriminant −239) is norm-Euclidean, with Euclidean minimum 8/9, in the table of S. Cavallar and
    F. Lemmermeyer, *The Euclidean algorithm in cubic number fields* (Number Theory, Eger 1996; de Gruyter 1998;
    arXiv:1202.6019), as read by the literature check of the 1611.a1 session; not re-checked here.
  - A sibling session reports an ideal-free proof of (H) for 1611.a1, whose field is not norm-Euclidean, by a
    division with multipliers {1, 2} and a one-unit descent
    ([hypothesis_H.md](../triangular_tetrahedral/hypothesis_H.md), 2026-10-15; not checked here).
  - So plain Euclidean division may give Uchiyama's §1 without ideals here, and make Theorem B unconditional for
    T(n) = P(m). This is a lead, not a result.
  - *Update 2026-10-26.* The lead is carried out: [hypothesis_H.md](hypothesis_H.md) proves (H) for all integral
    points by a Euclidean division in Z[λ] with an exact certificate of 116 boxes, and a one-unit descent, so Theorem
    B is unconditional. The Cavallar–Lemmermeyer row was re-read in review round 3 (section 6).
- Not found is not the same as does not exist. An elementary proof, the method of Avanesov 1971, and the details of
  Finkelstein–London 1972 were not found or not read.

## 5. Not accessed
- Avanesov, Volzh. Mat. Sb. 8 (1971) 3–6 (no online copy found) and its review MR 46 #8967 (MathSciNet returned
  HTTP 403).
- Avanesov, Izv. AN Arm. SSR Mat. 3 (1968) 160–165 (arar.sci.am record 110726, access restricted).
- Avanesov–Gusev, Eurasian Math. J. 2007, no. 3; and Avanesov–Gusev–Borzhenskaya, *Computers in number theory*
  (Pyatigorsk 2008).
- Finkelstein–London 1972, full text (robot check not passed); Thomas, Amer. Math. Monthly 75 (1968) 1018 and its
  solution.
- Krawczyk 1980; Deza–Deza 2012; Wells; Roberts; Gardner, Sci. Amer. 231 (July 1974) (known only through Guy).
- Uchiyama's MathSciNet review MR 86i:11010; Berwick 1927 and 1932; Billevich 1956; the Delone–Faddeev table.
- Mordell, *Diophantine Equations* (1969), Chs. 27–28 (cited by Uchiyama for finiteness only).

## 6. Hypothesis (H): the arithmetic of Z[λ], λ³ = λ − 3 (added 2026-10-26)
For [hypothesis_H.md](hypothesis_H.md). The literature agent of review round 3 ran this check with web search, page
fetches and the zbMATH Open API; the author edited it. PDFs marked [full] or [pages] were read through their text
layer, with table rows rebuilt from word positions. "Exact" marks an identification checked in sympy by
[lit_H_lemmermeyer.py](scripts/lit_H_lemmermeyer.py) or [lit_H_lmfdb.py](scripts/lit_H_lmfdb.py) (outputs beside them).

**6.1 The Euclidean property.**
- S. Cavallar, F. Lemmermeyer, *The Euclidean algorithm in cubic number fields*, Number Theory (Eger 1996), de Gruyter
  1998, 123–146; arXiv:1202.6019v1 [full, arXiv version]. The table of complex cubic fields (§7, p. 16) has the row
  "−239 E 8/9": norm-Euclidean, with Euclidean minimum M(K) = 8/9 and no second minimum given. The round-3 referee read
  the same row.
  - Method: a floating-point covering with k ≤ 0.99c, then exact integer arithmetic at the exceptional points.
  - The authors trust the results through agreement with earlier work; they do not present them as proofs.
- F. Lemmermeyer, *Euclidean Rings* (Diplomarbeit, Heidelberg 1989), arXiv:2608.23216 [pages 5–7, 91–95, 101, 145–147,
  246]. Row 25 of his table: −239, E, x³ − x + 3, h = 1, the unit 5 + 3λ = −ε² (exact), M(K) = 8/9.
- **First proof.** E. M. Taylor, *Euclid's algorithm in cubic fields with complex conjugates*, PhD thesis, Royal
  Holloway 1975 [abstract]; J. London Math. Soc. (2) 14 (1976) 49–54, Zbl 0358.12002 [review: lists for −D,
  0 < D < 5000].
  - The attribution is Lemmermeyer's (1989, p. 101): Godwin (1967) proved the fields with −152 ≤ d ≤ −23
    norm-Euclidean, and Taylor (1976) determined all norm-Euclidean fields with |d| < 680.
  - Taylor's text was not read: the scan is 17.7 MB, above the fetch limit, and the journal page returned 403.
- Not about −239: J. R. Smith (1969, 1971) and Godwin–Smith (1993), totally real fields; L. E. Clarke (1951);
  D. A. Clark (1994, 1996), weighted norms; V. G. Cioffari (1979), pure cubic fields.
- −239 is in the conjectured list of the 58 norm-Euclidean complex cubic fields: Cavallar–Lemmermeyer, *Euclidean
  windows*, LMS J. Comput. Math. 3 (2000), p. 352 [pages]; OEIS A187773 [full].
- **Verdict.** Remark 2.5 of hypothesis_H.md proves the lower bound 8/9 exactly, at the point ξ₀ fixed by ε, and it
  agrees with the recorded minimum. The table of hypothesis_H.md §2 proves M ≤ 32733/32768 < 1 in exact arithmetic. No
  exact certificate of Taylor's theorem was found in print.

**6.2 Class number and unit.**
- LMFDB 3.1.239.1 [full, 2026-10-26]: defining polynomial x³ − x − 3, whose root is −λ (exact); class number 1;
  fundamental unit x² − x − 1, which is ε (exact). The page does not say whether its class number is proven.
- Tables: Delone–Faddeev 1940, Uchiyama's source for ε; I. O. Angell, Bull. London Math. Soc. 5 (1973) 37–38. Neither
  was accessed; both are known through Uchiyama and Lemmermeyer.
- K. K. Billevich, Mat. Sb. 40 (1956) 123–136, gives a method for units [record only].
- hypothesis_H.md needs none of these: Sections 2–3 give the divisibility theory, and (4.1) proves that ε is
  fundamental.

**6.3 Earlier treatments of (H).**
- The only factorization found is Uchiyama's §1 (section 1 above): for the points with 3 | x, with ideals, h_K = 1
  without a source, and the tabulated ε. His congruences modulo 3, 4 or 8, which remove k = 1 and the wrong sign
  (section 1.3), do the work of the norm condition (5.5) and the step modulo 4 (5.7) of hypothesis_H.md.
- No source was found that uses Taylor's theorem, or any other ideal-free argument, for 3Y² = X³ − X + 3,
  y² = x³ − 9x + 81 or T(n) = P(m). Avanesov 1971 and Finkelstein–London 1972 were not read (section 5); the references
  of Finkelstein–London (Avanesov 1968, Billevich 1956, Delone 1964, Thomas 1968, per Crossref) point to tabulated
  invariants. Krawczyk 1980 is known only from its record ("reduction to Pell equation").
- **Verdict.** No earlier ideal-free, table-free proof of (H) was found in the sources checked. The parts are
  classical: Taylor's theorem, here certified exactly; the divisibility theory of a Euclidean ring; a unit descent;
  and the condition at 2 of the 2-descent. No novelty is claimed.
- Not accessed: Taylor's texts; P. Lezowski, Math. Comp. 83 (2014), and his online tables (the server refused
  connections on 2026-10-26); Lemmermeyer's updated survey (both hosts refused connections); Godwin 1967; Angell 1973;
  MathSciNet.
