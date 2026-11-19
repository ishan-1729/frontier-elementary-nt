# Literature check: y² = x³ + 17 (10404.f2), and R1 on it

Written 2026-10-15 by a Claude subagent for the R1 session, and edited by the author. Section 1 lists the published
determinations of the
integral points, section 2 compares their mechanisms with R1, section 3 gives the provenance of hypothesis (H),
section 4 the verdict, section 5 corrections to project files, and section 6 what was not read. Section 7 (added
2026-10-26) covers the arithmetic of the ring behind [hypothesis_H.md](hypothesis_H.md). Each source is marked
with how it was read: [full], [pages], [review], [abstract], [snippet] or [not accessed]. "Not found" means not found
in the sources listed here. The computations quoted in sections 2 and 3 are checks made for this note, not
literature: [lit_checks_17.py](scripts/lit_checks_17.py) → [lit_checks_17.out.txt](scripts/lit_checks_17.out.txt)
(Appendix A).

**Conventions.** One curve, many names:
- y² = x³ + k with k = 17: LMFDB, Gebel–Pethő–Zimmer, Bennett–Ghadermarzi;
- y² − k = x³ with k = 17: Mordell, Nagell, Hemer, London–Finkelstein, Ljunggren;
- Y² = X³ − D with D = −17: Cassels, Selmer;
- x² − 17 = y³: Cunningham; x² + D = yⁿ with D = −17, n = 3: Barros;
- x² + k = y³ with k = −17: this project's catalogue (hence the slug mordell_k_minus17).

The curve y² = x³ − 17 (Gerono 1877; no integral points) is a different curve.

## 1. Published determinations of the integral points

### 1.1 Complete determinations
1. **T. Nagell (1930)** [review]. *Einige Gleichungen von der Form ay² + by + c = dx³*, Avh. Norske Vid.-Akad. Oslo I,
   Mat.-Naturv. Kl. 1930, no. 7, 15 pp. JFM 56.0877.02, review by A. Brauer, read through the zbMATH Open API.
   - The review says that no general method is known for ay² + by + c = dx³, and that the paper solves
     y² − 9 = x³, y² − 17 = x³ and 3y² − 2 = x³.
   - The last one extends Nagell's 1925 theorem on x³ + Dy³ = 1: at most one solution with y ≠ 0, and then
     x₁ + y₁∛D is the fundamental unit, or its square for D = 19, 20, 28.
   - The review says nothing about how y² − 17 = x³ is solved.
   - This is the first complete determination found. OEIS A029728 cites it, and Hilliker (1982, below) credits the
     complete solution to Nagell, by algebraic number theory.
   - The National Library of Norway has digitized it (URN:NBN:no-nb_digibok_2016052308051), readable only from
     Norway. Not read.
   - Its method is known here only by inference (2.1).
2. **O. Hemer (1952, 1954).**
   - *On the Diophantine equation y² − k = x³*, Dissertation, Uppsala 1952, 101 pp., Zbl 0049.31004, review by
     B. Stolt [review].
     - Part 1 makes y² − kf² = x³ equivalent to binary cubic forms, by comparing irrational parts in
       ±y + f√k = ε(a + b√k)³, with ε a unit of Q(√k). The forms have negative discriminant when k > 0.
     - For k > 0 the results follow from Nagell's and Delone's theory of cubic fields with exactly one fundamental
       unit.
     - Part 2 solves 0 < kf² ≤ 100 numerically; some conclusions were wrong and were corrected in 1954.
   - *Notes on the Diophantine equation y² − k = x³*, Ark. Mat. 3 (1954) 67–77, DOI 10.1007/BF02589282 [snippet].
     - Project Euclid asked for a human verification and Springer is paywalled. The zbMATH record has no review.
     - A search engine's summary of the Project Euclid PDF: corrections and additions to the thesis; all equations
       with 0 < k ≤ 100 are solved; the start is the reduction to F(u, v) = 1 with F a binary cubic form.
     - So it covers k = 17. Whether it re-proves k = 17 or quotes Nagell was not seen.
   - The method, read in Hemer, *On some diophantine equations of the type y² − f² = x³*, Math. Scand. 4 (1956)
     95–107, DOI 10.7146/math.scand.a-10458 [pages 95–97, 106–107]:
     - F(u, v) = 1, with F of negative discriminant, becomes u + vρ = μᵢεⁿ, with ε the fundamental unit of R(ρ). A
       solution is an exponent n at which one coordinate of μᵢεⁿ vanishes: one exponent, one marker.
     - The Delone–Nagell theorem allows at most three solutions, except for discriminants −23, −31 and −44 (his
       Lemma 2).
     - Single exponents are excluded by a finite congruence (his Lemma 3). Let p be an odd prime and α = aρ² + bρ + c
       an integer of R(ρ) prime to p, with a ≡ b ≡ 0 (mod pᵏ) and ta + sb ≢ 0 (mod p²ᵏ). Write αⁿ = aₙρ² + bₙρ + cₙ.
       Then modulo p^(2k+h), where pʰ ∥ n, aₙ ≡ ncⁿ⁻¹a and bₙ ≡ ncⁿ⁻¹b, so taₙ + sbₙ ≠ 0 for n ≠ 0. It
       generalizes Lemma 8 of his 1954 note; his next lemma generalizes one of Delone.
     - That lemma is completion-free and has one exponent. Proposition 17 of
       [lifting.binomial.multi](../../methods/lifting.binomial.multi/proof.md) is a lemma of the same kind in two
       exponents.
3. **H. London, R. Finkelstein**, *On Mordell's equation y² − k = x³*, Bowling Green State Univ., 1973, 146 pp.,
   MR 340172 [not accessed]. Gebel–Pethő–Zimmer (1998) cite it for the complete solution of |k| ≤ 100. Its method was
   not seen. (Gebel–Pethő–Zimmer print the authors as "J. London, M. Finkelstein"; zbMATH and Hilliker give Hymie
   London and Raphael Finkelstein.)
4. **J. Gebel, A. Pethő, H. G. Zimmer**, *On Mordell's equation*, Compositio Math. 110 (1998) 335–367, MR 1602064
   [full: the authors' re-typeset copy dated 2009, on Pethő's page; the journal's page numbers were not checked].
   - They solve all |k| ≤ 10⁴, and most |k| ≤ 10⁵, by the Lang–Zagier method:
     - a Mordell–Weil basis;
     - David's lower bound for linear forms in elliptic logarithms;
     - de Weger's reduction.
   - They state that their ranks are unconditional for |k| ≤ 10⁴: Cremona's 2-descent, and a 3-descent for
     k = −7954 and 8206.
   - k = 17 is in their Table 9 (Hall's conjecture): x = 5234, √x/|k| ≈ 4.26.
5. **M. A. Bennett, A. Ghadermarzi**, *Mordell's equation: a classical approach*, LMS J. Comput. Math. 18 (2015)
   633–646, DOI 10.1112/S1461157015000182, arXiv:1311.7077 [full].
   - They solve all 0 < |k| ≤ 10⁷.
   - Their Proposition 2.1 is Mordell's correspondence. Integral points correspond bijectively to SL₂(Z)-classes of
     triples (F, x, y), with F = ax³ + 3bx²y + 3cxy² + dy³ of discriminant −108k and F(x, y) = 1.
   - The cubic Thue equations are solved in Magma by the method of Tzanakis–de Weger: Baker's lower bounds and
     lattice reduction.
   - Their introduction sets out both classical routes:
     - Q(√k) leads to cubic Thue equations, and Q(∛k) to quartic ones.
     - For k > 0 the cubic fields have one fundamental unit, and Skolem's p-adic method treats them rather easily.
     - For k < 0 one meets cubic or quartic fields with two fundamental units, citing Finkelstein–London 1970,
       Hemer 1954 and Ljunggren 1963.
   - k = 17 is in their Table 5 (Hall's conjecture): X = 5234, ratio 4.26.
6. **Software and databases.**
   - K. Wildanger, *Über das Lösen von Einheiten- und Indexformgleichungen in algebraischen Zahlkörpern mit einer
     Anwendung auf die Bestimmung aller ganzen Punkte einer Mordellschen Kurve*, Dissertation, TU Berlin 1997
     [not accessed; only its title and a 2010 seminar abstract by Pohst were seen].
   - The KANT/KASH manual entry IntegralPoints [full] uses y² = x³ + 17 as its example.
     - IntegralPoints(17) returns the 16 points, by the algorithm of a reference keyed "Wi3" (presumably
       Wildanger; not confirmed).
     - Its other form, IntegralPoints(a, b), reduces y² = x³ + ax + b to quartic Thue equations (keys "Stroe1",
       "TzWe1"). That is a quartic route solved by Baker-type methods. The example does not use it.
   - Magma's IntegralPoints and Sage's integral_points (programs in OEIS A029728), and LMFDB 10404.f2 [full]:
     - 16 points, rank 2, generators (−2, 3) and (4, 9);
     - the LMFDB page states no method and no proof status for the rank.
   - OEIS A029728 and A029727 [full].
     - They give the x and y values, with references to Mordell (1969, p. 246), Nagell (1930) and Silverman–Tate
       (1992).
     - A comment by H. Cohen describes a completeness proof by elliptic logarithms (Tzanakis–de Weger, as in Smart's
       book). LLL reduces the coefficient bound to about 200, then 20, then 12.
   - OEIS A081119 [full] gives 16 solutions at n = 17, citing Gebel–Pethő–Zimmer and Bennett–Ghadermarzi.
7. **Families that contain the theorem; their methods were not seen.**
   - C. F. Barros, *On the Lebesgue–Nagell equation and related subjects*, PhD thesis, Warwick 2010 [abstract].
     - It treats x² + D = yⁿ (n ≥ 3) for −100 ≤ D ≤ −1, which includes x² − 17 = y³.
     - How n = 3 is handled was not seen; the repository has no full text.
   - R. von Känel, B. Matschke, arXiv:1605.06079 [abstract]. It gives algorithms for Mordell equations through the
     Shimura–Taniyama conjecture, without linear forms in logarithms. Whether k = 17 is in their tables was not
     checked.

### 1.2 Statements and searches, not proofs
- **A. Cunningham (1905)**, solution of Question 15697, Math. Questions and Solutions from the Educational Times (2) 8,
  53–54 [not accessed; known through Gauthier–Lê].
  - The question asked for all integral solutions of x² − 17 = y³. Cunningham answered by a search, eliminating
    values by last digits and congruences.
  - Gauthier–Lê, who write the equation as y² − 17 = x³, say he listed the five solutions with x < 500 (y > 0).
    Exactly five have 0 < x < 500 (x = 2, 4, 8, 43, 52), so presumably x > 0 was meant; that reading is ours.
  - This is the earliest appearance of the problem found.
- **S. Gauthier, F. Lê**, *On the youthful writings of Louis J. Mordell on the Diophantine equation y² − k = x³*, Arch.
  Hist. Exact Sci. 73 (2019) 427–468 [full, authors' postprint]. A historical study. It is the source here for
  Cunningham, for Gerono, and for the remark (their §4, p. 33 of the postprint) that the integral-solution papers
  cited in Mordell's 1969 chapter start from Mordell's reduction to binary cubic forms, Nagell 1930, Hemer 1954 and
  Ljunggren 1963 among them.
- **L. J. Mordell**, *Diophantine Equations*, Academic Press 1969, Ch. 26, pp. 238–254 [abstract: the ScienceDirect
  chapter page, with its summary and part of its reference list].
  - The statement of the 16 points is at p. 246 according to OEIS, and at p. 250 according to Brown–Myers. Neither
    page was seen.
  - The references shown include Nagell 1930, Hemer 1952/1954, Ljunggren 1942 and 1963, Lal–Jones–Blundon and
    Fueter 1930.
- **D. Zagier**, *Large integral points on elliptic curves*, Math. Comp. 48 (1987) 425–436 [full; a scan with a text
  layer, on Zagier's MPIM page].
  - pp. 430–431 list the eight points with y > 0. Their coordinates in the basis (−2, 3), (−1, 4) are read off from
    elliptic logarithms.
  - There is no completeness proof. On p. 431 Zagier notes that Baker's bound is far too large, and points to
    elliptic-logarithm bounds of Masser–Wüstholz and Wüstholz.
  - Table 1 (p. 433), entry (c), lists 5234 among the large points of Odlyzko's search.
- **D. L. Hilliker**, *An algorithm for solving a certain class of Diophantine equations I*, Math. Comp. 38 (1982)
  611–626 [abstract and reference list on the AMS page; snippet of the PDF from a search engine].
  - It uses Mordell's equation as motivation and credits the complete solution of y² − 17 = x³ to Nagell (1930),
    by algebraic number theory.
  - Its own method computes bounds for a class of equations, and its references point to Runge's theorem. No
    re-proof of k = 17 was seen.
- **E. Brown, B. T. Myers**, *Elliptic curves from Mordell to Diophantus and back*, Amer. Math. Monthly 109 (2002)
  639–649 [full]. It quotes Mordell's book (p. 250) for exactly sixteen integer points, and generates them from
  (−1, 4) and (−2, 3). No proof.
- **J. H. Silverman, J. Tate**, *Rational Points on Elliptic Curves* (1992) [not accessed]; cited by OEIS A029728.

### 1.3 Negative checks (other curves)
- Gerono 1877, through Gauthier–Lê: x³ = y² + 17, that is y² = x³ − 17.
- W. J. Ellison, F. Ellison, J. Pesek, C. E. Stahl, D. S. Stall, J. Number Theory 4 (1972) 107–117 [abstract]: Baker's
  method for y² = x³ − 28 only.
- R. Finkelstein, H. London, J. Number Theory 2 (1970) 310–321 [abstract]: y² + 18 = x³ only.
- J. Gebel, A. Pethő, H. G. Zimmer, Acta Arith. 68 (1994) 171–192 [full, authors' copy]: its example is Mestre's
  rank-6 curve.
- K. Conrad, *Examples of Mordell's equation* (expository note) [full]: does not treat k = 17.

## 2. Mechanisms, compared with R1

### 2.1 The cubic route (Nagell 1930 by inference; Hemer; Bennett–Ghadermarzi)
- **Reduction.** Mordell's correspondence with binary cubic forms of discriminant −108·17 = −1836 that represent 1,
  or, in Hemer's form, the factorization ±y + √17 = ε(a + b√17)³ in Q(√17).
- **Fields.** The cubic fields have negative discriminant, signature (1,1) and one fundamental unit. Each Thue equation
  is a problem in one exponent, u + vρ = μεⁿ.
- **Exclusion.** Classically by the Delone–Nagell theory and finite congruences (Hemer's Lemma 3 above). Skolem's
  p-adic method also applies, as Bennett–Ghadermarzi remark; they themselves use Baker's bounds and lattice reduction.
- **The four equations for k = 17** (Appendix A, part A):
  - The 16 points fall into four GL₂(Z)-classes of such forms. They carry X ∈ {−2, 8, 5234}, {−1, 43}, {2, 52}
    and {4}.
  - The first class contains x³ + 6xy² + 6y³. It represents 1 at (1, 0), (1, −1) and (−23, 26), giving X = −2, 8
    and 5234: three representations, the most that the Delone–Nagell theorem allows at this discriminant.
  - The four fields have discriminants −1836, −459, −1836 and −204, all −51 times a square. They are pairwise
    non-isomorphic, and none is K = Q(θ), θ³ = −17 (discriminant −867).
  - No other representation of 1 has |v| ≤ 1000. This is a search, not a proof.
- **Why Nagell 1930 is placed here (inference).**
  - Gauthier–Lê's remark on the introductions (1.2).
  - Hemer's thesis review: for k > 0, the one-unit theory of Nagell and Delone.
  - Hilliker: algebraic number theory.
  - Dates: the two-unit Skolem arguments recorded in this project begin with Skolem in 1933 (x⁴ − 2y⁴ = ±1), after
    Nagell's paper.
  - Nagell's own normalization of the forms was not seen.

### 2.2 Elliptic logarithms (Gebel–Pethő–Zimmer 1998; Magma, Sage; Cohen's comment in OEIS)
A Mordell–Weil basis with generators (−2, 3) and (4, 9), David's lower bound for linear forms in elliptic logarithms,
and LLL reduction of the coefficient bound, then a search over small coefficients.

### 2.3 Cubic Thue equations with Baker's method (Bennett–Ghadermarzi 2015)
The route of 2.1, with every Thue equation solved by the method of Tzanakis–de Weger: lower bounds for linear forms in
complex logarithms, then lattice reduction.

### 2.4 R1 on this curve: the quartic route
- **Descent.** The 2-descent over K = Q(θ), the field Cassels used for the rank (section 3). It has three non-trivial
  classes, which give eight quartic Thue equations G = c.
- **Fields.** Three quartic S4 fields of signature (2,1), with discriminants −499392, −7803 and −55488
  ([00_fields.out.txt](scripts/00_fields.out.txt)). Each has cubic resolvent field K (Appendix A, part B). Each has
  two units, so every equation has two exponents.
- **Coverage** by the cone certificate of Lemma 3.2. **Exclusion** by Proposition 17 and a congruence sieve.
- **No field in common.** The four fields of the cubic route are S3 fields whose quadratic resolvent is Q(√−51),
  since a cubic field's discriminant is the resolvent's discriminant times a square. R1's fields are S4 fields with
  cubic resolvent K. The two routes share no Thue field and no equation.
- **Quartic route for k = 17 in print: not found.**
  - Bennett–Ghadermarzi name the Q(∛k) route only in general.
  - KANT implements a quartic route for general Weierstrass curves with Baker-type solving. That is software, not a
    published proof for k = 17.
  - A two-index Skolem argument on this curve was not found either.
- **Precedents of R1's ingredients** are as in
  [triangular_tetrahedral/literature_check.md](../triangular_tetrahedral/literature_check.md) §2:
  - two-index Skolem arguments, from Skolem (1933–34) to Avanesov (1966/67) on 1611.a1;
  - Uchiyama's completion-free congruences (1984) on 8604.a1 (his printed argument has gaps:
    [triangular_square_pyramidal/literature_check.md](../triangular_square_pyramidal/literature_check.md) §1.7);
  - Koecher's reduction theory for the coverage lemma.
- **One precedent specific to Mordell's equation.** The one-exponent finite congruence of Hemer (1956, Lemma 3,
  generalizing Lemma 8 of his 1954 note) is the one-unit analogue of Proposition 17.

## 3. Hypothesis (H) and the rank
- **Cassels's theorem is (H).** J. W. S. Cassels, *The rational solutions of the diophantine equation Y² = X³ − D*,
  Acta Math. 82 (1950) 243–273, DOI 10.1007/BF02398279 [pages 243–244, 249–250, 261–262 and 268–273; a scan with a text
  layer, page images checked; the scan was read in the session's scratch folder and is not stored in the
  repository].
  - §22 covers D ≡ ±1 (mod 9), cube-free, with odd class number h_D of Q(∛D). There the argument of §12 leaves only
    the classes μ = 1, η, λ, λη for x − δt², with δ³ = D. The classes λ and λη occur exactly when 3 | y.
  - For D = −17 (δ = θ), Cassels finds both η and η₁λ realized, by (−1)³ − 4² = (−2)³ − 3² = −17. So the group has
    two generators of infinite order (p. 262), and Table 1 (p. 269) lists (−1, 4) and (−2, 3).
  - Table 2 (p. 270) gives, for the field Q(∛17), h = 1 and η = 18 − 7∛17, which is 18 + 7θ.
  - So, for every rational point, x − θ lies in one of four classes. That is (H).
- **The classes match** (Appendix A, part C).
  - R1's class (1,1), μ = −1 − θ, is Cassels's η: (−1 − θ)η is a square in K.
  - Classes (1,0) and (0,1) are λ and λη, in an order not fixed here; the product of the three μ is a square.
  - This agrees with Cassels's criterion: among the known integral points, 3 divides y exactly in classes (1,0) and
    (0,1).
- **Its inputs.** §§11–12 use ideal factorization in K, the class number (h odd) and the fundamental unit η. Cassels's
  Lemma 6 (if h is odd, then η is not a quadratic residue of 4) rests on class field theory, or on Hilbert's theory
  of relative quadratic fields, but it is not used for the list at D = −17: §22 takes the list "by the argument of
  § 12", and §12 (pp. 250–251) uses only the ideal factorization (his Lemma 7), the class group and the units. In the
  paper Lemma 6 serves other D, for example D = +17 through Theorem VIII, Corollary 1 (p. 262). (Corrected in round 2
  of the review; the first version of this note listed Lemma 6 among the inputs.)
- **Other sources for the rank.**
  - G. Billing, Nova Acta Reg. Soc. Sci. Upsal. (4) 11 (1938) no. 1 [not accessed]: Cassels says it gives a
    complete basis for |D| < 25 and agrees with his results except at D = −15.
  - E. S. Selmer, Math. Scand. 4 (1956) 281–286 [review]: tables for |D| ≤ 100.
  - Gebel–Pethő–Zimmer (1998): unconditional ranks for |k| ≤ 10⁴.
  - A corrigendum to Cassels, Acta Math. 84 (1951), was not checked.
- **Consequence.** Theorem B holds unconditionally once Cassels's §22 is admitted. That input is ideal-theoretic, not
  elementary in the project's sense. "LMFDB does not call the rank proven" is true of the LMFDB page, but (H) itself is
  a published theorem.

## 4. Verdict
- **Theorem.** Classical.
  - First complete solution: Nagell (1930).
  - Tables: Hemer (1952/1954, 0 < k ≤ 100) and London–Finkelstein (1973, |k| ≤ 100).
  - Computer proofs: Gebel–Pethő–Zimmer (1998, elliptic logarithms) and Bennett–Ghadermarzi (2015, cubic Thue
    equations with Baker's method). Wildanger's algorithm (1997, not read) is KANT's IntegralPoints(k), whose manual
    uses k = 17 as its example.
- **Mechanisms in print.**
  - The algebraic proofs (Hemer by its review, Nagell by inference) take the cubic route: four cubic Thue equations,
    in four fields with one unit each, one exponent each.
  - The computer proofs use elliptic logarithms or Baker's bounds.
- **R1 here.** The quartic route:
  - the 2-descent over K;
  - eight quartic Thue equations in three S4 fields of signature (2,1) with cubic resolvent K;
  - two units, so two exponents;
  - a finite Jacobian criterion with a sieve.

  No published treatment of this curve by this route, or by any two-exponent Skolem-type (unit) argument, was
  found.
- **Novelty, by the project's criterion.** The same proof, or a rewording of it, was not found for:
  - this theorem;
  - x² − 17 = y³;
  - any family of Mordell equations containing k = 17.

  So R1 gives a proof of this theorem that is not in the sources checked. This is subject to the unread primary
  sources of section 6, above all Nagell 1930, Hemer 1952/1954, London–Finkelstein 1973 and Ljunggren 1963.
- **Owner's call.** The cubic route and R1 have the same overall shape: descent, then Thue equations, then units, then
  congruences in the exponents. They share no field and no equation, and they differ in the number of units and of
  exponents. Whether that shared shape makes R1 a rewording is Ishan's decision.
- **What is not new.**
  - The two-index mechanism: Avanesov on 1611.a1, and Uchiyama completion-free on 8604.a1 (with gaps in print).
  - The one-exponent congruence lemma for Mordell's equation: Delone; Hemer 1954, 1956.
  - The descent and its four classes: Cassels 1950.

  The claim must not be "a new method". It can be "the first proof of this theorem by the quartic, two-unit route,
  in the sources checked".
- **Elementary proof.** None found. Every published proof uses unique factorization or class numbers, together with
  unit groups of quadratic or cubic fields, or else linear forms in logarithms. Theorem A of R1 is elementary.
  Theorem B needs (H), which is known only through Cassels's ideal-theoretic proof.

## 5. Corrections to project files
- **[ghost_screen_2026-10-14.md](../../targets/ghost_screen_2026-10-14.md) line 34 (row 4).** Its "Existing proofs"
  entry, "Elliptic logarithms plus lattice reduction", is incomplete.
  - The first proof is Nagell 1930, by algebraic number theory. The screen's own source, OEIS A029728 (line 62),
    cites it.
  - Later proofs: Hemer 1952/1954 (cubic forms, one-unit cubic fields), London–Finkelstein 1973, Gebel–Pethő–Zimmer
    1998 (elliptic logarithms) and Bennett–Ghadermarzi 2015 (cubic Thue equations with Baker's method).
  - The rest of the row agrees with the literature: discriminant −108·17, one unit per cubic Thue equation, Nagell's
    16 points.
  - *Addressed (2026-10-15).*
- **[literature_screen_2026-10-15.md](../../targets/literature_screen_2026-10-15.md) line 33.** Its "Known proofs"
  entry, "Elliptic logarithms with lattice reduction", has the same omission. *Addressed (2026-10-15).*
- **[README.md](README.md) and [proof.md](proof.md).** The drafts said only that (H) rests on an mwrank value that
  LMFDB does not call proven. (H) is a published theorem: Cassels 1950, §22 (section 3). It is ideal-theoretic, so for
  this card it stays an input that is not proved elementarily. *Addressed:* the card and proof.md §6.1 now cite
  Cassels.
- **[rank2_signatures_2026-10.md](../../targets/rank2_signatures_2026-10.md) lines 124–126.** It says that for the
  comparison curves the rank is LMFDB's (mwrank). For 10404.f2 the rank and the class list are classical (Cassels
  1950; Billing 1938). This is an addition, not an error. *Addressed (2026-10-15).*
- **No error was found** in the other lines that mention this curve:
  [rank2_signatures_2026-10/literature.md](../../targets/rank2_signatures_2026-10/literature.md) lines 27, 438 and
  458,
  [GENERALIZATIONS.md](../../GENERALIZATIONS.md), [ghost.unit.md](../../methods/ghost.unit.md) and
  [lifting.binomial.multi.md](../../methods/lifting.binomial.multi.md).
- **Outside the project.**
  - Bennett–Ghadermarzi's reference [28] gives Ljunggren's paper as "The diophantine equation y² = x³ − k", Acta
    Arith. 8 (1961) 451–465. The IMPAN and zbMATH records give *On the diophantine equation y² − k = x³*, Acta Arith.
    8 (1963) 451–463.
  - The page of Mordell's statement is given as 246 by OEIS and as 250 by Brown–Myers; this is unresolved.

## 6. Not read
- **Nagell 1930**: the full text (readable only from Norway); only the JFM review was read.
- **Hemer**: the 1952 thesis (not online) and the 1954 note (Project Euclid needed a human check; Springer is
  paywalled).
- **London–Finkelstein 1973**; Ljunggren, Acta Arith. 8 (1963) 451–463 (the matwbn server refused the fetch tool).
- **Mordell 1969**, Ch. 26 (pp. 238–254), including p. 246 or p. 250.
- **Wildanger 1997**; Barros 2010; von Känel–Matschke beyond the abstract.
- **Books**: Silverman–Tate 1992; Smart's *The Algorithmic Resolution of Diophantine Equations* (1998); Cohen,
  *Number Theory I* (2007); Tzanakis, *Elliptic Diophantine Equations* (2013); Skolem, *Diophantische Gleichungen*
  (1938).
- **Papers**: Stroeker–Tzanakis, Acta Arith. 67 (1994); Lal–Jones–Blundon, Math. Comp. 20 (1966); Hall 1971; Billing
  1938; Cassels's 1951 corrigendum; Selmer beyond the review; Delone 1930 and Nagell 1925 and 1928 (the one-unit
  theory); Cunningham 1905 and Gerono 1877 in the original.
- **Not searched**: MathSciNet; the Educational Times volumes; theses beyond those named.

## 7. Hypothesis (H): the arithmetic of Z + Zθ + Zω, θ³ = −17 (added 2026-10-26)
For [hypothesis_H.md](hypothesis_H.md). Same check, marks and conventions as
[triangular_square_pyramidal/literature_check.md](../triangular_square_pyramidal/literature_check.md) §6 (the
literature agent of review round 3, edited by the author). "Exact" marks an identification checked in sympy by
[lit_H_lemmermeyer.py](scripts/lit_H_lemmermeyer.py) or [lit_H_lmfdb.py](scripts/lit_H_lmfdb.py).

**7.1 Class number and unit.**
- LMFDB 3.1.867.1 [full, 2026-10-26]: defining polynomial x³ − x² + 6x − 12, whose root is ω (exact); class number 1;
  fundamental unit (7/2)a² − (7/2)a − 4, which is −η (exact). The page does not say whether its class number is
  proven.
- E. S. Selmer, *Tables for the purely cubic field K(∛m)*, Avh. Norske Vid.-Akad. Oslo I 1955, no. 5 (1956),
  Zbl 0070.03603 [review]: class numbers and fundamental units for the 84 cube-free m ≤ 100.
- Cassels 1950, Table 2 (section 3): h = 1 and η. The round-3 referee and the author checked the page: the entry for
  D = 17 is 18 − 7δ with δ = ∛17, without the star that marks units not known to be fundamental (note 2, p. 271).
- Lemmermeyer 1989, arXiv:2608.23216 [pages 91, 101, 149]: units of Q(∛m), m ≤ 70, go back to Markov (1892);
  Nakamula's table (1988) gives h = 1 for m = 17; row 107 lists −867 with x³ + 7x² + 5x + 2, h = 1 and the unit η⁻¹
  (exact).
- hypothesis_H.md needs none of these: Sections 2–3 give the divisibility theory, and (4.1) proves that η is
  fundamental.

**7.2 Not norm-Euclidean; the minimum.**
- V. G. Cioffari, *The Euclidean condition in pure cubic and complex quartic fields*, Math. Comp. 33 (1979) 389–398,
  Zbl 0399.12001 [abstract, through Crossref; the AMS page returned 403]: Q(∛d) is norm-Euclidean exactly for
  d = 2, 3, 10. This is the first exclusion of m = 17 found. Whether Taylor's lists of 1975–76 (|d| < 5000) contain
  −867 was not seen.
- Lemmermeyer 1989 [pages 32, 91–94] follows Cioffari, with sharper bounds. For m = 17 he gives the unit 18 − 7∛17,
  which is η, and the point x = (−376 + 466∛17 − 19α)/1028 with α = (1 − ∛17 + ∛17²)/3 = ω, where M(K, x) = 1115/1028.
  Exactly, x = −ξ* of Remark 2.5.
- Cavallar–Lemmermeyer 1998, arXiv:1202.6019v1 [full], §7 table, p. 17: "−867 N 1115/1028", that is, class number 1,
  not norm-Euclidean, M(K) = 1115/1028. The round-3 referee read the same row.
- *Euclidean windows* (2000) [pages]: −867 is neither in the conjectured list of norm-Euclidean fields nor in Table 3,
  the fields with a weighted-norm window.
- **Verdict.** The recorded minimum equals the lower bound of Remark 2.5, at the same point up to sign;
  1028 = |N(η + 1)| fits a point fixed by −η (exact). Remark 2.5 re-proves Lemmermeyer's value by an exact finiteness
  argument. The upper bound M ≤ 1115/1028 and the second minimum are not proved in the card, and are not needed.

**7.3 Division with a multiplier.**
- Not found for this ring: a published Dedekind–Hasse certificate, with integer multipliers or others; a two-stage
  division (the complex cubic examples found are −199 and −351, Lemmermeyer 1989, p. 7, and −307 and −335, Lezowski
  2014, see [triangular_tetrahedral/literature_check.md](../triangular_tetrahedral/literature_check.md) §5.2); a
  weighted norm (*Euclidean windows*, Table 3, has −199, −327, −351 and −367).
- Under GRH, a ring of integers with infinitely many units is Euclidean for some function exactly when it is a
  principal ideal domain (Weinberger 1973, as reported in the zbMATH review of D. A. Clark, Math. Comp. 65 (1996),
  Zbl 0853.11084 [review]). That is conditional and gives no division.
- A finite Dedekind–Hasse procedure is in print for definite quaternion orders: A. Cardoso, A. Machiavelo, *The
  Dedekind–Hasse criterion in quaternion algebras*, arXiv:2506.22651 (v3, 2026) [abstract; pages 5–10]. It does not
  treat number fields.

**7.4 Earlier treatments of (H).**
- Cassels 1950, §§11–12, 22 (section 3) proves the list for all rational points, with ideal factorization, the class
  number and the tabulated η.
- E. S. Selmer, Math. Scand. 4 (1956) 281–286, Zbl 0077.05001 [review], works with elements whose principal ideals are
  squares of ideals. Billing 1938 [not accessed].
- **Verdict.** No earlier ideal-free proof of (H) was found, and no published division algorithm for this ring that
  could give one. The descent is Cassels's, written for elements; h = 1 and η are classical; the {1, 2} table and the
  unit descent are certificates of known facts. No novelty is claimed.
- Not accessed: the texts of Taylor and Cioffari; Lezowski 2014 and his tables; Lemmermeyer's updated survey; Markov
  1892; Nakamula 1988.

## Appendix A. Computations for this check
[lit_checks_17.py](scripts/lit_checks_17.py) → [lit_checks_17.out.txt](scripts/lit_checks_17.out.txt) (Python 3.13.5,
sympy 1.14.0; about 30 s; exits non-zero if a check fails). Exact arithmetic only: integer roots, the round-two maximal
order, field isomorphism, and factorization over K.
- **A.** The cubic route: for each integral point P = (X, Y), the form F_P = x³ − 3Xxy² + 2Yy³ of discriminant
  −108·17, all its representations of 1 with |v| ≤ 1000, and the points they carry. The 16 points fall into four
  GL₂(Z)-classes, carrying X ∈ {−2, 8, 5234}, {−1, 43}, {2, 52} and {4}; their cubic fields have discriminants −1836,
  −459, −1836 and −204, are pairwise non-isomorphic, and none is K (discriminant −867).
- **B.** The cubic resolvent field of each R1 class form is K.
- **C.** Cassels's classes: η = 18 + 7θ has norm 1 and is not a square in K; (−1 − θ)η is a square, so R1's class
  (1,1) is Cassels's η; the product of the three μ is a square.

The script was written by the literature subagent in a scratch folder and moved here unchanged by the author.
