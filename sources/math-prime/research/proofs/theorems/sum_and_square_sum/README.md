# Integers with a prescribed sum and sum of squares: Σyᵢ=s, Σyᵢ²=N
Status: WRITTEN-AUDITED (one external audit, below)     Methods: not yet classified (four-square parity splitting; its four-square input is Lagrange's multiplier descent, as in [descent.classical](../../methods/descent.classical.md))
Generalizes: none     Generalized by: none

## Statement
For d ≥ 1 and an integer s, write s = da + r with 0 ≤ r < d and put M_d(s) = da² + 2ar + r. This is the least
value of Σyᵢ² over integer d-tuples with Σyᵢ = s; it is attained by r entries a+1 and d−r entries a.
1. **Theorem 1.** Let d ≥ 8 and let s, N be integers. There are integers y₁,…,y_d with Σyᵢ = s and Σyᵢ² = N
   if and only if N ≥ M_d(s) and N ≡ s (mod 2).
2. **Theorem 5 (box threshold).** Let d ≥ 8, let L ≤ U be integers with dL ≤ s ≤ dU, and put
   B = min{(L−1)² + M_{d−1}(s−L+1), (U+1)² + M_{d−1}(s−U−1)}. If N ≡ s (mod 2) and M_d(s) ≤ N < B, a
   representation exists and every representation has all yᵢ in [L, U]. B is the least squared norm of an integer
   vector of sum s that leaves the box, so the threshold is sharp for "every representation". Nothing is claimed
   about boxed existence at or above B.
3. **Corollary 6.** If s ≥ 0, d ≥ 8, N ≡ s (mod 2) and M_d(s) ≤ N < 1 + M_{d−1}(s+1), there is a representation
   with every yᵢ ≥ 0. With T_y = y(y+1)/2 it gives ΣT_{yᵢ} = (N+s)/2 together with Σyᵢ = s.

In Hoffmann's notation (below), S_m(n) is the set of sums T of m integers whose squares add up to n, and T*_m(n) is its
maximum. Theorem 1 says that for m ≥ 8, S_m(n) = {T ≡ n (mod 2) : M_m(T) ≤ n}. So T*_m(n) is the largest
T ≡ n (mod 2) with M_m(T) ≤ n.

## Proof 1: four-square parity splitting (branch `cloud-agent/elementary-connections`, round 2)
- **Files:**
  - [source/proof.md](source/proof.md): Theorem 1, Lemmas 2–4, Theorem 5, Corollary 6, and Appendix A, the
    four-square descent;
  - [source/PROVENANCE.md](source/PROVENANCE.md);
  - [source/experiment.py](source/experiment.py) with [source/results.json](source/results.json): 59,409
    constructions for d = 8…16, 6,190 box cases, and 18,343 dynamic-programming comparisons at d = 8;
  - [source/freeze.json](source/freeze.json): the hashes at submission.
- **Steps:**
  1. *Minimum.* Sum the inequality (z−a)(z−a−1) ≥ 0, which holds for every integer z.
  2. *Parity in four squares.* For j ∈ {1,2,3}, a representation of 4n+j as four squares has exactly j odd
     entries. Writing these as 2aᵢ+1 and the even entries as 2bᵢ gives n = Σ_{i≤j} aᵢ(aᵢ+1) + Σ bᵢ².
  3. *Affine pairs.* A pronic term a(a+1) becomes the pair (a+1, −a), with sum 1 and squared norm 2a(a+1)+1.
     A square b² becomes the pair (b, −b), with sum 0 and squared norm 2b². For residual sum 4, one square pair is
     recentred as (c+1, 1−c). So eight coordinates realize every (r, R) with 0 ≤ r ≤ 4, R ≥ r and R ≡ r (mod 2).
  4. *Translation.* Shift by an integer a nearest to s/d. When r = |s−da| > 4, add r−4 coordinates equal to 1;
     the total r+4 is at most d because d ≥ 8. Pad with zeros, and change sign if s−da < 0.
  5. *Box.* h(v) = v² + M_{d−1}(s−v) has difference 2(v − ⌊(s−v−1)/(d−1)⌋). This is negative below L and
     positive above U.
- **Lean:** none.
- **Audit:**
  - [source/review/review.md](source/review/review.md), with [check.py](source/review/check.py) and
    [results.json](source/review/results.json): a Codex reviewer distinct from the author. PASS for every statement.
    Its own exhaustive dynamic programming made 43,441 comparisons (d = 8…16, N ≤ 100), and it ran 159,805
    constructions (d = 8…32).
  - Consolidation check, 2026-10-26: the proof was read in full.
    [checks/check_affine.py](checks/check_affine.py) (output [checks/check_affine.json](checks/check_affine.json)) is
    an exhaustive bitset dynamic programme, independent of both branch scripts. It compared the criterion with the
    truth for every d ≤ 30, every N ≤ 220 and every sum:
    - no mismatch for 5 ≤ d ≤ 30;
    - mismatches for d ≤ 4, as expected (for example d = 4, s = 0, N = 28);
    - Hoffmann's Example 4.14 is reproduced: S₁₂(25) = {odd T : |T| ≤ 15}.
  - This is one external audit. The result has not had the two reviews that AGENTS.md asks of a lane draft.
- **Provenance:** the branch claims no novelty. The ingredients are classical:
  - Lagrange's four-square theorem (Appendix A is the classical descent);
  - parity counts in four-square representations;
  - Cauchy's lemma, which prescribes both moments for four nonnegative integers using Gauss's three-triangular
    theorem (Nathanson, Proc. AMS 99 (1987), as inspected by the branch's source audit).
- **Novelty** (web searches 2026-10-26, recorded in the consolidation notes):
  - *Where the question is.* D. W. Hoffmann, *Sums of integers and sums of their squares*, Acta Arith. 194 (2020)
    295–313, doi:10.4064/aa190219-15-10, arXiv:1902.07109. Hoffmann proves:
    - S_m(n) is a parity interval for m ≥ 5 (Prop. 4.1);
    - S_m(n) is "full" for 5 ≤ m ≤ 8 (Thm 3.4(iv) and Thm 4.7(i));
    - S′_m(n) = {T : T² ≤ m(n−1)+1, T ≡ n (mod 2)} for 9 ≤ m ≤ 11 and n > m (Thm 4.7(ii)).
    His Example 4.14 shows that this bound fails for m = 12. As **Problem 4.15** he leaves open an explicit
    description of T*_m(n).
  - *What Theorem 1 adds.* It answers Problem 4.15 for every m ≥ 8.
    - For 8 ≤ m ≤ 11 it agrees with Hoffmann's Thm 4.7 and Cor. 4.8, by a different, elementary proof. This was
      checked by hand: for m = 9 the criterion makes S₉(n) non-full exactly when 9n−2 is a square, as in his
      Cor. 4.8(i), and Examples 4.9–4.11, 4.13 and 4.14 come out right.
    - For m ≥ 12 it is new in the sources found.
  - *Small m.* For 5 ≤ m ≤ 7 the same criterion follows from Hoffmann's fullness; his proof uses Mordell's theory of
    binary forms as sums of squares of linear forms, and was not inspected here. For m = 4 the criterion is false:
    T ∈ S₄(n) exactly when T ≡ n (mod 2) and 4n−T² is a sum of three squares (Goldmakher–Pollack, Amer. Math. Monthly
    125 (2018)).
  - *Later work.* Four searches found no later paper answering Problem 4.15. OpenAlex lists one paper citing
    Hoffmann, on weighted sums of polygonal numbers (Acta Arith. 2022, doi:10.4064/aa200630-1-6), which does not.
  - *Confidence:* medium. A MathSciNet or zbMATH citation search was not available.
- **Related, left out:**
  - *The branch's first-campaign return theorem* (`cloud_agent/convex_geometry/proof.md` on branch
    `cloud-agent/elementary-connections`, not imported): a rational vector with integer sum s and squared norm
    N ≡ s (mod 2) returns to an integer one when d ≤ 7.
    - It is a special case of F. Dacar's affine Aubry–Davenport–Cassels theorem (arXiv:1310.2093, Theorem 1),
      applied to half the norm form of A_{d−1}. The squared covering radius ⌊d²/4⌋/d of A_{d−1} is below 2
      exactly when d ≤ 7.
    - For 5 ≤ d ≤ 7 it also follows from Hoffmann's fullness. A rational point forces N ≥ s²/d, and with the parity
      this gives N ≥ M_d(s), because M_d(s) − s²/d = r(d−r)/d < 2 for d ≤ 7.
  - *Nonnegative representations in four coordinates with both moments prescribed:* Cauchy's lemma.

## Open questions
- Representation counts with both moments prescribed.
- Existence inside a box at or above the threshold B (the branch's NEXT_TARGETS §1).
- d = 3, the Eureka shape: n = T_a + T_b + T_c with a+b+c = s is the case d = 3, N = 2n−s. Hoffmann's Thm 3.4(ii)
  decides it through Mordell's results, which are not elementary. Nothing here reaches d < 5.
