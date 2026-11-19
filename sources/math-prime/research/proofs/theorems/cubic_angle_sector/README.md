# A cubic angle-division form: min |U³−3aU²V−3UV²+aV³| = 8a on the sector U ≥ 2V > 0
Status: WRITTEN-AUDITED (one external audit, below)     Methods: [approx.normsign.all_degree](../../methods/approx.normsign.all_degree.md)
Generalizes: none     Generalized by: none

## Statement
Let Φ_a(U,V) = U³ − 3aU²V − 3UV² + aV³ = Re((U+iV)³) − a·Im((U+iV)³). The roots of Φ_a(1,t) are
tan((arctan(1/a) + kπ)/3).
1. For every integer a ≥ 1000 and all integers U > 3aV > 0: |Φ_a(U,V)| > 6a²/125³.
2. For every a ≥ 1000 there are no integers U ≥ 2V > 0 with |Φ_a(U,V)| ≤ 3.
3. For every a ≥ 2,604,167, the minimum of |Φ_a(U,V)| over integers U ≥ 2V > 0 is 8a, attained only at
   (U,V) = (3a,1).

The other real-root sectors and the full equation Φ_a(U,V) = k are not classified. The thresholds are sufficient,
not optimal.

## Proof 1: `approx.normsign.all_degree` with m = 3 (branch `cloud-agent/elementary-connections`, round 2)
- **Files:**
  - [source/proof.md](source/proof.md): Theorems 1–3, with the cubic in §5;
  - [source/sources.md](source/sources.md) and [source/README.md](source/README.md);
  - [source/experiments.py](source/experiments.py) with
    [source/experiment_results.json](source/experiment_results.json) and
    [source/experiment_output.txt](source/experiment_output.txt);
  - [source/freeze.json](source/freeze.json): the hashes at submission;
  - [source/first_campaign/](source/first_campaign/CLAIM.md): the sign lemma A1–A3, with its experiments and review.
- **Steps:**
  1. For m = 3 the ladder, multiplied by d₃(n) = 3^{v₃(n!)} ≤ √3ⁿ, gives integers P_n, Q_n with P_n/Q_n bracketing
     ρ = tan(arctan(1/a)/3). For a ≥ 1000 they satisfy Q_n ≤ Bⁿ and |ρQ_n − P_n| < (3a)⁻¹E^{−n}, with
     B = (125/6)a and E = a/(4√3); then E²/B = a/1000.
  2. Suppose U > 3aV and |Φ_a| ≤ H = 6a²/125³.
     - The mean value theorem gives |Uρ − V| < H/(2aU²).
     - Integer separation at the least N with U ≤ aE^N forces N = 1.
     - The first remainder 2(U − 3aV) ≥ 2 then leads to 2 < 2/125³ + 4/27, a contradiction.
  3. On 2V ≤ U ≤ 3aV, the identity Φ_a(U,V) + 8aV³ = (U − 3aV)(U² − 3V²) gives Φ_a ≤ −8aV³. Equality holds only at
     V = 1, U = 3a.
- **Lean:** none.
- **Audit:**
  - [source/review/REVIEW.md](source/review/REVIEW.md), with [verdict.json](source/review/verdict.json),
    [check.py](source/review/check.py), [results.json](source/review/results.json) and
    [commands.md](source/review/commands.md): a Codex reviewer distinct from the author. PASS for Theorems 1–3, with
    97,370 cubic probes.
  - The first-campaign review of A1–A3 is in [source/first_campaign/review/](source/first_campaign/review/REVIEW.md):
    PASS after a correction of presentation.
  - Consolidation check, 2026-10-26: §§4–5 of the proof and the identity in step 3 were checked by hand.
    [checks/check_angle.py](checks/check_angle.py) (output [checks/check_angle.json](checks/check_angle.json)),
    in exact arithmetic, confirms:
    - the normalization for m ∈ {3,4,5,6,7} and n ≤ 24, least possible for odd m;
    - the A1 signs for m ∈ {3,5}, n ≤ 20 and a ∈ {1,2,3,10,1000};
    - the outer sector for a ∈ {1000, 1001, 5000, 2604166, 2604167, 10⁷}. By Legendre's criterion every (U,V) there
      with |Φ_a| ≤ H is a continued-fraction convergent of ρ. No convergent with U < 10²⁰⁰ falls below H; the least
      value found at a = 1000 is about 3·10⁶, against 8a = 8000.
- **Provenance:** the branch claims no novelty.
  - The approximants are the classical hypergeometric (Padé) approximants to an m-th root (Thue, Siegel).
  - The denominator count is the one in
    [PLUS_K_PROOF.md](../../../shared/prior/ghost_parametrization/team_work/plus_k_completion/PLUS_K_PROOF.md).
  - The sign lemma is the repository's degree-4 norm-sign lemma ([ljunggren](../ljunggren/README.md), proof 2),
    carried to degree 3.
- **Novelty** (four web searches, 2026-10-26, recorded in the consolidation notes): not established.
  - The family was not found.
  - The closest published families are the "simple" (cyclic) cubic, quartic and sextic families of Thomas and of
    Lettl–Pethő–Voutier (*Simple families of Thue inequalities*, Trans. Amer. Math. Soc., 1999; known here only
    through secondary sources). Their quartic member X⁴ − tX³Y − 6X²Y² + tXY³ + Y⁴ is the degree-4
    analogue Re((X+iY)⁴) − (t/4)·Im((X+iY)⁴).
  - Φ_a is not of that kind: its discriminant is 108(a²+1)² = 3·(6(a²+1))², so its field is not Galois over Q (the
    quadratic resolvent is Q(√3)).
  - General effective measures for roots of this type may exist and were not read: Voutier, *Thue's
    Fundamentaltheorem, I*, arXiv:0805.3176, extends the method to imaginary quadratic fields.
  - The mechanism is classical; Ishan counts approximation proofs as classical in spirit
    ([approx.normsign](../../methods/approx.normsign.md)).

## Open questions
- The other two real-root sectors, and the full Thue equation Φ_a(U,V) = k.
- A lower threshold. The finite evidence above suggests that the minimum 8a already holds from a = 1000.
