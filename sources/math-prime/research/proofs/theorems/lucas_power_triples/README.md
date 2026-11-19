# Unit Lucas sequences: no three r-th powers at index spacing divisible by r
Status: WRITTEN-AUDITED (one external audit, below)     Methods: not yet classified (a gap principle on r-th roots: archimedean, in the spirit of the [approx](../../methods/README.md) family)
Generalizes: none     Generalized by: none

## Statement
Let U₀ = 0, U₁ = 1 and U_{n+2} = A·U_{n+1} − Q·U_n, with either Q = −1 and A ≥ 2, or Q = 1 and A ≥ 3. Examples are
the Pell numbers, (A,Q) = (2,−1), and the even-index Fibonacci numbers F_{2n}, (A,Q) = (3,1).
- **D1.** For all r ≥ 2, n ≥ 1 and k ≥ 1, the integers U_n, U_{n+rk} and U_{n+2rk} are not all r-th powers.
- **D2 (fixed sieves).** Fix such a sequence and r ≥ 2. Fix finitely many moduli, with least common multiple M, and
  finitely many congruences on the index, with least common multiple q. Let a ≥ 1 be an index with U_a an r-th
  power; a = 1 always works. Let T be the period of (U_{n+1}, U_n) modulo M, and put L = lcm(r, T, q). Then:
  - every U_{a+jL} is ≡ U_a (mod M), so it passes every r-th-power residue test at the chosen moduli;
  - every index a+jL satisfies the index congruences that a satisfies;
  - among any three consecutive j, at least one U_{a+jL} is not an r-th power.

D1 restricts the pattern of power indices. It gives no finiteness and classifies nothing: an infinite set of
indices can avoid all three-term progressions.

## Proof 1: root recurrence and monotonicity (branch `cloud-agent/elementary-connections`, first campaign)
- **Files:**
  - [source/proof_draft.md](source/proof_draft.md): D1 is in §§2–4, D2 in §5, and the attempts in §6;
  - [source/experiment.py](source/experiment.py) with
    [source/experiment_output.json](source/experiment_output.json);
  - [source/COMMANDS.md](source/COMMANDS.md) and [source/source_status.md](source/source_status.md).
- **Steps:**
  1. U_j = (α^j − β^j)/(α − β) with α > 2 and |β| < 1. The traces T_j = α^j + β^j are integers, and α^k is
     irrational.
  2. With c = (α−β)^{−1/r}: |U_j^{1/r} − c·α^{j/r}| ≤ (4/3)·c·α^{−(2−1/r)j}.
  3. Suppose U_n, U_{n+h}, U_{n+2h}, with h = rk, are x^r, y^r, z^r. The integer E = z − T_k·y + Q^k·x is a
     combination of the three errors, so |E| < 8/(9√2) < 1. Hence z = T_k·y − Q^k·x.
  4. The step-h recurrence then gives (T_k − st)^r + s^r·t^r = T_{rk}, with t = x/y in (0,1) and s = Q^k. This
     equation has exactly one root in (0,1), namely α^{−k}, which is irrational. Contradiction.
  5. *D2.* The state map (u,v) ↦ (Au − Qv, u) is a permutation modulo M, because its determinant is ±1. So the
     states are purely periodic, and D1 applies to every three consecutive terms of the progression.
- **Lean:** none.
- **Audit:**
  - [source/review/review.md](source/review/review.md), with [verdict.json](source/review/verdict.json),
    [check.py](source/review/check.py) and [output.json](source/review/output.json): a Codex reviewer distinct from
    the author. The same report also reviews the branch's polynomial_methods claims. PASS for D1 and D2. It checked
    959,364 triples (A ≤ 12, r ≤ 12, indices ≤ 300).
  - Consolidation check, 2026-10-26: the proof was read in full.
    [checks/check_dynamics.py](checks/check_dynamics.py) (output
    [checks/check_dynamics.json](checks/check_dynamics.json)) lists every r-th power U_n with n ≤ 400 and
    r ≤ 12, for 27 parameter pairs with A ≤ 15. It finds no forbidden triple.
  - Powers are rare. Only nine (sequence, r) pairs have two power indices, for example indices 1 and 7 for Pell
    squares (P₇ = 169), and indices 1 and 6 for squares in F_{2n} (F₁₂ = 144). So the finite checks are weak
    evidence; the proof carries the claim.
- **Provenance:** the branch claims no novelty. The closed formula, the periodicity and the exponential
  approximation are standard.
- **Novelty** (three web searches, 2026-10-26, recorded in the consolidation notes):
  - No source for D1 was found.
  - For particular sequences D1 follows from complete classifications of powers or squares: Fibonacci and Lucas
    numbers (Bugeaud–Mignotte–Siksek, Ann. of Math. 163 (2006)), and squares in Lucas sequences with P and Q odd
    (Ribenboim–McDaniel). For general A and r ≥ 3 no classification is known, so D1 is not implied in general.
  - Each such sequence has only finitely many perfect powers (Pethő; Shorey–Stewart; both by linear forms in
    logarithms, and cited from a secondary source, not read). That bound is effective but huge, and it says nothing
    about index patterns. D1 is a different, explicit statement.
  - The mechanism is a gap principle, classical in spirit: an integer of absolute value below 1 must vanish.
  - Confidence that the statement is new: medium to low. It is a weak statement.
- **Relation to the programme:** D2 makes rigorous, for r-th power markers, the dead end "a finite congruence sieve
  alone: large indices imitate small ones" ([GENERALIZATIONS.md §3.6](../../GENERALIZATIONS.md)). The `lifting.*`
  engines avoid it with exact p-adic drift laws.

## Open questions
- An arithmetic rule that forces power indices into a pattern D1 forbids, or that bounds the largest power index
  (the branch's NEXT_TARGETS §5).
- Gaps not divisible by r: the proof needs the root ratio α^k to be a unit with integer trace.
