# The gcd invariant of (a+r)/b + (b+r)/a = k, and Conjecture 4.2 of Miller–Nikolakopoulos–Srinivasan
Status: WRITTEN-AUDITED     Methods: [involution.vieta](../../methods/involution.vieta.md)
Generalizes: MNS Theorem 1.2, the case r = 1 (not a card)     Generalized by: none.
Related: [vieta_conics](../vieta_conics/README.md) gives the set of k for each r.

## Statement
Let r ≥ 1. A *solution* is a pair of positive integers (a,b) with a² + b² − kab + r(a+b) = 0 for an integer k,
that is (a+r)/b + (b+r)/a = k. Write g = gcd(a,b), a = ga′, b = gb′, and Q(a,b) = r(a+b)/g². The Vieta jumps are
σ₁(a,b) = (kb − a − r, b) and σ₂(a,b) = (a, ka − b − r); a solution is fundamental if neither jump lowers a
coordinate.

**Theorem.**
(a) Q = k·a′b′ − a′² − b′². In particular Q is a positive integer, and g = r(a′ + b′)/Q.
(b) Conversely, if a′, b′ are coprime and Q := k·a′b′ − a′² − b′² > 0 divides r(a′ + b′), then (ga′, gb′) with
    g = r(a′ + b′)/Q is a solution with invariant Q.
(c) Jump rule. Put e = gcd(k + 2, a′ + b′). Then Q(σ₁(a,b)) = Q(σ₂(a,b)) = (k + 2)Q/e², and the new e is
    (k + 2)/e. So on every orbit Q takes the two values Q and (k + 2)Q/e² alternately.
(d) Every orbit contains a fundamental solution. A fundamental solution has |a − b| < r/2 and
    max(a,b) < (3 + √5)r/2 < 2.62r, and 3 ≤ k ≤ 2r + 2 for every solution.
(e) The set V(r) of values of Q is the union, over the fundamental solutions, of {Q, (k + 2)Q/e²}: a finite
    computation for each r, tabulated for r ≤ 100.
    - r = 2: for every positive solution of (a+2)/b + (b+2)/a = k, 2(a+b)/gcd(a,b)² ∈ {1,2,3,4,5,8}. The values
      are 1 and 5 for k = 3, 2 and 3 for k = 4, and 4 and 8 for k = 6, each on infinitely many solutions. This is
      Conjecture 4.2 of MNS (arXiv 2605.19083).
    - r = 1 gives MNS Theorem 1.2: the values {1,5} for k = 3 and {2,3} for k = 4.
(f) For each divisor d of 2r, (2r/d, 2r/d) is a fundamental solution with k = d + 2, contributing
    {d, d(d+4)/gcd(d,2)²}.

## Proof 1: `involution.vieta`
- **Origin:** branch `origin/claude/math-prime-cross-domain-g1vwsn`, tip `26c718ac`, folder
  `cloud_agent_claude/13_gcd_invariant_vieta/` (track 13, the programme lead's own work, 18 October 2026). The
  files are copied byte for byte into [source/](source/).
- **Files:** [T1_gcd_invariant.md](source/13_gcd_invariant_vieta/theorems/T1_gcd_invariant.md); the script
  [gcd_invariant.py](source/13_gcd_invariant_vieta/scripts/gcd_invariant.py), with its output
  [gcd_invariant.txt](source/13_gcd_invariant_vieta/evidence/gcd_invariant.txt): symbolic identities, V(r) for
  r ≤ 100, and a brute force over all solutions with max(a,b) ≤ 20,000 for r ≤ 20.
- **Steps:**
  1. (a): substitute a = ga′, b = gb′ and divide by g.
  2. (c): with c = (k + 1)b′ − a′ and d = a′ + b′, one has a* : b = c : d and gcd(c,d) = e. The identity
     kcd − c² − d² = (k + 2)(k·a′b′ − a′² − b′²) then gives the new invariant, and a second application returns Q.
  3. (d): a jump that lowers a + b exists until the solution is fundamental, and a fundamental solution satisfies
     (a − r)(b − r) < 5r²/4.
  4. (e) for r ≤ 2, by hand: every fundamental solution is diagonal, (t,t) with t = 2r/(k − 2).
- **Lean:** none.
- **Audit:**
  - The branch's review [13__T1_gcd_invariant.md](source/reviews/13__T1_gcd_invariant.md): PASS, fidelity MATCH.
    The referee's own C enumerator found every solution with max(a,b) ≤ 10⁶ for r ≤ 30 and with max(a,b) ≤ 10⁵
    for r ≤ 100; 4,000 random solutions had r up to 2.5·10²⁹; V(r) was recomputed in PARI for r ≤ 100. The
    author's changes after the review ([REVISIONS.md](source/13_gcd_invariant_vieta/notes/REVISIONS.md)) are
    presentational and were not reviewed again.
  - Our check, 2026-10-26 (rev1): CONFIRMED. A brute force over all a, b ≤ 6000 and a divisor search over
    a ≤ 3·10⁵ reproduce the r = 2 value sets, and the general-r formula agrees with brute force for every r ≤ 12
    ([t13_out.txt](../../../../research/shared/consolidation_2026-10-26/review_checks/rev1/t13_out.txt),
    [t13_general_out.txt](../../../../research/shared/consolidation_2026-10-26/review_checks/rev1/t13_general_out.txt)).
    The consolidation re-checked (a) and (c) symbolically, and the jump rule on all 2,756 solutions with
    a, b < 400 and r ≤ 40
    ([t13_identities.py](../../../../research/shared/consolidation_2026-10-26/a_tracks/work/checks/t13_identities.py)).
- **Provenance and novelty:**
  - MNS (v1 18 May 2026, v2 27 May 2026, the latest) state Conjecture 4.2 as open. No proof was found up to
    26 October 2026: the branch's logs L1–L6, its referee's W1–W11, and ours
    ([novelty_web_log.md](../../../../research/shared/consolidation_2026-10-26/a_tracks/work/novelty_web_log.md),
    Q2: 18 queries, OpenAlex counts no citing work). Alekseyev–Khomovsky, arXiv 2610.00808 (September 2026), use
    the same conic for coprime pairs only, with no k-sets and no gcd invariant.
  - The referee's observation: the k = 3 and k = 4 parts of the conjecture reduce to MNS Theorem 1.2 by the scaling
    (a,b,r) ↦ (2a,2b,2r). The new part is the k = 6 orbit of (1,1), with the values {4, 8}.
  - The jump rule (c) and V(r) for general r: no match. Remark 4 of the source reads (c) as multiplication by 1 + θ,
    θ² = kθ − 1, an element of norm k + 2; that explains the rule but does not make it known.
  - MNS prove Theorem 1.2 through Fibonacci and Lucas closed forms. The proof here uses one identity, valid for
    every r and k.

## Open questions
- A closed form for V(r).
- Three variables: for (a+r)/b + (b+r)/c + (c+r)/a = k, Vieta jumping need not preserve integrality (MNS §4.2),
  and no analogue of (c) is known.
