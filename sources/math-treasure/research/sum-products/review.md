# Independent review receipt

Date: 2026-10-26. Board: docs/goals/archive-research/state.yaml.
Proof package: T026/R09. Independent Judge: T026/R13, /root/outcome_audit.

The coordinating PM relayed the R13 review: the Judge read the complete proof
and verification files, visually inspected all five original images, read all
six local source notes, checked all eleven source hashes, and approved the
written mathematical proofs and exact verification ranges. This is independent
written review, not formal verification or an exhaustive novelty assessment.

Reviewed claims include weighted indexed subsets and cardinality recurrence,
the coefficient 64, completeness versus uniqueness, formal and analytic
infinite-product hypotheses, mixed radix, the Omega inverse kernel and corrected
Chebyshev cosine formula, and the positive-irrational Beatty partition proof.
The exact verifier covers 210 subset/prefix cases, 120 mixed-radix cases,
2470 inverse entries, and 437 cosine traces, with ranges in the verifier and
saved record. It also checks the Omega reduction for every integer
0≤b≤a≤n with 2≤n≤20.

R13 requested one provenance correction: integrality of a,b and the finite
integer k range are reconstructed restrictions, not explicit hypotheses in
the original Omega image. README.md now says so and records the review's
real-parameter counterexample n=2, a=b=1/2 (left side 1/2, right side −3/4).
The proved theorem and verifier were unchanged by this clarification.

Passing commands for the proof package, rerun/checked in the independent
review as reported by the PM:

- `python3 experiments/sum-products/verify.py`
- `git diff --check -- research/sum-products experiments/sum-products`

After this provenance repair, the scoped diff check passed again. The numerical
verifier was not rerun for this documentation-only change, per the PM's scope.

The linked Kainui/Brilliant statements, broader Euler reflection connections,
prime/quadratic-form, Goldbach and Landau bridges, and undefined general
complete-sequence multiplication remain open. Review of this package does
not discharge them or establish archive-wide completion.
