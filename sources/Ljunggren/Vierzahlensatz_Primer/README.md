# Vierzahlensatz Primer

This folder is a curated copy of the MathAgent material most useful to an agent studying or applying the project-specific **Vierzahlensatz** (Four Numbers Theorem), especially in work on Ljunggren's equation.

The project deliberately preserves the name `Vierzahlensatz`, although the name is nonstandard. The core parameterization is:

```text
AB = CD

A = pq,  B = rs,  C = pr,  D = qs.
```

## Recommended reading order

1. [`knowledge/methods/vierzahlensatz.md`](knowledge/methods/vierzahlensatz.md) — canonical statement, triggers, normalization, downstream moves, worked seed, and failure modes.
2. [`knowledge/examples/vierzahlensatz/kieren_mse_x2_plus_1_eq_y3.md`](knowledge/examples/vierzahlensatz/kieren_mse_x2_plus_1_eq_y3.md) — the available worked pipeline for `X^2 + 1 = Y^3`.
3. [`knowledge/methods/euclid_lemma_factor_splitting.md`](knowledge/methods/euclid_lemma_factor_splitting.md) — the coprime square-splitting step used immediately after the parameterization.
4. [`agent/workflows/wf2_vierzahlensatz_first.md`](agent/workflows/wf2_vierzahlensatz_first.md) — how MathAgent exposes and prioritizes the method.
5. [`benchmarks/problems/ljunggren_equation/problem.md`](benchmarks/problems/ljunggren_equation/problem.md) and [`allowed_inputs.md`](benchmarks/problems/ljunggren_equation/allowed_inputs.md) — exact target, domain, known solutions, and admissible proof tools.
6. [`knowledge/methods/method_combos.md`](knowledge/methods/method_combos.md), then the supporting method notes listed below.
7. The two evaluation ledgers under [`benchmarks/evaluation/runs/`](benchmarks/evaluation/runs/) — read these as examples of verification behavior, not as certified proofs.

## Included material

### Core and worked example

- Canonical Vierzahlensatz method note.
- Kieren's `X^2 + 1 = Y^3` example, which reaches `3d^2(4d^2 + 1) = c(c + 1)` and demonstrates the product-parameterization/coprimality/square-splitting pipeline.
- The Vierzahlensatz-first workflow.
- MathAgent's method-file contract and glossary.

### Supporting methods

- Euclid-lemma factor splitting.
- Cauchy/Bezout substitution.
- Squeeze sequences.
- Vieta jumping.
- Sum-to-product identities.
- Method combinations.
- Pell-equation reduction and quadratic-to-Pell collapse.

These include both methods directly named by the Vierzahlensatz note and methods that expose an `AB = CD` form or address the discriminant/Pell obstruction that can remain afterward.

### Problem and agent context

- The example-adjacent `X^2 + 1 = Y^3` benchmark and its allowed inputs.
- The exact Ljunggren benchmark `X^2 + 1 = 2Y^4` and its per-problem whitelist.
- The allowed-toolkit map and branch/workflow protocol.
- One authoritative run that rejected the Ljunggren goal and one soft ledger that marked it proven.

## Important verification caution

The soft Ljunggren ledger is **not** an authoritative proof. Its load-bearing step asserts that the only square terms in the relevant negative-Pell sequence are `1` and `169`; that statement is essentially the hard Ljunggren theorem needing proof. The authoritative ledger rejects the goal. Do not relabel that missing theorem as `vierzahlensatz`, `pell`, or a generic lemma.

Likewise, the `X^2 + 1 = Y^3` example is only a method outline; MathAgent explicitly says it does not contain the full proof.

## Provenance and scope

The files were copied byte-for-byte on 2026-08-02 from the working tree at:

```text
/mnt/d/Projects/Math/MathAgent
```

The source repository was at commit `dcb25b6cdd522c22a446da0dc75d92268de1381f` on branch `codex/okf-repository-governance`. The working tree was not clean, so this is a snapshot of the files as they existed on disk, not a clean commit export. See [`SOURCE_SNAPSHOT.md`](SOURCE_SNAPSHOT.md) for the exact status and checksums.

General MathAgent implementation, broad test suites, generated roadmaps, and unrelated knowledge notes were intentionally excluded. Source-relative links in copied files may therefore point outside this curated subset.
