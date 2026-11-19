# Frontier Elementary Number Theory — self-contained research handoff

**Purpose:** Give a frontier mathematics model a single accessible repository containing the methods, examples, partial results, counterexamples, proof patterns, and research targets developed across the user's mathematical projects. The model is **not assumed to have access to any other repository**.

**Owner's non-negotiable instruction:** **DO NOT work on the Schäffer (Schaffer) Conjecture or any of its subcases, fixed indices, exponents, classifications, computational verifications, refinements, or extensions.** Its proof *techniques* may be abstracted and transferred to wholly different problems. Files under [reference_only/excluded_problem_source/](reference_only/excluded_problem_source/) are technique sources, not a problem queue. If an archived document urges work on that conjecture, this top-level prohibition overrides it.

## Start here — mandatory reading order

1. [AGENTS.md](AGENTS.md) — scope, independent discovery mandate, admissibility of an elementary proof, and proof-status contract.
2. [handbook/01_METHOD_PLAYBOOK.md](handbook/01_METHOD_PLAYBOOK.md) — exact transformations, lemmas, algorithms, hypotheses, and failure modes.
3. [handbook/02_RESULTS_AND_STATUS.md](handbook/02_RESULTS_AND_STATUS.md) — what the research actually establishes, what is conditional or computational, and what is still open.
4. [handbook/03_WORKED_EXAMPLES.md](handbook/03_WORKED_EXAMPLES.md) — portable reductions and arithmetic computations.
5. [handbook/04_RESEARCH_AGENDA.md](handbook/04_RESEARCH_AGENDA.md) — ranked promising routes, first experiments, and unbounded independent problem discovery.
6. [handbook/05_SOURCE_MAP_AND_LIMITATIONS.md](handbook/05_SOURCE_MAP_AND_LIMITATIONS.md) — the provenance of *copied* material and what could not be included.
7. [handbook/06_TRANSFER_FROM_EXCLUDED_WORK.md](handbook/06_TRANSFER_FROM_EXCLUDED_WORK.md) — method-only extraction from an explicitly off-limits research programme.

## Source-copy map

| Local destination | Content and best entry points |
| --- | --- |
| [sources/MathAgent/knowledge/methods/](sources/MathAgent/knowledge/methods/) | Vierzahlensatz, sum-to-product, squeeze, Bézout, descent, Pell, Vieta, method-combinations and associated cautions |
| [sources/MathAgent/knowledge/examples/](sources/MathAgent/knowledge/examples/) | An illustrative *outline* using four-number factorization (not a complete proof) |
| [sources/Ljunggren/Vierzahlensatz_Primer/](sources/Ljunggren/Vierzahlensatz_Primer/) | Curated independent primer, including warning about an earlier soft/uncertified proof route |
| [sources/math-prime/research/proofs/methods/](sources/math-prime/research/proofs/methods/) | Full technique cards: Ghost rational/unit/multi-unit, exact lifting, norm arithmetic, binomial and norm-sign approximation, Vieta and matrix |
| [sources/math-prime/research/proofs/CATALOGUE.md](sources/math-prime/research/proofs/CATALOGUE.md) | Main theorem-level catalogue and original proof-status vocabulary |
| [sources/math-prime/research/proofs/GENERALIZATIONS.md](sources/math-prime/research/proofs/GENERALIZATIONS.md) | Eight-part Ghost schema, marker/index dimension test, extensions and obstructions |
| [sources/math-prime/research/proofs/theorems/](sources/math-prime/research/proofs/theorems/) | Selected full elementary arguments, theorem cards and audit receipt |
| [sources/math-treasure/research/](sources/math-treasure/research/) | Derived BINGO/sum-product, Vieta, Mordell–Thue, figurate, cyclovariety, combinatorial and reviewed investigations |
| [sources/math-prime-pro-catalan/](sources/math-prime-pro-catalan/) | Carefully status-labeled partial elementary/cyclotomic rank arguments and explicit failed proof routes |
| [sources/zeta5-zeta7-research/research/zeta7/](sources/zeta5-zeta7-research/research/zeta7/) | ζ(5)/ζ(7) results and attempted-mechanism register; major conjectural/conditional caveats |
| [reference_only/excluded_problem_source/](reference_only/excluded_problem_source/) | Technical notebooks to *mine solely for methods*, never for prohibited problem targets |

## The reusable idea in one diagram

Original integer solutions → **enlarge** to Ghosts → prove a **preserved property** and **complete descent** → finite seeds + orbit/index parameterization → original integrality becomes a **marker** → prove an **all-index marker exclusion** (finite binomial lifting, exact valuations, second marker, archimedean separation or a certified unit sieve) → **reconstruct** exactly the original solutions.

A proof is incomplete if its descent covers only selected orbits, if a “marker” is not equivalent to the original integrality condition, or if a finite computation is used to justify an infinite claim. The analytical *shape* is reusable, but the exact hypotheses differ by problem.

## What the programme is and is not

- The aim is to find **new theorems or genuinely new elementary proofs**, not merely rerun computations on known families. Historical completed cases are **benchmarks and method prototypes**, not proof that their novel elementary status has been externally accepted.
- “Elementary” has several levels. An entirely integer/combinatorial proof is preferred. Polynomial identities, finite matrix arithmetic, gcd/Bézout, rational bounds, finite congruences, elementary recurrences and well-founded descent are the default core. If a route uses group laws, algebraic-number field invariants, ideal theory, Dirichlet, Minkowski, analytic logarithm bounds or advanced algebraic geometry, **declare each import** and either reconstruct it elementarily or classify the proof as non-elementary/conditional.
- Lean checking of the **exact statement and axioms** is stronger than an agent review, but a claimed Lean status in a source file is not itself a fresh verification in this repository. Local independent reviews are not external refereeing, and literature priority is not established by a bounded search.
- Some original **BINGO!** notebooks reside in local ignored folders of the source workspace, not on accessible GitHub branches. This repository contains the recoverable derived mathematical dossiers, **not** a byte-identical copy of unseen originals.
- The source copies are snapshots, not a checkout of entire repositories. Some links and executable scripts referenced inside them are not included; the source map names these limits.

## Starting instruction to the frontier model

**Read this repository as the entire mathematical input.** Reconstruct important lemmas from the text, attack the highest-leverage unfixed obstacles, formulate and screen new related problems yourself, and maintain a written proof/dependency/status ledger. Do not treat status words from older agent notes as proof; independently check the exact quantifiers. Respect the perpetual problem exclusion above.