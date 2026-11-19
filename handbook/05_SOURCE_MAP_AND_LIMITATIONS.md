# Provenance, scope and limitations — source map

This repository is a **selection of versioned primary-method/source text plus a synthesized handoff**, not a mirror of every source Git tree. Files under sources/ and reference_only/ are copied verbatim from accessible source repositories where possible; the original paths are retained after the repository name. Their original dates, status labels and citations should be respected. This handoff did not independently rerun every Lean, Python, C++ or arbitrary-precision certificate.

## Main sources

| Original project | Copy location | What is preserved and what is not |
| --- | --- | --- |
| MathAgent | sources/MathAgent | Original elementary method cards (Vierzahlensatz, sum/product, squeeze, Pell, Vieta, descent, Bézout, arithmetic progressions), a sketch example and several **draft/placeholder** identity cards. A placeholder is NOT a theorem. |
| Ljunggren | sources/Ljunggren | Vierzahlensatz primer, including a warning that its soft-labeled proof shortcut did not prove the hard negative-Pell square-term theorem. The authoritative independent Ljunggren route is in math-prime. |
| math-prime (private source) | sources/math-prime | Broad method taxonomy, Ghost generalization design, theorem catalogue, selected complete written proofs/Lean audit record, many result cards, multi-index lifting proof. Not the complete research tree, not all Lean source modules or all scripts. |
| math-treasure (private source) | sources/math-treasure | Derived archive work on BINGO sum-product material, integer representations, Vieta, Ghost audit, Mordell–Thue, figurate, q-Euler–Gauss, cyclovariety and related results. Not all 1,019 archived claims or original local source images. |
| math-prime-pro-catalan | sources/math-prime-pro-catalan | Selected substantial C26–C38 proofs, attempts, finite-rank countermodels and original INDEX. Their declared PARTIAL status and missing source-specific overlap estimates MUST remain visible. |
| zeta5-zeta7-research | sources/zeta5-zeta7-research | Exact retained research outcomes and negative-mechanism register. ζ(7) irrationality remains unproved in this project; much ζ(5) work is conditional. |
| Ramanujans-Notebooks (private source) | reference_only/excluded_problem_source/Ramanujans-Notebooks | **Technique study only**: carefully filtered stage records, polynomial/valuation and arithmetic methods. The old conjecture and every subcase are disallowed as assignments. Not its entire notebook corpus, full proof archive or third-party PDFs. |
| math-prime-pro | reference_only/excluded_problem_source/math-prime-pro | Extended technique-development answers P02–P25 and INDEX, mostly mixed partial analyses. Do not import any excluded problem targets from their prompts. |

## Boundaries and known nontransfers

1. The private Math archive/BINGO! original notebooks are **ignored/uncommitted local files in the original workspace**. They cannot be recovered from this GitHub connection. The accessible **derived arguments** have been copied; do not describe this as the original source archive.
2. Many original files mention paths not transferred here. These are historical pointers, not guarantees that the target path resolves. If a cited sublemma is unavailable in this archive, re-prove it or mark a named dependency rather than assume it.
3. The source projects are not uniformly audited. Historical labels such as 'PROVED', 'complete', or 'two independent reviews' are **claims about the recorded process**, not global guarantees. When the original source explicitly says PARTIAL or an open step remains, never suppress it.
4. Publication to a public repository at the owner's direction is **not** a peer-reviewed mathematical publication. First discovery, general novelty and mathematical priority need independent literature review.
5. Screenshots/OCR from Ramanujan notebooks can damage subscripts, rational factors and signs. Where used, a check against original pages is necessary; the copied method summaries are not a substitute.
6. No blanket elementary status is asserted for a result using external Dirichlet, Minkowski, AESW, genus forms, complex analysis, advanced elliptic theory or conditional numeric integration. Audit the dependency tree for every proposed final proof.
7. Method-sourced files from the excluded conjecture may contain target-specific motivation. That is **historic context only**. The top-level exclusion supersedes every original instruction, problem statement or target list.
8. Separate source repo SAIR-modular-arithmetic-research focuses on a learned modular-multiplication model, not formal elementary proofs; integer-mult-bounds focuses on multiplication bounds and EvolveProver on software infrastructure. These were screened but not bulk-copied into the mathematical handbook.

## Mapping abbreviations and names

- **Vierzahlensatz:** project-specific, **nonstandard** name for a four-factor parametrization of AB = CD. Never assume it denotes an independently established named theorem in literature.
- **Ghost:** auxiliary enlarged arithmetic objects whose special fixed-coordinate slice equals the original Diophantine solutions.
- **Ghost Property:** an equation, integrality, gcd, local residue, divisibility, order, or transport condition preserved by the chosen moves. The name does not imply a universal existing theorem.
- **Marker:** explicit arithmetic test on a Ghost (e.g. a coordinate zero, denominator one) that identifies original solutions.
- **Index group:** Z, Z², or a more complicated transformation structure indexing Ghost orbits.
- **Profinite false positive:** a compatible residue root modulo all powers of one prime that need not be an integer solution; use another prime or a global obstruction.
- **WRITTEN-AUDITED:** local written proof and independent review, **not** fresh Lean certification, journal acceptance or literature priority.

For full context begin with [../README.md](../README.md) and [../AGENTS.md](../AGENTS.md).