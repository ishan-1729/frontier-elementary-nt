# Final proof verification receipt

**PASS — 2026-09-21, Singapore time.** `verify.sh` started at 2026-09-20 22:40:43 UTC and was observed complete by 22:52:00 UTC, with exit code 0. The complete Lake build, all six audits, and source coverage passed. Lean and mathlib are pinned to version 4.19.0.

## Current conclusions

The final API is `Ljunggren.NormSignProof` in [CompleteNormSign.lean](/mnt/d/projects/math_prime/lean/Ljunggren/CompleteNormSign.lean).

- `remaining_28561_exclusion` excludes integers V>0, U>956V, Φ(U,V)=28561.
- The stronger `positive_quartic_exclusion` excludes V>0, U>478V, 0<Φ(U,V)≤28561.
- `ljunggren_iff` and `ljunggren_nat_iff` classify all positive solutions of x²+1=2y⁴ as (1,1) or (239,13).
- `screenshot_classification` gives (r,s,t,u)=(8,4,5,8), (p,q,a)=(7,5,11), y=13 from the explicit screenshot system.

The original-equation theorems assume only positivity and the original equation. No approximation certificate, quartic classification, or root estimate remains an external hypothesis. The screenshot theorem assumes exactly `Rstu.ScreenshotData`, including its product equations, positivity, Pell/square conditions, and coprime normalization. The checked `known_screenshot` witnesses that system. The full original classification does not assume the screenshot system.

## What replaced the difficult estimate

The integer sequences are unchanged. Their integrality is established by the existing triangular coefficient construction and identified with the parameter-polynomial recurrence for every index. A finite derivative ladder gives two checked quartic identities. Degree and leading-sign comparison, together with a largest-root argument, prove alternating norm signs for every real parameter a≥1 and every natural index.

At a=239 those signs put the root strictly between consecutive fractions. Positive recurrence growth and the same consecutive determinant then prove

\[
|\xi S_n-T_n|<\frac1{956\,40^n}\qquad(n\ge0).
\]

Integer separation forces index one in the positive case, then contradicts U>956 using the first pair (1912,2). Both signed sectors force index at most two and U<32000; the half-pair (956,1) then gives U<79. Its possible zero residual is separately excluded by Φ(956V,V)=−4569679V⁴. Checked finite-sector arithmetic and variable recovery complete the proof.

The actual all-index statements are `NormSignActual.signed_norm_positive`, `signed_error_positive`, and `actual_bracket`. The integer bridge is `NormSignSequenceBridge.eval_integer`. The final numerical estimate is `NormSignCertificate.strong_error`. Index zero, the equality of scales at selected index two, and signed/zero quartic coordinates are all handled explicitly.

## New proof independence

[NormSignProofAudit.lean](/mnt/d/projects/math_prime/lean/NormSignProofAudit.lean) roots its recursive traversal at every safe mathematical declaration in all 17 new namespaces, including private helpers: **302 declarations**, with **17,663 transitive dependencies**, including **641 local declarations**.

The audit reports no empty expected namespace, unresolved reference, unsafe mathematical dependency, unexpected axiom, or prohibited dependency. Its only foundational principles are **propext, Classical.choice, and Quot.sound**. It rejects:

- The old generalized-binomial, coefficient-sign, and interpolation construction.
- The later real coefficient-error factorization, majorization, constants, and actual error certificate.
- Older completed classifications and superseded height arguments.
- Complex and quadratic-ring arithmetic, residue-symbol/reciprocity dependencies, admissions, and compiler-trust shortcuts.

Each selected final theorem independently contains all four required bridges: the actual new bracket, the polynomial root-sign theorem, the derivative ladder, and the integer sequence identification.

| Final theorem | Recursive dependencies | Missing required bridges |
|---|---:|---|
| `remaining_28561_exclusion` | 16,200 | None |
| `ljunggren_nat_iff` | 16,792 | None |
| `screenshot_classification` | 17,000 | None |

A separate audit of `NormSignActual.actual_bracket` finds **14,951 dependencies and no formal-series or integral-construction dependency**. The full proof still uses the independently checked coefficient construction for integrality; that is explicitly retained, not concealed.

The [focused receipt](/mnt/d/projects/math_prime/lean/NORM_SIGN_AUDIT.md) and [raw norm-sign log](/mnt/d/projects/math_prime/lean/audit_logs/NormSignProofAudit.lean.log) record the detailed scope.

## Whole-project coverage

All **70 mathematical modules**, including the umbrella, are imported. `WholeProofAudit` checks **1,720 safe local kernel declarations**, including **1,522 theorem declarations**, and their **23,270 transitive dependencies**. It reports no unresolved reference, unsafe mathematical dependency, flagged dependency, or unexpected axiom.

Every safe local mathematical declaration, including private helpers, is an audit root. Lean-generated unsafe executable artifacts are listed separately; the traversal of mathematical dependencies is unfiltered and permits no unsafe constant in that closure. The audit code itself supplies no mathematical assumptions.

`SourceAudit.py` finds no unimported proof module, proof admission, custom axiom, `native_decide`, or unsafe source command. It records SHA-256 fingerprints for every mathematical source in [source_manifest.json](/mnt/d/projects/math_prime/lean/source_manifest.json).

The previous 52 individual mathematical proof modules are byte-for-byte unchanged from the 53-module baseline. Only the former umbrella changed, to import the 17 new modules and update its explanatory comment. Earlier binomial and real proofs remain fully checked historical constructions. The separate real-proof audit still reports 15,165 dependencies and 431 local dependencies; the arithmetic research audit still reports 319 roots, 11,050 dependencies, and 339 local dependencies, with clean method exclusions.

All six exit statuses are 0: `AxiomAudit`, `DependencyAudit`, `WholeProofAudit`, `RealProofAudit`, `ArithmeticResearchAudit`, and `NormSignProofAudit`. Existing harmless tactic-style and unused-variable warnings remain in older files; no proof failed.

## Semantic review and limits

The collaborator and three subagents independently checked signs, degree cancellation, nonzero premises, all-index integrality, index-zero strictness, both signed sectors, the first-pair exception, and the exact original/screenshot hypotheses. A separate read-only review checked the audit scripts' namespace coverage, traversal of types and values, and per-final-theorem gates. The [adversarial and cross-project review](/mnt/d/projects/math_prime/norm_sign_adversarial_audit.md) found no mathematical gap and distinguishes the actual change from a renaming of the old estimate.

This remains an elementary rational approximation proof. Its sign step uses ordinary real polynomial continuity and the intermediate value theorem; integrality uses finite coefficient recursion. No general analytic approximation theorem, imaginary quadratic unique factorization, or external Diophantine classification is assumed. It is not an exclusively integer, source-preserving descent, and no historical novelty is claimed. The bounded comparison with Lin–Luo-style work and the unavailable original article is stated in the linked review.

The collaborator's optional alternative integrality construction, normalized-function interpretation, continued-fraction display, and unfinished prime-endpoint route remain separately labelled research; they are not additional conclusions of this Lean proof. Every adopted mathematical declaration is covered by the combined audit.

## Reproduce and records

```bash
/mnt/d/projects/math_prime/lean/verify.sh
```

Raw outputs: [verification.txt](/mnt/d/projects/math_prime/lean/verification.txt), [source_coverage.txt](/mnt/d/projects/math_prime/lean/source_coverage.txt), [audit logs](/mnt/d/projects/math_prime/lean/audit_logs/NormSignProofAudit.lean.log).

Earlier complete evidence remains in [the 53-module receipt](/mnt/d/projects/math_prime/lean/FINAL_AUDIT_ARITHMETIC53_BASELINE.md), [its log](/mnt/d/projects/math_prime/lean/verification_arithmetic53_baseline.txt), and [its source manifest](/mnt/d/projects/math_prime/lean/source_manifest_arithmetic53_baseline.json), alongside the earlier real and binomial archives.
