# Norm-sign argument: bounded adversarial and provenance audit

Date: 2026-09-21. This is a read-only review of the current `NormSign*`
sources, their proposed final assembly, and the closest material in the
other local project. It is not a claim of historical novelty or the final
combined Lean verification receipt.

**Verdict.** No mathematical gap was found in the parameter-polynomial
sign argument or its conversion into the new base-40 bounds. The integer
approximants are unchanged. The decisive proof of their error estimate has
changed: the new route obtains a bracket from polynomial norm signs and
then bounds the error by a consecutive determinant. It does not invoke the
previous coefficient-error factorization. No equivalent parameter-derivative
proof was found in the bounded cross-project search described below.

## 1. What is the same, and what changes

The identification of sequences is exact, not merely a resemblance.
[RealCertificate.lean](../../../../../lean/Ljunggren/RealCertificate.lean:11)
and
[NormSignCertificate.lean](../../../../../lean/Ljunggren/NormSignCertificate.lean:13)
both define their integer sequences as
`IntegralCoefficientConstruction.S 239` and
`IntegralCoefficientConstruction.T 239`. Their common initial pairs are
\((1,0)\) and \((1912,2)\), with the same order-two recurrence.

The older entirely real route proves the actual error estimate from
`err*J*G = xi`, two square identities, and coefficient majorization;
see
[RealCoefficientApproximation.lean](../../../../../lean/Ljunggren/RealCoefficientApproximation.lean:102)
and
[RealCoefficientError.lean](../../../../../lean/Ljunggren/RealCoefficientError.lean:166).
The new route instead proves, for every index and every real parameter
\(a\ge1\), the alternating sign of the quartic norm of the recurrence
polynomials. A negative secant converts that sign into a root bracket at
\(a=239\). Elementary recurrence growth and determinant bounds then give
the required numerical error estimate. The change concerns a necessary
proof step, not a renaming of the old estimate.

This remains an elementary approximation strategy using the same
approximants. It is not a new descent preserving the original square
condition, and no conclusion about publication priority follows.

## 2. Sign, degree, and boundary audit

Write, at index \(n\ge1\),
\[
 W=S_n^2+T_n^2,\quad A=S_n^4-6S_n^2T_n^2+T_n^4,
 \quad B=4S_nT_n(S_n^2-T_n^2),
\]
\[
 F=A-aB,\qquad \Psi=aA+B,\qquad
 \Delta=S_{n-1}T_n-S_nT_{n-1},\qquad
 c_n=(16n^2-1)/(2n).
\]
The recurrence and derivative ladder yield
\[
 (1+a^2)WF'-[2(1+a^2)W'+aW]F=4c_n\Delta\Psi. \tag{1}
\]
The source also records the shorter companion identity
\[
 F'\Psi-F\Psi'=4c_n\Delta W^3. \tag{2}
\]
Both are explicit consequences of the same two derivative equations;
see
[NormSignAlgebra.lean](../../../../../lean/Ljunggren/NormSignAlgebra.lean:26).
The present actual assembly uses (1). Equation (2) is not an additional
arithmetic hypothesis or a different sequence.

The audited degree argument retains every needed nonzero premise:

- \(\deg S_n=n\), \(\deg T_n=n-1\), and their leading coefficients
  are positive. Thus \(\deg W=2n\), \(\deg\Psi=4n+1\), both with
  positive leading coefficients, and \(\deg F\le4n\).
- If \(d=\deg F\), the coefficient of degree \(d+2n+1\) on the
  left of (1) is
  \((d-4n-1)\operatorname{lc}(W)\operatorname{lc}(F)\).
  Since \(d\le4n\), that coefficient cannot vanish. The right side is
  nonzero because \(c_n>0\), \(\Delta\ne0\), and \(\Psi\ne0\).
  Degree comparison therefore gives \(d=2n\); leading signs give
  \(\operatorname{sign}(\operatorname{lc}F)=-\operatorname{sign}\Delta\).
- On \(a\ge1\), recurrence positivity gives \(S_n>0\),
  \(T_n\ge0\), and \(4aT_n\le S_n\). Consequently \(A>0\),
  \(B\ge0\), \(W>0\), and \(\Psi>0\). At a zero of \(F\),
  (1) forces \(F'\) to have the sign of \(\Delta\).
- A largest real zero in \([1,\infty)\) cannot have this crossing
  direction while the polynomial has the eventual opposite sign. The
  endpoint 1 is covered by the same right-hand argument. The zero
  polynomial is excluded before taking a largest zero. Index zero is
  handled separately by \(F_0=1\).

The determinant indexing is correct: the source's `delta n` belongs to
the pair at indices \(n,n+1\), and has sign \((-1)^n\). Accordingly
the norm at index \(n+1\) has sign \((-1)^{n+1}\). No parity case is
dropped. These steps are explicit in
[NormSign.lean](../../../../../lean/Ljunggren/NormSign.lean:134)
and
[NormSignActual.lean](../../../../../lean/Ljunggren/NormSignActual.lean:91).

At \(a=239\), the ratio \(T_n/S_n\) lies in the interval where the
quartic secant is negative. The factor connecting the norm to
\(\xi S_n-T_n\) is strictly positive, so alternating norm signs give
alternating errors and hence adjacent bracketing. The argument does not
assume beforehand that these pairs approximate the root.

This sign proof uses elementary real-polynomial differentiation,
continuity, the intermediate value theorem, and finiteness of the roots of
a nonzero polynomial. It should not be described as exclusively integer
divisibility.

## 3. Integrality and termination are separate obligations

Division by \(n(n+1)\) in the recurrence does not prove integrality.
The new route keeps the actual integer coefficient construction and proves
its equality with the evaluated polynomial recurrence for every integer
parameter; see
[NormSignSequenceBridge.lean](../../../../../lean/Ljunggren/NormSignSequenceBridge.lean:26).
This supplies the integer determinant needed for separation. Formal series
remain in this construction, but no convergence estimate or old error
identity is being assumed there.

After bracketing, the checked generic estimates give
\[
 S_{n+1}\ge1912\,2868^n,\qquad
 |\Delta_n|\le2\,64^n,\qquad
 |\xi S_n-T_n|<1/(956\,40^n).
\]
The positive quartic branch then reaches the first-index contradiction;
the full signed quartic reaches \(U<32000\), followed by \(U<79\).
The parallel first-pair exception is excluded by the exact quartic identity
\(\Phi(956V,V)=-4569679V^4\). The constants, strict inequalities,
and endpoints are proved in
[NormSignTermination.lean](../../../../../lean/Ljunggren/NormSignTermination.lean)
and
[FullNormSignTermination.lean](../../../../../lean/Ljunggren/FullNormSignTermination.lean).

Alternating signs of the approximants do **not** imply that the target's
integer remainders \(UT_n-VS_n\) alternate forever. Their finite sign-exit
theorem remains only a structural result. The new exclusion uses integer
separation and the quartic error upper bound. Similarly, the positive
\(28561\) exclusion retains its cone restriction; it does not exclude
the outside-cone solution \((U,V)=(2,3)\).

## 4. Bounded cross-project comparison

The search covered Markdown, TeX, text, and Python files in the other
project's `Elementary_Proof`, `docs`, and `Attempts`, followed by reading
the closest matches. Search terms included Wronskians, parameter
derivatives, norm signs, largest roots, degree/leading-coefficient sign
arguments, and the distinctive recurrence/formula notation. It was not
an exhaustive search of unextracted PDFs or external literature.

| Closest existing material | Its actual mathematical role | Comparison |
|---|---|---|
| [Wronskian-Gaussian height selector](/mnt/d/Projects/Ljunggren/Elementary_Proof/wronskian_gaussian_height_selector.md:53) | Discrete identities \(H_{ad}W_{bc}=7WH+6Z\), selected gcd allocations, and Gaussian norm factors | No parameter derivative or all-degree approximant norm sign; a divisor-size selector remains unproved. |
| [Mixed Somos height gate](/mnt/d/Projects/Ljunggren/Elementary_Proof/mixed_somos_height_descent_gate.md:37) | Products and signed determinants of the two mutations | These determine products and allocations, not a smaller neighbor; they are not (1) or (2). |
| [Wronskian/conjugate synchronization gate](/mnt/d/Projects/Ljunggren/Elementary_Proof/wronskian_conjugate_square_transfer_gate.md:48) | Tests whether those discrete identities force the missing conjugate square | It explicitly identifies a missing cross-layer coupling, absent from the new approximation argument. |
| [Earlier quartic proof, §1.2](/mnt/d/Projects/Ljunggren/Elementary_Proof/proof.md:109) and [finite-descent presentation, §2](/mnt/d/Projects/Ljunggren/Elementary_Proof/proof_finite_descent.md:230) | Differentiate the Padé remainder in its argument, then integrate along a complex arc | Same established approximation family, but a different necessary error proof. Its finite quartic/continued-fraction tail is not the new base-40 argument. |

No exact parameter-derivative root-sign argument was located in the
searched local material. This negative search result supports only that
bounded statement, not a claim that the argument is historically new.

The new chain contains no moving quadratic-residue character, Pell-index
factor-base sieve, or Euclidean parity routine for reciprocity. The other
project's
[deep-descent audit](/mnt/d/Projects/Ljunggren/Elementary_Proof/proof_deep_descent_audit.md:254)
explicitly describes such a moving Jacobi calculation. Thus the present
proof has a different decisive mechanism from that Lin–Luo-style route.
The qualifications about the unavailable full original Lin–Luo paper in
the earlier [source comparison](../../../../../Archive/ljunggren_working_2026_09/lin_luo_comparison.md:17)
still apply; this bounded pass did not reopen that literature audit.

## 5. Verification boundary

A read-only traversal of the current `CompleteNormSign` source imports
found 26 local modules. It excludes `RealCertificate`,
`RealCoefficientApproximation`, `RealCoefficientError`,
`RealErrorFactorization`, `BinomialConstruction`, `BinomialError`,
`CoefficientSigns`, and every older completed-classification module.
It retains the integer construction, generic separation and recurrence
lemmas, root/secant facts, and the original/Vierzahlensatz reductions.
An import graph is evidence of source separation, not a substitute for
checking the resulting theorem terms.

The six structural modules completed under this agent's ownership contain
40 declarations. They compiled, and fresh dependency checks returned only
subsets of `propext`, `Classical.choice`, and `Quot.sound`.
The current actual norm-sign assembly was reviewed mathematically here;
the parent is compiling and auditing its combined theorem dependencies.
**The global final receipt remains pending at the time of this report.**
No all-project completion claim is made by this audit.

## 6. Assembly and collaborator update

After the preceding review, the assembled
[CompleteNormSign.lean](../../../../../lean/Ljunggren/CompleteNormSign.lean)
compiled cleanly. It includes the remaining positive \(28561\) exclusion,
the complete original classification, and the normalized and screenshot
endpoint theorems. The current
[main lemma](elementary_norm_sign_lemma.md) and
[proof](elementary_norm_sign_proof.md) describe
this replacement. The whole-project `verify.sh` run and all six audits subsequently passed. The historical pending paragraph above records the status at the initial review. The definitive [final receipt](../../../../../lean/FINAL_AUDIT.md) covers70 modules and1,720 local declarations; the focused audit covers302 new declarations and17,663 recursive dependencies, with all gates clean.

The [Pro collaborator's completed review](../../../../../Archive/ljunggren_working_2026_09/collaborator_norm_sign_review.md)
independently agrees with the sign, degree, endpoint, bracketing, and
base-40 termination arguments, including the signed-sector constants.
Its discussion of integrality correctly emphasizes that arbitrary-index
integrality is needed before concluding that the selected index is one.
The review's symbolic checks and finite regression tests are separate
from its all-index argument and from the local Lean verification.

The optional normalized-function interpretation, continued fraction, and
alternative two-series integrality construction in its Sections 4–5 are
not adopted. The actual proof retains the already checked integer
coefficient construction. The completed replacement is still an
approximation proof; a purely integer descent preserving the source
equations remains unproved. This distinction is also reflected in the
updated [arithmetic research receipt](../../../../../Archive/ljunggren_working_2026_09/arithmetic_replacement_research.md).
