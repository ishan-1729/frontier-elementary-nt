# approx.normsign.all_degree: the norm-sign ladder in every degree m ≥ 2
Parent: [approx.normsign](approx.normsign.md)     Children: none

## What it does
It carries the Wronskian sign certificate of the second Ljunggren proof from degree 4 to every degree m ≥ 2. It also
supplies the integer normalization, which degree 4 does not need.
1. **Ladder.** S₀ = 1, S₁ = 2mx, T₀ = 0, T₁ = 2, and Y_{n+1} = α_n·x·Y_n + β_n·Y_{n−1}, with
   α_n = 2m(2n+1)/(n+1) and β_n = 4(m²n²−1)/(n(n+1)).
2. **Sign (A1).** Put A + iB = (S_n + iT_n)^m, F_n = A − xB, Ψ_n = xA + B and W = S_n² + T_n². Then
   F′Ψ − FΨ′ = m·c_n·Δ_{n−1}·W^{m−1}, where c_n = 2(m²n²−1)/(mn) and Δ_n = S_nT_{n+1} − S_{n+1}T_n is a constant
   of sign (−1)ⁿ. Hence deg F_n = n(m−2), and (−1)ⁿF_n(a) > 0 for every real a ≥ 1.
3. **Brackets (A2).** The fractions T_n(a)/S_n(a) alternate around ρ_m(a) = tan(arctan(1/a)/m), with
   |ρS_n − T_n|·S_{n+1} < |Δ_n|.
4. **Closed form and denominators** (round 2, Theorem 1).
   S_n + iT_n = mⁿ Σ_j binom(n+1/m, n−j)·binom(n−1/m, j)·(x+i)^{n−j}(x−i)^j.
   Multiplying by d_m(n) = ∏_{p|m} p^{v_p(n!)} makes both polynomials integral. For odd m this multiplier is the least
   possible, and d_m(n)^{1/n} → c_m = ∏_{p|m} p^{1/(p−1)}.
5. **Measure** (round 2, Theorem 2). Put B = c_m·2m(a + √(a²+1)), E = 3a/(4m·c_m) and κ = log B / log E. For
   integers a ≥ 1 with E > 1, every U ≥ 1 and every V satisfy
   |ρ − V/U| > 1/(2B²U·max(1, 2U/(ma))^κ).

## Ingredients (with the file that proves each)
All in [theorems/cubic_angle_sector/source/](../theorems/cubic_angle_sector/source/proof.md):

| Ingredient | Proof | Check |
|---|---|---|
| Derivative ladder, Δ_n, the Wronskian identity, the sign | [first_campaign/CLAIM.md](../theorems/cubic_angle_sector/source/first_campaign/CLAIM.md) A1; [proof.md](../theorems/cubic_angle_sector/source/proof.md) §3 | first_campaign/experiments.py; first_campaign/review/check.py; [check_angle.py](../theorems/cubic_angle_sector/checks/check_angle.py) |
| Alternating brackets and the determinant error bound | CLAIM.md A2; proof.md §§3–4 | review/check.py |
| Closed binomial form; exact denominators | proof.md §§1–2 | experiments.py; check_angle.py |
| Explicit approximation bound | proof.md §4 | review/check.py |

## Requirements: when it applies
- The target must reduce to approximating tan(arctan(1/a)/m) for an integer a ≥ 1. Equivalently, it concerns a form
  Re((U+iV)^m) − a·Im((U+iV)^m).
- Integer separation needs the normalized ladder, and the cost c_m of normalization must be paid in E. So E > 1
  needs a > 4m·c_m/3.

## Limits: where it stops, with examples
- Without normalization the ladder is not integral, already at m = 3: T₃ = 180x² + 140/3 (A3). For m = 4 no
  multiplier is needed (Ljunggren). For even m, d_m(n) suffices but is not shown to be least.
- To give a Thue-type result, the exponent 1 + κ of the measure must not exceed the degree m. For m = 3 that is
  E² ≥ B, which holds for a ≥ 1000. So a must be large: [cubic_angle_sector](../theorems/cubic_angle_sector/README.md) uses a ≥ 1000, and gets the exact minimum
  only from a ≥ 2,604,167.
- Only the sector of the small root is reached; the other roots of the form are not treated.
- The method is classical in spirit (hypergeometric or Padé); see the parent's owner's note.

## Theorems using it
[cubic_angle_sector](../theorems/cubic_angle_sector/README.md).

## Provenance
- Branch `cloud-agent/elementary-connections` (an OpenAI Codex cloud run):
  - the first-campaign A1–A3 (`cloud_agent/angle_division/CLAIM.md`);
  - the round-2 Theorems 1–2 (`cloud_agent/round2/angle_integrality/proof.md`).
  Each had one external review, PASS.
- The case m = 4 is [approx.normsign](approx.normsign.md) ([ljunggren](../theorems/ljunggren/README.md), proof 2).
- The denominator count is the one in
  [PLUS_K_PROOF.md](../../shared/prior/ghost_parametrization/team_work/plus_k_completion/PLUS_K_PROOF.md). There it is
  used for m = 3, with a 9ⁿ multiplier and a different normalization.
- No priority is claimed.
