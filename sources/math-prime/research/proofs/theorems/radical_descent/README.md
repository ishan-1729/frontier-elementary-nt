# Classical descents: Fermat's theorems, n²∣2ⁿ+1, and representations
Status: n²∣2ⁿ+1 is LEAN; the other results are written, with an independent final review
Methods: [descent.classical](../../methods/descent.classical.md)

## Statements
- Fermat's square-area theorem; x⁴+y⁴=z² has no nontrivial solutions; the rational points of Y²=X³−X.
- Fermat's two-square theorem for primes, included as a benchmark.
- Rational-to-integral representations by X²+DY² for D=1 and 2, and for D=3 after the change of coordinates
  (X,Y)→(X−Y,2Y). D=7 shows that representations need not return to integer coordinates.
- n²∣2ⁿ+1 if and only if n∈{1,3}. This is proved in Lean as `Radical.Orders.square_novak_iff`.
  Further results: n∣2ⁿ−1, the two variable-base families n^a∣aⁿ±1, and a constructive Novák characterization.
- t₁²+⋯+tₙ² is not a sum of fewer than n squares of rational functions.

## Proof files
- **Overview:** [radical_descent_research/README.md](../../../shared/prior/radical_descent_research/README.md).
- **Report:** [report.tex](../../../shared/prior/radical_descent_research/report.tex), with a PDF next to it.
- **Candidate and novelty ledger:** [CANDIDATES.md](../../../shared/prior/radical_descent_research/CANDIDATES.md).
- **Theorem map:** [verification/THEOREM_MAP.md](../../../shared/prior/radical_descent_research/verification/THEOREM_MAP.md).
- **Provenance:** the conclusions are classical, and historical novelty is unverified.
