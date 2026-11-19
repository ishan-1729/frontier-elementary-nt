# Figurate programme: polygonal reduction and exact return

T026/R244, independently reviewed in R256; dependency composed in R262.

The full Dolan–Macmillan paper is reconstructed in [proof.md](proof.md):
for every integer m>=3, every n>=0 is a sum of at most m ordinary m-gonal
numbers, with at most five summands exceeding one. The original derivation
is retained as an implication from Gauss's three-triangular theorem; the
[reviewed geometric proof](../three-triangular-proof/proof.md) now supplies
that premise. [review.md](review.md) checks the exact dependency match and
records the remaining explicit standard Dirichlet import.
The proof includes the Euler integer transformation, nonnegative Pollock
window, exact N=3..7 table, complete interval filling and all zero cases.

The same dossier gives an exact sphere/plane/orthant interpretation and a
quadratic reformulation retaining congruences and nonnegative-index bounds.
Those return conditions explain why a mixed ternary form over all integers
does not automatically settle the linked triangular–square–pentagonal question.

The owner-requested descent for 8n+3 and broader extensions remain open in
[claims.json](claims.json). A geometric existence proof does not discharge
the owner's requested descent method.
See [sources.md](sources.md) and [sources.json](sources.json) for actual reading,
current paths and hashes; [audit.md](audit.md) separates finite checks from proof.

Run `python3 experiments/figurate-polygonal/verify.py` from the repository root.
Normal mode compares [verification.json](../../experiments/figurate-polygonal/verification.json)
without writing it; only `--record` initializes or regenerates that record.
