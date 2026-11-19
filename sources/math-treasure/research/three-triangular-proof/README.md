# Three triangular numbers: an audited geometric route

T026/R255 — independently approved without corrections in R258;
documentation finalized in R262. See [review.md](review.md).

[proof.md](proof.md) proves that every n>=0 is a sum of three nonnegative
triangular numbers, via Ankeny's lattice argument. It includes squarefree
reduction, an explicit auxiliary prime class, integer lattice norm and binary
norm return, and a terminating exact construction.

[arithmetic.md](arithmetic.md) supplies the required reciprocity, Minkowski and
two-square proofs. Dirichlet's progression theorem remains an explicit standard
import with actual source-reading scope in [sources.md](sources.md) and hashes
in [sources.json](sources.json). Original Ankeny OCR was read fully; the original
scan was unavailable and is not falsely claimed inspected.

This supplies the existence theorem needed by the
[Dolan–Macmillan reduction](../figurate-polygonal/proof.md); the exact
[dependency composition](../figurate-polygonal/review.md) is now recorded.
It **does not** supply the owner's requested descent proof or broader extensions.

Run `python3 experiments/three-triangular-proof/verify.py` from the repository
root. Normal mode compares [verification.json](../../experiments/three-triangular-proof/verification.json)
read-only; only `--record` initializes it. See [audit.md](audit.md) for exact
finite ranges and their distinction from the universal proof.
