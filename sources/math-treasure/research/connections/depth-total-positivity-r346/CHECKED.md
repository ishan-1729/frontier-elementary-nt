# Full shifted-power total-positivity theorem — checked status

2026-10-27. The complete R346 proof has passed independent
[R362 review](R362-receipt.json) and complementary [R363 review](R363-receipt.json).
Both normal verification programs reproduced their retained records without
rewriting them. Historical author files retain their original pending-review
wording and exact hashes.

For every even n>=2, the infinite sequence ((k+t)^n) for k>=0 has a Toeplitz
matrix with all minors nonnegative exactly when t belongs to [-c_n,1+c_n],
where

    0 < c_2 < c_4 < c_6 < ... < 1 - 1/sqrt(2).

This proves the full numbered [Conjecture5.8](https://arxiv.org/html/2607.06207v1)
and strengthens its nesting and upper bound. The endpoint limit exists,
but **its value remains unproved**. Other questions in the paper remain open.

The [interval argument](interval-proof.md) proves strict parameter monotonicity
of every normalized positive maximum. The [nesting argument](nesting-proof.md)
propagates all nonouter maxima and bounds the entire outer alternating tail;
eight exact Sturm certificates cover its four finite base comparisons.

AESW is used through its checked primary statement and hypotheses; its
original1952 proof was not separately audited. Local review does not establish
priority, exhaustive novelty, publication, external acceptance or formal
verification. This package is separate from the fixed R361 register batch.
