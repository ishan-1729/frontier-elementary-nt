# lifting.binomial: finite-binomial zero count
Parent: lifting     Children: [lifting.binomial.multi](lifting.binomial.multi.md)

## Statement
For integers b₀, b₁, … put f(n)=Σ_{j≤n} 4^j·C(n,j)·b_j. Suppose some b_j with j∈{0,1,2} is nonzero
and has ν₂(b_j)≤s. Then f(n)=0 for at most s+2 values n≥0. This is the three-coefficient return bound.

In matrix form: ℓ(I+4K)ⁿv=0 for at most s+2 values n≥0. In
[mordell_negative_k](../theorems/mordell_negative_k/README.md) this gives at most ν₂(mΔ₀)+2 returns,
counted over all signed n.

## Proof idea
1. For a finite set of zeros, build the integer polynomial G_N with G_N(n)=N!·f(n) for 0≤n<N.
2. Remove its 2-content and reduce mod 2. Evaluations at 0, 1 and 2 bound the content, while every coefficient of
   degree above s+2 stays even.
3. Each zero gives a monic linear factor, so there are at most s+2 zeros.

No completion and no analytic interpolation is used.

## Sources
[main.tex §returns](../../shared/prior/elementary_diophantine_paper/main.tex). The Lean statements
are `Paper.ZeroTwo.sequence_zero_bound_small` and `Paper.ZeroTwo.matrix_zero_bound_small`, in
[elementary_diophantine_paper](../../shared/prior/elementary_diophantine_paper/README.md).

## Generalization to two indices (done, conditional)
[lifting.binomial.multi](lifting.binomial.multi.md) counts the common zeros of two such sums in two indices. It
is proved when the two stable reductions mod p are coprime, or after a saturation step, and has a frame criterion
giving at most 2 zeros per orbit.

## Possible generalization (not done)
The same argument with d coefficients and any prime would describe the zero sets of integer linear recurrences on residue classes. That is a
statement of Skolem–Mahler–Lech type. Hansel (1986) may already give an elementary proof of it; this is unchecked.

## Theorems using it
[mordell_negative_k](../theorems/mordell_negative_k/README.md).
