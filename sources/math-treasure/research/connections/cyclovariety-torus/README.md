# Every cyclovariety meets the unit torus

T026/R281, independently reviewed in R282 and R285, 26 October 2026.
**The complete geometric target has a checked proof. Priority is unestablished.**
This is a local research result; no external referee or author acceptance is
claimed. The original screening document remains frozen with its historical
candidate label; this document records the subsequent review.

For every integer k≥1, the cyclovariety C_k defined by Cornelissen, Hokken and
Ringeling meets the unit torus T^(φ(k)+1). This answers the exact question
posed in §10(iii) of their [paper](https://arxiv.org/abs/2507.09303), whose
latest arXiv version is v3, 19 June 2026. The
[journal publication](https://doi.org/10.1112/plms.70198) appeared on
25 August 2026. A bounded search on 26 October found no verified earlier
resolution or equivalent theorem; it does not establish that this proof is
the first solution.

The [complete proof](../breadth-screening/cyclovariety/R281-candidate.md)
actually shows that F_k has a zero on T^φ(k) for every k>1, where F_k is the
sum of the Laurent monomials representing **all** kth roots of unity in the
integral cyclotomic basis. Since F_k(1)=k, continuity of its absolute value
then supplies a point with |F_k|=1. Taking x0=−F_k gives the required point
of C_k. The zero of F_k alone would have x0=0 and would not suffice. For
k=1 the point is explicitly (−1,1).

The key lemma concerns a nonempty finite exponent set A in a lattice L,
preserved by an integral automorphism U with U−I invertible over Q, and
with A∩2L empty. Its character sum must vanish somewhere on Hom(L,S¹).
If it never vanished, a logarithm on the covering vector space would have
an integer winding homomorphism. Invariance under U makes that homomorphism
zero. The resulting logarithm on the torus, normalized at the identity,
commutes with complex conjugation and is therefore real at every sign
character. All those character sums would be positive. Their average is
zero because no exponent lies in 2L, a contradiction.

Multiplication by ζ_k gives the required U, and its orbit of the primitive
lattice vector 1 consists of all kth roots. This checks the hypotheses for
every k>1, without a restriction on φ(k)/k. The proof uses a real intermediate
value theorem only on |F_k| and never treats a complex-valued path as real.

The derivation was found by the breadth scout, motivated by the archive's
finite-character cancellation mechanisms, and independently checked by the
depth researcher and coordinator. The posed question and cyclovariety
definitions belong to the cited authors; this argument is not attributed
to Ishan's original notes. See the [review and source-search limits](review.md).

The mathematical result concerns this geometric intersection question.
It does not resolve the paper's other open questions. Its equivalence with
non-expansiveness of the associated dynamical action has not been audited
here and is not imported into the approved result.
