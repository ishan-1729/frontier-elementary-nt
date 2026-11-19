# Unrestricted positive offsets: quotient and solution classification

T026/R197, 2026-10-26. **Independently approved in R198; see [review record](review.md).** This package addresses ITEM-VIETA-POSITIVE-OFFSET under the approved R195 scope. It does not mark the shared register or archive goal complete.

For fixed positive integers a,b,p and any positive integer q, all positive integral quotients of (ax²+by²)/(abp xy+q) are two explicit infinite square-multiple families together with a finite computable mixed-boundary list. There is no q≤abp restriction, no coprimality requirement, and k=1 is included. The same classification is proved for rational q>0; irrational q permits no solutions.

- [Classification and completeness proof](proof.md): uniform mixed product bound, exact axis divisibility, both positive orientations, nonprimitive labels and fully evaluated q=3 and q=8 examples.
- [Every fixed-k orbit](orbits.md): bilateral matrix-word coverage, formal generating functions, finite/parabolic/hyperbolic regimes and reversible weighted norm lattice.
- [Source and Mills comparison](comparison.md): original restrictions and preserved corrections, plus a bounded primary-paper comparison that does not transfer negative-offset type attribution.
- [Source audit](source-audit.json): current unit paths, hashes, actual reading scopes, dependencies, retained primary artifacts and open boundaries.

The direct proofs reproduce every needed descent and return argument. Finite scans are implementation controls, not universal completeness evidence. Reflections, axis constructions and Chebyshev mechanisms are attributed to Ishan's source work and the reviewed Vieta dossier. No novelty or unique-orbit claim is made.

## Executable classifier

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 experiments/positive-offset-classification/classify.py 1 1 1 3
PYTHONDONTWRITEBYTECODE=1 python3 experiments/positive-offset-classification/classify.py 1 1 1 8 --k 2 --orbit-radius 3
PYTHONDONTWRITEBYTECODE=1 python3 experiments/positive-offset-classification/classify.py 2 3 1 7/2
```

[The classifier](../../experiments/positive-offset-classification/classify.py) accepts exact rational offsets and returns both labelled symbolic axis families and every finite ordered mixed record with both positive seeds. `--k` adds every compatible fixed-k seed. `--orbit-radius R` gives only the stated finite sample j=−R,…,R, preserving labels and duplicate descriptions; the complete mathematical parameter domain is all integer j. It rejects nonpositive parameters and never uses floating-point arithmetic.

For q=3 and a=b=p=1 the exact quotient set is {3t²:t≥1}∪{1,5}; the mixed records are precisely (k,u,v)=(1,1,1),(5,1,2),(5,2,1). The q=8,k=2 example keeps mixed gcd2 witnesses (6,2),(2,6), as well as other labels and axes; primitive filtering would lose valid solutions.

## Verification and limits

```sh
PYTHONDONTWRITEBYTECODE=1 python3 experiments/positive-offset-classification/verify.py
git diff --check -- research/positive-offset-classification experiments/positive-offset-classification
```

Normal verification compares [certificates.json](../../experiments/positive-offset-classification/certificates.json) and [record.json](../../experiments/positive-offset-classification/record.json) without rewriting. Explicit `--write-record` is for initial generation. It checks exact symbolic reflection/boundary/norm/matrix identities, original and retained-primary hashes, authored local links and whitespace. The arithmetic is integer or rational and the runtime needs only Python's standard library. PDF rendering is an optional source-inspection aid.

Finite controls use a,b=1,…,3,p=1,2 and reduced q=r/s with r=1,…,12,s=1,…,3. They compare every global boundary record with the independent fixed-k rectangle algorithm for k=1,…,20; independently descend every integral-quotient point with x,y=1,…,30 without a quotient cutoff; and reverse its word. All labelled fixed-k seeds are checked at j=−6,…,6 in both orientations, including signed intermediates and denominator return. Separate valuation controls use r≤24,s≤5,w≤12,z≤40. Matrix controls use a,b≤4,p≤2,k≤5,j=−12,…,12 and generating-function coefficients through degree20 in both directions. The saved record retains exact counts and examples.

R198 independently approved this precise scope without corrections; see the [review record](review.md). Nonintegral reflections, canonical positive-level representatives, parameter uniqueness and broader Primer combinations remain separate open obligations. Neither bounded source reading nor this classification certifies the whole containing project.
