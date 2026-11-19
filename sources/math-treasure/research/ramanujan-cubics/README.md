# Ramanujan cubic cluster — R164 proof dossier

Independently approved within the documented scope by complementary reviews
R176 and R177. Task T026, worker package R164 (acceptance scope R163),
source recovery R159. See the [review record](review.md).

- [Core cubic proof](core.md): affine normalization, both discriminant
  choices, exceptional cases, compatible complex radicals, real root
  branches and the precise field-automorphism theorem.
- [Real RCP classification](rcp.md): complete zero/nonzero constant cases,
  corrected parameter symmetry, radical substitutions and seventh-angle
  families.
- [All examples](examples.md): exact cyclotomic factorizations and identities,
  minimal degrees, real radical selection and numerical-comparison bound.
- [Both quadratic proofs](quadratic.md): corrected denominator and cross-term
  errors, with positive-root restrictions.
- [Source ledger](ledger.md), [machine ledger](ledger.json) and
  [source/dependency audit](sources.md): all numbered claims and contextual
  boundaries, including the owner's separate unfinished programmes.

Reproduce from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 experiments/ramanujan-cubics/verify.py --write-record
PYTHONDONTWRITEBYTECODE=1 python3 experiments/ramanujan-cubics/verify.py
git diff --check -- research/ramanujan-cubics experiments/ramanujan-cubics
```

Python3 with SymPy is required. The verifier is inspectable exact arithmetic;
it uses no numerical tolerance or network dependency. It regenerates the
universal and example records, checks five original/snapshot hashes, ledger
coverage, local links and authored whitespace. The normal invocation compares
all mathematical outputs to the retained record. Page rendering is optional;
the retained render script used the installed Windows pypdfium2 runtime
because that library was unavailable in WSL. It is not a proof dependency.

The corrected source statements include Dresden's translated-cube constant;
Wituła's product sign, parameter symmetry, missing r=0 coverage, repeated
parameter representations and (40) scale convention; and Shevelev's two
incorrect denominators and omitted +2 cross term. Refuted printed assertions
remain identifiable in the ledger and proof text. The exact scope of each
repair, rather than blanket source correctness, is what the complementary reviews approve.
