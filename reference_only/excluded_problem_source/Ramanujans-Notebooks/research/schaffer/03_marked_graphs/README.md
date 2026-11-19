# Schäffer marked-graph research bundle

Read `Schaffer_marked_graphs.md` for proofs and scope.

This is NOT a proof of the full conjecture. It contains an all-k/all-m finite signed-coefficient cover with a unitary factorization of the source root; graph gain and divisibility lemmas; a norm-cycle integrality obstruction; and a complete graph-descent control on the known k=5 square family.

Run with Python 3.10+ and SymPy:

```bash
python verify_graph_framework.py
```

`verification_output.txt` is the actual saved output. `diagnostics.json` contains the normalization data, extracted source states and controls. The default checks are finite diagnostics, not substitutes for the general proofs or a nonexistent all-source termination certificate.
