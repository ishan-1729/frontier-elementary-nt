# Results and dependency index

| Stage | Principal contents | Status boundary |
|---|---|---|
| 00 | Initial k=13 cubic trial | Retains the unresolved local cube-root-alignment gap; superseded for that case by stage 01. |
| 01 | k=13 cubic complete argument, unit certificate, continued fractions | Historical completed-case proof with Delaunay–Nagell input; exact checker rerun passed. |
| 02 | Arbitrary-k Faulhaber reductions; k=15 cubic proof | Fixed-case cubic argument completed; uniform reductions do not solve all cases. Two checkers passed. |
| 03 | Exact marked factor graphs, cycle descent, norm cycles | Algebraic graph framework and conditional transitions; no all-source closing theorem. Checker passed. |
| 04 | Kummer reduction graphs, critical groups, gonality | Graph calculations passed; arithmetic applications retain semistable-model assumptions not independently constructed by the checker. |
| 05 | Prym quotient, modular homology, matroid determinants | Conditional on the preceding arithmetic models; finite calculations passed. No automatic source divisor divisibility. |
| 06 | Cubic theta classes, 691, local-lifting countertest | Exact audit of earlier proposed routes; does not turn local witnesses into a global source. Checker passed. |
| 07 | S15 squares, S26 cubes, S35 squares; infinity profiles; cubic transports | Three full fixed-case elementary computer-assisted exclusions. Effective infinity criteria and exact transports do not classify all remaining sources. |

The archive does not establish the whole Schäffer conjecture. No claim of historical novelty is made for any fixed case or standard theorem. Original notes preserve the wording and limitations from their respective stages; this index is not an independent endorsement of every geometric inference.

## Additional completed exponent subfamilies

For positive n>1:
- S13(n) is not an m-th power whenever 3 divides m (stage 01).
- S15(n) is not an m-th power whenever 2 or 3 divides m (stages 02 and 07).
- S26(n) is not an m-th power whenever 3 divides m (stage 07).
- S35(n) is not an m-th power whenever 2 divides m (stage 07).

Each follows by reducing an m-th power to the indicated square or cube. Other prime target exponents are not excluded by these statements.

## Reproduction

The three new stage-07 checkers require only Python's standard library. Some historical symbolic checkers and the optional certificate generator require SymPy. `verify_all.py` reruns checks in temporary copies so the original archived outputs are not overwritten. Runtime logs from the present investigation are retained in `verification_reruns/` and stage 07.
