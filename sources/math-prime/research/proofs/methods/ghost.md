# ghost: ghost parametrization with descent
Parent: none     Children: [ghost.rational](ghost.rational.md), [ghost.unit](ghost.unit.md)

## What it does
The equation is enlarged to a family of auxiliary records, the ghosts, that contains every
original solution. The method then:
1. Records a ghost property that every allowed move preserves: a fixed norm, a homogeneous
   equation, or primitive denominators.
2. Proves an exhaustive descent. Outside finitely many seeds, every ghost has a parent of strictly smaller
   positive integer height. Reversing the moves puts every ghost on finitely many trees.
3. Turns the trees into explicit indices: orbits of one move, or one signed index.
4. Recognises the original solutions by a fixed-value marker: a coordinate equal to a
   constant, or a reduced denominator equal to 1.
5. Excludes the marker at every index except the known ones (`lifting.*`), then reconstructs the
   solutions.

Coverage, marker exactness and reconstruction are three separate obligations
([main.tex, §1.1](../../shared/prior/elementary_diophantine_paper/main.tex)).

## Requirements
- **A terminating descent:** the enlargement must admit a descent in a positive integer height.
- **A fixed-value marker:** a marker such as "is a square" can lift to every prime power, so finite
  congruences cannot exclude it. See "Why these are hypotheses" in the
  [September target note](../../shared/prior/ghost_parametrization/ghost_parametrization_handover/inputs/ghost_parametrization_targets.md).
- **Few indices:** the trees must reduce to a bounded number of indices. That number is one for `ghost.rational`, in the
  rank-1 case, and one for `ghost.unit`, when there is a single unit.

## Theorems using it
See the children.

## Sources
- The framework: [main.tex, §1.1](../../shared/prior/elementary_diophantine_paper/main.tex), and §1 of the
  [ghost-descent Ljunggren proof](../theorems/ljunggren/ghost_descent/ljunggren_elementary_proof_descent.pdf).
- The September 2026 programme that developed it:
  [ghost_parametrization/](../../shared/prior/ghost_parametrization/README.md).
