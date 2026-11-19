# lifting.halving: exclusion by halving onto the bounded real component
Parent: lifting     Sibling: [lifting.somos](lifting.somos.md)     Children: none yet

## What it does
It finds every index at which the marker B = 1 (an integral point) occurs, with no prime and no search bound. It
applies when the rational points are Q_n for one signed index n, the doubling map acts by n ↦ 2n, and the real curve
has two components with the odd indices on the bounded one (the *egg*).
1. **Integrality survives halving.** For the doubling map D, the reduced denominator B(R) divides B(D(R)). So if
   Q_{2m} is integral, so is Q_m, and an integral Q_n with n = 2^a·m (m odd) gives integral points Q_{2^b·m} for
   every 0 ≤ b ≤ a.
2. **The egg has few integral points.** An integral point on the egg has its abscissa between two real roots, so
   only a few integers are possible. This fixes the odd part m.
3. **One non-integral double stops each chain.** If Q_{2^c·m} is not integral, no Q_{2^a·m} with a ≥ c is. For
   y² + y = x³ − x these are Q_8 (B = 5) and Q_12 (B = 29).

## Ingredients (with the file that proves each)
All in [theorems/mordell_1963/ghost/](../theorems/mordell_1963/ghost/proof.md):

| Ingredient | Proof | Check |
|---|---|---|
| B(R) divides B(D(R)), and doubling cancels at most 37² | [proof.md](../theorems/mordell_1963/ghost/proof.md) (2.4) | [01_moves.py](../theorems/mordell_1963/ghost/scripts/01_moves.py), 2.3 |
| D maps every point to the branch; S switches the components | (2.5) | (sign argument in the text) |
| Q_n is on the egg exactly when n is odd | (4.6) | [04_index.py](../theorems/mordell_1963/ghost/scripts/04_index.py) |
| The integral points on the egg; Q_8 and Q_12; the final list | (6.1)–(6.3) | [06_exclusion.py](../theorems/mordell_1963/ghost/scripts/06_exclusion.py) |

## Requirements: when it applies
- **One index and doubling.** The rational points must be Q_n for one signed index (from
  [ghost.rational](ghost.rational.md)), with D(Q_n) = Q_{2n}. Halving then needs coverage: every point of the branch
  must be a double, as in Lemma (3.7) of the proof.
- **An egg carrying the odd indices.** The cubic must have three real roots, so that the real curve has two
  components, and the generator must lie on the bounded one. Then D maps everything to the unbounded component, and
  the odd indices stay on the egg.
- **Denominators that divide.** The marker must be B = 1 in a model with integer coefficients, so that B(R) divides
  B(D(R)); for y² + y = x³ − x the proof is the identity N₄ − A⁴ = B⁴(2A² − 2AB² + B⁴) of (2.4).

## Limits: where it stops, with examples
- **No egg.** If the cubic has one real root, the real curve is connected and nothing bounds the odd indices; another
  exclusion is needed. x² + 7 = y³ and x² + x + 1 = y³ are of this kind (j = 0, one real root) and use
  [lifting.somos](lifting.somos.md).
- **Generator on the unbounded component.** Then every Q_n is on the unbounded component, and step 2 gives nothing.
- **Other markers.** The method uses that the marker is B = 1. A marker such as "a coordinate equals a constant" does
  not descend under halving in general.
- **Rank 2** has not been tried.

## Theorems using it
- [mordell_1963](../theorems/mordell_1963/README.md) (proof 2): y(y+1) = x(x+1)(x+2).

## Provenance
Classical; no novelty is claimed. See the [literature check](../theorems/mordell_1963/ghost/literature_check_ghost.md).
- *The lemma that an integral multiple forces an integral point* is credited by Bremner–Stroeker–Tzanakis (1997) to
  Lutz (1937) and Zimmer (1970); it is Ayad's Lemma A (1992) and Silverman's AEC Ex. 9.12. Ayad (1993) and
  Bremner–Stroeker–Tzanakis use it to find all integral points of rank-one curves, splitting off the power of 2.
- *The combination with the non-identity real component,* for y² + y = x³ − x and this theorem, is Silverman's AEC
  (2nd ed., 2009) Ex. 9.13, credited there to an unpublished letter of Tate to Serre. A 2016/17 Mainz exam repeats it.
- *What the card adds* is the record of the requirements and limits, and an elementary proof of each ingredient for
  this curve (the identity of (2.4), the sign argument of (2.5)).
