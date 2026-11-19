# Methods

IDs are dotted, `family.method.variant`. A child either refines its parent or generalizes it to a
new setting. Each card says what the method needs, where it stops, and which theorems use it.
To add a method or a generalization, follow [../README.md](../README.md#adding-to-the-catalogue).

```
ghost                          ghost parametrization with descent (framework)
├── ghost.rational             ghosts = primitive points of an enlarged genus-1 equation
│   └── ghost.rational.multi_index   OPEN: two or more independent indices (rank >= 2)
└── ghost.unit                 ghosts = integer points of a norm form with one unit
    └── ghost.unit.multi_unit        two or more units: reduced-cone coverage, kernel-lattice sieve
lifting                        all-index exclusion of a fixed-value marker
├── lifting.somos              exact p-adic lifting on Somos-4 sequences
├── lifting.halving            halving onto the bounded real component (no prime needed)
└── lifting.binomial           finite-binomial zero count (three-coefficient bound)
    └── lifting.binomial.multi       two-index count: coprime or saturated reductions; frame criterion
approx                         explicit rational approximation
├── approx.normsign            integer ladder with a Wronskian sign certificate
│   └── approx.normsign.all_degree   every degree m ≥ 2: the sign lemma, exact odd-degree denominators
└── approx.binomial            binomial (Padé) approximants
norm                           arithmetic of integer pairs and triples
├── norm.pairs                 Euclidean division, Bézout and cube extraction in pairs
└── norm.triples               Euclidean division, units by descent and square allocation in a cubic ring
    └── norm.triples.multiplier      division with multipliers {1,2} (Dedekind–Hasse) when not norm-Euclidean; one-unit descent
descent                        classical descents
└── descent.classical          valuation, auxiliary-multiplier, denominator and order descents
involution                     descents and counts by involutions on solution sets
├── involution.vieta           Vieta jumps as reflections in conjugate diameters; fundamental arc; gcd invariant
└── involution.parity          parity of fixed points of side involutions on hyperbolic polygons
matrix                         matrix rings
└── matrix.cayley_hamilton     Cayley–Hamilton linearisation of powers in M₂(Z)
```

| ID | Card | Theorems |
|---|---|---|
| `ghost` | [ghost.md](ghost.md) | all `ghost.*` theorems |
| `ghost.rational` | [ghost.rational.md](ghost.rational.md) | ljunggren (proof 1), mordell_1963 (proof 2), mordell_k7, quartic_t_family, x2_x_1_y3 |
| `ghost.unit` | [ghost.unit.md](ghost.unit.md) | mordell_negative_k |
| `ghost.unit.multi_unit` | [ghost.unit.multi_unit.md](ghost.unit.multi_unit.md) | triangular_tetrahedral, triangular_square_pyramidal, mordell_k_minus17 |
| `lifting.somos` | [lifting.somos.md](lifting.somos.md) | ljunggren (proof 1), mordell_k7, quartic_t_family, x2_x_1_y3 |
| `lifting.halving` | [lifting.halving.md](lifting.halving.md) | mordell_1963 (proof 2) |
| `lifting.binomial` | [lifting.binomial.md](lifting.binomial.md) | mordell_negative_k |
| `lifting.binomial.multi` | [lifting.binomial.multi.md](lifting.binomial.multi.md) | triangular_tetrahedral, triangular_square_pyramidal, mordell_k_minus17 (Proposition 17) |
| `approx.normsign` | [approx.normsign.md](approx.normsign.md) | ljunggren (proof 2), quartic_t_family (proof 2: odd Pell index of t=6 and t=12) |
| `approx.normsign.all_degree` | [approx.normsign.all_degree.md](approx.normsign.all_degree.md) | cubic_angle_sector |
| `approx.binomial` | [approx.binomial.md](approx.binomial.md) | mordell_positive_k |
| `norm.pairs` | [norm.pairs.md](norm.pairs.md) | mordell_small_k, mordell_k7, x2_x_1_y3 |
| `norm.triples` | [norm.triples.md](norm.triples.md) | mordell_1963 (proofs 1 and 2); triangular_square_pyramidal ((H)); triangular_tetrahedral, mordell_k_minus17 ((H), via its child) |
| `norm.triples.multiplier` | [norm.triples.multiplier.md](norm.triples.multiplier.md) | triangular_tetrahedral, mordell_k_minus17 (hypothesis (H)); triangular_square_pyramidal (the one-unit descent only) |
| `descent.classical` | [descent.classical.md](descent.classical.md) | radical_descent |
| `involution.vieta` | [involution.vieta.md](involution.vieta.md) | vieta_conics, vieta_gcd_invariant |
| `involution.parity` | [involution.parity.md](involution.parity.md) | parity_genus_forms |
| `matrix.cayley_hamilton` | [matrix.cayley_hamilton.md](matrix.cayley_hamilton.md) | matrix_cubes |

The open node `ghost.rational.multi_index` is a generalization target. It is described on its parent's card and gets
its own card once work on it starts. `ghost.unit.multi_unit` has had its own card since 2026-10-15.
