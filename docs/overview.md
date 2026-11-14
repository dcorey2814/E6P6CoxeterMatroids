# Where each result of the abstract is computed

Notebook cells are referred to by the heading they appear under.

## Section 2: Coxeter (E6,P6)-matroids

| statement | code | where |
|---|---|---|
| 2_21 has 27 vertices; edges are the skew pairs | `Gosset_lines`, `meets` | notebook 1, "The Gosset polytope and the 27 lines" |
| N(y), lozenge(y), antipodes, 45 tritangent planes (Example 2.2) | `neighbours`, `lozenge`, `antipode`, `tritangent_pairs`, `tritangent_planes` | notebook 1 |
| Theorem 2.4: matroid = criterion (iii) = tritangent exchange | `is_matroidal`, `is_E6P6_matroid`, `has_tritangent_exchange` | notebook 1, "Three descriptions of matroids" |
| Theorem 2.5: 5486 orbits of subsets, 230 orbits of matroids | `scripts/classify.jl` | `data/subset_orbit_reps.dat`, `data/matroid_orbit_reps.dat` |
| Theorem 2.5: 1 335 373 matroids, table by dimension | `orbit_size` | notebook 1, "The classification" |
| the 7-regular matroids Q1, Q2, Q3 and the 5160 triads | `valency`, `orbit_size` | notebook 1, "The graph Gamma_Q and triads" |

## Section 3: Matroidal strata

| statement | code | where |
|---|---|---|
| the cubic D and the quadrics D_y (equation (2)) | `D`, `I` | notebook 2, "Equations" |
| Theorem 3.1 (big cell) | `isolate_vars`, `verify_big_cell` | notebook 2, "The big cell" |
| Proposition 3.2: supports are matroids | | notebook 2, "Supports are matroids" |
| equation (4): the stratum in a chart | `E6P6_realization_space` | notebook 2, "The equations of a stratum" |
| exchange graph, leaf order (Lemma 3.4) | `chart_matchings`, `exchange_graph`, `leaf_order`, `good_chart` | notebook 2 |
| Theorem 3.5: only the orbit of Q1 has no leaf order; its equations are 2x2 minors | `good_chart` | notebook 2, "The main theorem" |
| Theorem 3.5: Qnr is not realizable; all other orbits are | `stratum_is_nonempty` | notebook 2; `scripts/strata_table.jl` |
| Theorem 3.5: nonempty fibre over every prime | `realizable_over_every_prime` | `scripts/strata_table.jl`, `data/strata.csv` |
| Table 1: dimensions of the strata | `stratum_dimension` | notebook 2; `data/strata.csv` |
| Remark 3.7: no stratification | `E6P6_stratum`, `test_Q2_meets_Q1bar`, `frontier_condition_holds` | notebook 2, last section |

## Section 4: Tropicalization and toric degenerations

| statement | code | where |
|---|---|---|
| starting cone from a matroid | `Q_start`, `gfan_ideal_input`, `gfan_symmetry_input` | notebook 3; `scripts/gfan_input.jl` |
| Theorem 4.1: the fan, its f-vector, simplicial, multiplicity one | gfan | `scripts/gfan_traverse.sh`, `data/gfan/cayley_plane_trop.dat` |
| lineality space | `lineality_space_V_rep` | notebook 3 |
| ray orbits and double sixes | `ray_orbits`, `same_ray` | notebook 3 |
| Theorem 4.1: matroid subdivisions, 51 to 69 cells | `induced_subdivision`, `is_matroid_subdivision` | notebook 3; `scripts/matroid_subdivisions.jl` |
| Theorem 4.2: 60 prime orbits | `in_w`, `is_prime_binomial` | `scripts/toric_degenerations.jl` |
| Theorem 4.2: 27 lattice points, volume 78, normal, h* = (1,10,28,28,10,1) | `toric_configuration_matrix`, `toric_invariants` | notebook 3; `data/toric_degenerations.csv` |
| Table 2: facets; 53 combinatorial types | `toric_invariants` | `scripts/toric_degenerations.jl` |
| Hibi degeneration and chain polytope | `order_polytope_points`, `chain_polytope_points`, `lattice_equivalence` | notebook 3 |

## Software

Only Julia, Oscar and gfan are used. Compared with the computations reported in the abstract:

* Macaulay2 (initial ideals, Gröbner bases over ZZ) is replaced by Oscar's `initial` and `groebner_basis`.
* Normaliz (lattice points, volume, h*-vector, normality) and polymake (f-vectors, combinatorial types) are replaced by Oscar's polyhedral geometry functions; `scripts/toric_degenerations.jl` calls polymake's `isomorphic` through Oscar for the combinatorial types. Oscar's backend for these functions is polymake, which may use its own bundled libraries internally.
* gfan is called from the shell instead of through Macaulay2.
