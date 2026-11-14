# Data for the website

These tables are read by `docs/index.html`. They are written from the files in `data/` by

```sh
python3 scripts/dashboard_data.py
```

which needs Python 3 and numpy.

Lines are named a1..a6, b1..b6, c12..c56 as in `data/README.md`. Orbits are W(E6)-orbits. W(E6) is generated from the six symmetries in `data/gfan/cayley_plane_trop.dat`.

## `matroid_orbits.csv`

One row per orbit of (E6,P6)-matroids, in the order of `data/matroid_orbit_reps.dat`.

| column | content |
|---|---|
| `orbit` | number of the orbit, as in `data/strata.csv` |
| `lines` | number of lines in Q |
| `dim_polytope` | dimension of conv(Q) |
| `dim_stratum` | dimension of the stratum, empty if the stratum is empty |
| `equations` | number of equations in the chosen chart |
| `realizable` | whether the stratum is nonempty |
| `nonempty_every_prime` | whether the stratum has a point over the algebraic closure of every prime field |
| `matroids_in_orbit` | number of matroids in the orbit |
| `representative` | the representative from `data/matroid_orbit_reps.dat`, by line names |

## `tropical_fvector.csv`

| column | content |
|---|---|
| `dim_mod_lineality` | dimension of the cone modulo the lineality space |
| `dim_projective` | dimension of the cone in R^27/R1 |
| `cones` | number of cones (the f-vector) |
| `cone_orbits` | number of orbits of cones |

## `ray_orbits.csv`

| column | content |
|---|---|
| `ray_orbit` | number of the orbit, by increasing size |
| `rays` | number of rays in the orbit |
| `description` | for three orbits, the set of lines whose indicator vector spans a ray of the orbit modulo the lineality space |

## `maximal_cone_orbits.csv`

One row per orbit of maximal cones, numbered from 1 in the order of `MAXIMAL_CONES_ORBITS` in the gfan output (`data/binomial_initial_ideals.txt` and `data/toric_summary.json` number them from 0).

| column | content |
|---|---|
| `orbit` | number of the orbit |
| `binomial_initial_ideal` | whether in_w(I) is binomial |
| `prime` | whether in_w(I) is prime |
| `known_polytope` | `order` or `chain` if the toric polytope is the order or chain polytope of the minuscule poset |

## `toric_degenerations.csv`

One row per prime orbit: the columns of `data/toric_degenerations.csv`, with

| column | content |
|---|---|
| `edges` | number of edges of the toric polytope |
| `quadric_generators`, `cubic_generators` | number of generators of in_w(I) of degree 2 and 3 |
| `known_polytope` | `Order polytope (Hibi)` or `Chain polytope`, if any |
