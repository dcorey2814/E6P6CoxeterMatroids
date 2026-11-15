# Data for the website

These tables are read by `docs/index.html`. They are written from the files in `data/` by

```sh
python3 scripts/dashboard_data.py
```

which needs Python 3, numpy and scipy, and takes a few minutes (it uses all cores).

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

## `matroids.json`

Read by the matroid explorer. `orbits` has one entry per orbit, in the order of `matroid_orbits.csv`, and `example` is the matroid of Figure 1 of the abstract in the chart a3. Each entry has

| key | content |
|---|---|
| `orbit`, `lines`, `size` | the orbit, its representative Q and the number of lines in Q |
| `dim_polytope`, `dim_stratum`, `nonempty_every_prime` | as in `matroid_orbits.csv`; `dim_stratum` is `null` for the empty stratum |
| `matroids_in_orbit`, `stabilizer_order` | the size of the orbit, and of the stabilizer of Q in W(E6) |
| `gamma_edges` | the edges of the graph Gamma_Q: the pairs of lines of Q that meet |
| `degree_min`, `degree_max` | the smallest and largest degree in Gamma_Q |
| `tritangent_planes` | the number of tritangent planes contained in Q |
| `kind` | `leaf`, or `Q1` (no chart admits a leaf order) or `Qnr` (the stratum is empty) |
| `chart` | the chart y: `good_chart(Q)` of `src/strata.jl`, or the first line of Q for `Q1` and `Qnr` |
| `meeting_chart`, `skew_to_chart` | the lines of Q meeting y and skew to y, in the order drawn |
| `units` | for each line v of Q meeting y, its coordinate p_v = g_v, as terms `[sign, u, w]` meaning sign * p_u * p_w |
| `equations` | the equations g_v = 0 for v in E_y, in a leaf order when `kind` is `leaf`, each with its terms and, along the leaf order, the `leaf` solved for and its `solution` p_leaf = (terms) / p_over |
| `free` | the coordinates left free by the leaf order, when `kind` is `leaf` |

In the chart y, with p_y = 1, the quadric D_{v*} is s (p_v p_y - g_v) with s = 1 or -1, for the cubic D = a^T C b + Pf(C) of Section 3 of the abstract. The script checks the number of equations and the dimension of every stratum against `data/strata.csv`, and the equations and leaves of the example against Figure 1.

## `cones.csv`

One row per orbit of cones of the tropical variety, in the order of `CONES_ORBITS` in `data/gfan/cayley_plane_trop.dat` (the lineality space itself left out), so that maximal cone orbit k is cone 13734 + k.

| column | content |
|---|---|
| `cone` | number of the orbit |
| `dim` | dimension of the cone modulo the lineality space |
| `rays` | its rays, as 0-based indices into `RAYS` of the gfan output |
| `ray_orbits` | the orbit of each ray, numbered as in `ray_orbits.csv` |
| `cones_in_orbit` | the number of cones in the orbit; for each dimension these add up to the f-vector |
| `maximal_cone_orbit` | for a maximal cone, its number in `maximal_cone_orbits.csv` |
| `cells` | the number of maximal cells of the subdivision of 2_21 it induces |

## `rays.csv`

The 12885 rays of `RAYS` in the gfan output: `ray` (0-based index), `ray_orbit` and `vector` (modulo the lineality space, coordinates a1..c56).

## `subdivisions.json` and `subdivisions/`

The subdivision of 2_21 induced by each cone in `cones.csv`, for all 14155 orbits. It is the regular subdivision for the heights w, the sum of the rays of the cone, on the upper hull (max convention). The script checks every maximal cell exactly (an integer normal vector of the lifted cell, with every other lifted vertex strictly below it), checks that it is a full-dimensional (E6,P6)-matroid polytope, and, for the maximal cones, that the volumes of the cells add up to the volume of 2_21.

`subdivisions.json` holds `chunk` (500), the list of `files` in `subdivisions/`, each with the cones `chunk` * k + 1 to `chunk` * (k + 1), and two sets of coordinates of the 27 vertices of 2_21: `gosset`, the integer coordinates used for the subdivisions, and `euclidean`, coordinates in which W(E6) acts orthogonally.

Each file in `subdivisions/` is a list of entries `[cone, cells, plane, rotation, mirrors]`:

| field | content |
|---|---|
| `cone` | the number of the cone, as in `cones.csv` |
| `cells` | the maximal cells, separated by spaces, each written `orbit.mask` in base 36: the matroid orbit of the cell (as in `matroid_orbits.csv`), and the bit mask of its vertices (bit k for line k, 0-based, a1..c56) |
| `plane` | a 6 x 2 matrix, row by row, whose columns span a plane in the `euclidean` coordinates preserved by part of the stabilizer of the cone |
| `rotation`, `mirrors` | the largest order of a rotation, and the number of mirrors, of the stabilizer acting on that plane |

The plane is chosen using only projectors and seeded random vectors, so these files come out the same with any version of numpy. The website draws the dual graph of a subdivision with each cell at the barycenter of its vertices, projected to `plane`. Its edges, the pairs of maximal cells whose common vertices span a 5-dimensional face, are computed in the page.
