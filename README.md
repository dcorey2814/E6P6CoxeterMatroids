# Coxeter matroids of type (E6,P6) and the E6-Grassmannian

Code and data for the extended abstract

> D. Corey, *Coxeter matroids of type $(E_6,P_6)$ and the $E_6$-Grassmannian* (`docs/corey_e6p6.pdf`).

The abstract has three parts, and so does this repository:

| | abstract | notebook | source |
|---|---|---|---|
| 1 | combinatorics and classification of $(E_6,P_6)$-matroids (Section 2) | `notebooks/1_matroids.ipynb` | `src/gosset.jl`, `src/matroids.jl` |
| 2 | matroidal strata of the $E_6$-Grassmannian (Section 3) | `notebooks/2_strata.ipynb` | `src/cayley_plane.jl`, `src/strata.jl` |
| 3 | tropicalization and toric degenerations (Section 4) | `notebooks/3_tropical_toric.ipynb` | `src/tropical.jl`, `src/toric.jl` |

`docs/overview.md` lists, for each statement of the abstract, the function, notebook cell or script that checks it.

## Requirements

* [Julia](https://julialang.org) 1.10 or later with [Oscar](https://www.oscar-system.org) 1.7 (and IJulia for the notebooks).
* [gfan](https://users-math.au.dk/jensen/software/gfan/gfan.html) 0.6.2, only to recompute the tropical variety.

Nothing else is used. Polyhedral computations go through Oscar (whose backend is polymake).

## Getting started

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl        # quick checks of the main statements
```

In a Julia session or a notebook, everything is loaded with

```julia
include("src/E6P6CoxeterMatroids.jl")
```

This defines all functions and data in the current session (it is a plain script, not a package). The $27$ lines are ordered $a_1,\dots,a_6,b_1,\dots,b_6,c_{12},\dots,c_{56}$, and a set of lines is a sorted vector of indices in `1:27`; `line_labels` and `abc` translate indices to names and to coordinates.

## Layout

```
src/         the code, one file per topic (see the table above)
notebooks/   three notebooks explaining the three parts of the abstract
scripts/     the long computations, which write their results to data/
data/        orbit representatives, gfan input and output, tables (see data/README.md)
docs/        the extended abstract and an overview of where each result is computed
test/        quick checks
```

## Scripts

| script | computes | time |
|---|---|---|
| `scripts/classify.jl` | orbits of subsets of the 27 lines, and of matroids (Theorem 2.5) | hours |
| `scripts/strata_table.jl` | dimensions and realizability of the 230 strata (Theorem 3.5, Table 1) | minutes |
| `scripts/gfan_input.jl`, `scripts/gfan_traverse.sh` | the tropical variety, with gfan (Theorem 4.1) | 5 minutes |
| `scripts/matroid_subdivisions.jl` | matroid subdivisions for all 14155 orbits of cones (Theorem 4.1) | long |
| `scripts/toric_degenerations.jl` | the 60 prime cones and their polytopes (Theorem 4.2, Table 2) | long |

Each script overwrites its output in `data/`, so `git diff data/` shows whether a run reproduces the stored results.

## Status

The code was reorganized from the notebooks and scripts used for the abstract. What has been run on the code in this repository, and what has not:

* The gfan pipeline was rerun with gfan 0.6.2 and reproduces `data/gfan/cayley_plane_trop.dat` exactly.
* The functions that only use index sets (matroid tests, orbits, exchange graphs, leaf orders, stratum dimensions, the gfan file parser, ray orbits, the poset and the comparison of point configurations) were run in plain Julia and reproduce the numbers of the abstract.
* **The functions that call Oscar have not been run in this form**, and the notebooks are stored without outputs. Run `test/runtests.jl` first.
* `data/strata.csv` (realizability columns) and `data/toric_degenerations.csv` hold the values obtained for the abstract with Macaulay2, Normaliz and polymake. The scripts recompute them with Oscar only.
