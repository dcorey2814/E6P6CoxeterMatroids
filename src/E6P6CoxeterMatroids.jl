# E6/P6 Coxeter matroids
#
# Code accompanying the extended abstract
#   D. Corey, "Coxeter matroids of type (E6,P6) and the E6-Grassmannian".
#
# Usage:  include("src/E6P6CoxeterMatroids.jl")
# Everything is defined in the current module, as in an interactive session.

using Oscar
const pm = Oscar.Polymake

include("gosset.jl")        # the 27 lines, the polytope 2_21, the group W(E6)
include("matroids.jl")      # (E6,P6)-matroids and their W(E6)-orbits
include("cayley_plane.jl")  # equations and charts of E6/P6
include("strata.jl")        # matroidal strata
include("tropical.jl")      # the tropical E6-Grassmannian (gfan input/output)
include("toric.jl")         # toric degenerations
