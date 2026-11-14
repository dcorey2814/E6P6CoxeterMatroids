# Prepare the input of gfan (Theorem 4.1).
#
# Writes data/gfan/cayley_plane_init.dat (the ideal I and the starting weight,
# the indicator vector of the matroid Q_start) and data/gfan/symmetries.dat
# (the six generators of W(E6) as signed permutations of the coordinates).
# Then run scripts/gfan_traverse.sh.
#
#   julia --project=. scripts/gfan_input.jl

include(joinpath(@__DIR__, "..", "src", "E6P6CoxeterMatroids.jl"))

w = matroid_index_to_weight_vector(Q_start)
write(joinpath(datadir, "gfan", "cayley_plane_init.dat"), gfan_ideal_input(w))
write(joinpath(datadir, "gfan", "symmetries.dat"), gfan_symmetry_input())
