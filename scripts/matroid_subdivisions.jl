# Every cone of the tropical E6-Grassmannian induces a subdivision of 2_21
# into (E6,P6)-matroid polytopes (last statement of Theorem 4.1).
#
# Checks a representative of each of the 14155 nonempty orbits of cones, and
# counts the maximal cells of the subdivisions induced by the maximal cones.
#
#   julia --project=. scripts/matroid_subdivisions.jl

include(joinpath(@__DIR__, "..", "src", "E6P6CoxeterMatroids.jl"))

bad = [c for c in cones_OP2 if !isempty(c) && !is_matroid_subdivision(vec_rep(c))]
println(length(cones_OP2) - 1, " orbits of cones, ", length(bad), " without a matroid subdivision")

cells = [length(induced_subdivision(vec_rep(c))) for c in max_cones_OP2]
for n in sort(unique(cells))
    println(n, " maximal cells: ", count(==(n), cells), " orbits of maximal cones")
end
