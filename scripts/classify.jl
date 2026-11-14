# Classification of (E6,P6)-matroids up to W(E6) (Theorem 2.5).
#
# Writes data/subset_orbit_reps.dat (the W(E6)-orbits of subsets of the 27
# lines of sizes 1..26) and data/matroid_orbit_reps.dat (the 230 orbits of
# matroids). This is the longest computation in the repository (hours).
#
#   julia --project=. scripts/classify.jl

include(joinpath(@__DIR__, "..", "src", "E6P6CoxeterMatroids.jl"))

# Suppose all orbit representatives of size k are known. Every (k+1)-subset
# containing none of them can be skipped, and a (k+1)-subset lying in an orbit
# already found is skipped as well. Only k <= 13 is needed: the orbits of
# larger subsets are obtained by taking complements.
W_Orbits_Polytopes = Dict()

to_check = []
for k in 1:13
    to_check_update = []
    for P_index in subsets(27, k)
        sort!(P_index)

        if length(to_check) > 0 && !any(is_subset(s, P_index) for s in to_check)
            continue
        end

        if any(P_index in W_Orbits_Polytopes[K] for K in keys(W_Orbits_Polytopes) if length(K) == k )
            continue
        end

        Orb = Set( on_sets(P_index, x) for x in W_in_S27);
        W_Orbits_Polytopes[P_index] = Orb;
        push!(to_check_update, P_index)
    end
    global to_check = to_check_update[:]
end

# The orbits of size k >= 14, by taking complements.
reps_leq_13 = collect(keys(W_Orbits_Polytopes))
reps_geq_14 = []
n27 = collect(1:27)
for a in reps_leq_13
    push!(reps_geq_14, sort!(setdiff(n27, a)))
end

# Record the lexicographically minimal representative of each orbit, and sort
# the orbits by size, then lexicographically.
orbit_reps = vcat(reps_leq_13, reps_geq_14)
minimal_orbit_reps = []
for a in orbit_reps
    O = Set( on_sets(a, x) for x in W_in_S27);
    push!(minimal_orbit_reps, minimum(O))
end
sort!(minimal_orbit_reps, by = x -> [length(x), x])

# The matroids among them, by the definition (edges parallel to roots),
# together with the full set of 27 lines.
matroidal_polytopes = []
for a in minimal_orbit_reps
    if is_matroidal(a)
        push!(matroidal_polytopes, a)
    end
end
push!(matroidal_polytopes, n27)

function write_sets(file, sets)
    open(joinpath(datadir, file), "w") do io
        for v in sets
            println(io, join(v, " "))
        end
    end
end

write_sets("subset_orbit_reps.dat", minimal_orbit_reps)
write_sets("matroid_orbit_reps.dat", matroidal_polytopes)
println(length(minimal_orbit_reps) + 2, " orbits of subsets (with the empty and the full set), ",
        length(matroidal_polytopes), " orbits of matroids")
