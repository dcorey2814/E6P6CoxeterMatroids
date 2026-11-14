# Matroidal strata of the 230 orbits (Theorem 3.5 and Table 1).
#
# Writes data/strata.csv with one row per orbit: size of Q, dimension of
# conv(Q), dimension of the stratum (-1 if empty), number of equations in the
# chosen chart, whether the stratum is nonempty, and whether it has a point
# over the algebraic closure of every prime field.
#
#   julia --project=. scripts/strata_table.jl

include(joinpath(@__DIR__, "..", "src", "E6P6CoxeterMatroids.jl"))

open(joinpath(datadir, "strata.csv"), "w") do io
    println(io, "orbit,size,dim_polytope,dim_stratum,equations,nonempty,all_primes")
    for (i, Q) in enumerate(E6_mats)
        y = something(good_chart(Q), first(Q))
        nonempty = stratum_is_nonempty(Q)
        row = (i, length(Q), dim(convex_hull(Gosset_lines[Q, :])),
               nonempty ? stratum_dimension(Q) : -1, length(exchange_graph(Q, y)),
               nonempty, nonempty && realizable_over_every_prime(Q))
        println(io, join(row, ","))
        println("orbit $i done")
    end
end
