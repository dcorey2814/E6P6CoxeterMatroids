# Quick checks of the main statements of the abstract (a few minutes).
#
#   julia --project=. test/runtests.jl

using Test
include(joinpath(@__DIR__, "..", "src", "E6P6CoxeterMatroids.jl"))

@testset "27 lines" begin
    @test length(W_perms) == 51840
    @test all(length(lozenge(y)) == 10 && length(neighbours(y)) == 16 for y in 1:27)
    @test length(tritangent_planes) == 45
    @test all(Set(opposite_cross_polytope(abc[y])) == Set(abc[lozenge(y)]) for y in 1:27)
end

@testset "matroids (Theorems 2.4, 2.5)" begin
    subs = read_sets(joinpath(datadir, "subset_orbit_reps.dat"))
    @test all(Q -> is_matroidal(Q) == is_E6P6_matroid(Q) == has_tritangent_exchange(Q), subs[1:50:end])
    @test filter(is_E6P6_matroid, subs) == E6_mats[1:229]
    @test sum(orbit_size, E6_mats) == 1335373
    @test minimal_rep(Qnr) == E6_mats[139]
end

@testset "charts and strata (Theorems 3.1, 3.5)" begin
    @test all(verify_big_cell, abc)
    dims = [stratum_dimension(Q) for Q in E6_mats if Q != E6_mats[139]]
    @test [count(==(d), dims) for d in 0:16] == [1, 1, 2, 3, 8, 12, 23, 34, 43, 38, 30, 17, 9, 4, 2, 1, 1]
    @test [i for i in 1:230 if isnothing(good_chart(E6_mats[i]))] == [findfirst(==(minimal_rep(Q1)), E6_mats)]
    @test !stratum_is_nonempty(Qnr)
    @test stratum_is_nonempty(Q1) && realizable_over_every_prime(Q1)
end

@testset "tropical variety (Theorem 4.1)" begin
    squeeze(s) = join(split(s))
    w = matroid_index_to_weight_vector(Q_start)
    @test squeeze(gfan_ideal_input(w)) == squeeze(read(joinpath(datadir, "gfan", "cayley_plane_init.dat"), String))
    @test gfan_symmetry_input() == read(joinpath(datadir, "gfan", "symmetries.dat"), String)
    @test [count(c -> length(c) == k, cones_OP2) for k in 1:10] == [10, 56, 224, 697, 1639, 2841, 3558, 3070, 1639, 421]
    @test sort(length.(ray_orbits())) == [27, 36, 36, 216, 240, 270, 540, 720, 4320, 6480]
    @test all(is_matroid_subdivision(vec_rep(c)) for c in max_cones_OP2[1:5])
end

@testset "toric degenerations (Theorem 4.2)" begin
    J = in_w(vec_rep(max_cones_OP2[88]))          # the Hibi degeneration
    @test is_prime_binomial(J)
    A, _ = toric_configuration_matrix(J)
    @test size(A) == (16, 27)
    @test !isnothing(lattice_equivalence(toric_points(A), order_polytope_points))
    @test lattice_volume(convex_hull(matrix(QQ, toric_points(A)))) == 78
    @test !is_prime_binomial(in_w(vec_rep(max_cones_OP2[5])))
end
