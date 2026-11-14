# Toric degenerations of E6/P6 from prime cones (Section 4.2).

# --- Initial ideals and primality -------------------------------------------

nu_0 = tropical_semiring_map(QQ);

"The initial ideal in_w(I), max convention (Oscar's `initial` uses min, hence -w)."
in_w(w) = initial(I, nu_0, -w)

"The exponent differences u - v of the binomials x^u - x^v in the reduced Gröbner basis of `I`."
function exponent_differences_from_binomial_ideal(I)
    diffs = Vector{Vector{Int}}()

    for f in groebner_basis(I, complete_reduction = true)
        iszero(f) && continue

        ce = collect(coefficients_and_exponents(f))
        exps = [Int.(e) for (_, e) in ce]

        if length(exps) != 2
            error("Generator $f has $(length(exps)) terms, not 2. It is not a binomial.")
        end

        push!(diffs, exps[1] .- exps[2])
    end

    return diffs
end

"""
    toric_configuration_matrix(I)

For a homogeneous binomial ideal `I`, returns `A, D`. The rows of `D` are the
exponent differences u - v. The columns of `A` are the 27 lattice points of the
toric variety: the rows of `A` are a basis of the integer vectors orthogonal to
all u - v and with first coordinate 0, which places the first point at the
origin and makes the polytope full-dimensional (in Z^16 for a prime cone).
"""
function toric_configuration_matrix(I)
    diffs = exponent_differences_from_binomial_ideal(I)

    if isempty(diffs)
        error("No nonzero binomial generators found.")
    end

    n = length(diffs[1])

    # Matrix whose rows are exponent differences u - v.
    D = matrix(ZZ, length(diffs), n, ZZ.(vcat(diffs...)))

    # Columns of K span the right kernel of D with first coordinate 0.
    e1 = matrix(ZZ, 1, n, [Int(j == 1) for j in 1:n])
    K = kernel(vcat(D, e1), side = :right)

    return transpose(K), D
end

"The columns of the integer matrix `A`, as a list of points."
toric_points(A) = [[Int(A[i, j]) for i in 1:nrows(A)] for j in 1:ncols(A)]

"The lattice polytope of the toric variety V(I), for a homogeneous prime binomial ideal `I`."
toric_polytope_from_binomial_ideal(I) =
    convex_hull(matrix(QQ, toric_points(toric_configuration_matrix(I)[1])))

"""
    is_prime_binomial(J)

Is the monomial-free ideal `J` prime over CC? By Eisenbud--Sturmfels this holds
iff `J` is binomial, its lattice of exponent differences is saturated, and `J`
is saturated with respect to the variables.
"""
function is_prime_binomial(J)
    is_binomial(J) || return false
    S = snf(toric_configuration_matrix(J)[2])
    all(S[i, i] <= 1 for i in 1:min(nrows(S), ncols(S))) || return false
    return J == saturate_all(J, abc)
end

"""
    toric_invariants(P)

Invariants of the lattice polytope `P` computed by Oscar's polyhedral
geometry: number of facets and of lattice points, f-vector, normalized volume,
h*-vector and normality.
"""
toric_invariants(P) = (
    facets = n_facets(P),
    lattice_points = length(lattice_points(P)),
    normalized_volume = lattice_volume(P),
    normal = is_normal(P),
    h_star = collect(coefficients(h_star_polynomial(P))),
    f_vector = f_vector(P),
)

# --- The minuscule poset: order polytope and chain polytope -----------------
# (index sets only)

# Order the 27 weights by x <= y when y - x is a nonnegative combination of
# simple roots. This is a distributive lattice.
simple_roots_in_R8 = [[ai(i) - ai(i+1) for i in 1:5]; [cij(1,2) - ai(3)]]

"`leq[x,y]` is true when x <= y in the weight poset."
leq = [x == y || Gosset_lines[y, :] - Gosset_lines[x, :] in simple_roots_in_R8 for x in 1:27, y in 1:27]
for z in 1:27, x in 1:27, y in 1:27       # transitive closure
    leq[x, y] |= leq[x, z] && leq[z, y]
end

"The 16 join-irreducible weights: those covering exactly one element."
join_irreducibles = [x for x in 1:27 if
    count(z -> z != x && leq[z, x] && !any(t -> !(t in (z, x)) && leq[z, t] && leq[t, x], 1:27), 1:27) == 1]

"The join-irreducibles below `x`; these are the 27 order ideals of the poset of join-irreducibles."
irreducibles_below(x) = [p for p in join_irreducibles if leq[p, x]]

"The maximal elements of a subset `S` of the poset."
maximal_among(S) = [p for p in S if !any(q -> q != p && leq[p, q], S)]

indicator(S) = [Int(p in S) for p in join_irreducibles]

"Vertices of the order polytope: indicators of the order ideals (the Hibi toric variety)."
order_polytope_points = [indicator(irreducibles_below(x)) for x in 1:27]

"Vertices of Stanley's chain polytope: indicators of the antichains."
chain_polytope_points = [indicator(maximal_among(irreducibles_below(x))) for x in 1:27]

# --- Comparing point configurations -----------------------------------------

"""
    quad_classes(pts)

The quadratic relations among 27 lattice points: the sets of at least two
pairs {a,b} with the same sum pts[a] + pts[b].
"""
function quad_classes(pts)
    sums = Dict{Vector{Int}, Vector{Set{Int}}}()
    for a in 1:27, b in a:27
        push!(get!(sums, pts[a] + pts[b], Set{Int}[]), Set([a, b]))
    end
    return Set(Set(c) for c in values(sums) if length(c) > 1)
end

"A permutation `s` of 1:27 carrying the quadratic relations `C` onto `T`, or `nothing` (backtracking)."
function find_bijection(C, T)
    points(c) = union(c...)
    profile(X, p) = sort([length(c) for c in X if p in points(c)])
    sort([profile(C, p) for p in 1:27]) == sort([profile(T, p) for p in 1:27]) || return nothing
    s = zeros(Int, 27)
    consistent() = all(Set(Set(s[collect(m)]) for m in c) in T
                       for c in C if all(s[p] > 0 for p in points(c)))
    function extend(p)
        p > 27 && return true
        for q in setdiff(1:27, s)
            profile(C, p) == profile(T, q) || continue
            s[p] = q
            consistent() && extend(p + 1) && return true
        end
        s[p] = 0
        return false
    end
    return extend(1) ? s : nothing
end

"""
    lattice_equivalence(pts1, pts2)

A bijection `s` such that `pts1[i] -> pts2[s[i]]` extends to an isomorphism of
the affine lattices spanned by the two configurations (so the two polytopes
are unimodularly equivalent), or `nothing` if the search finds none.
"""
function lattice_equivalence(pts1, pts2)
    s = find_bijection(quad_classes(pts1), quad_classes(pts2))
    isnothing(s) && return nothing
    A1 = matrix(ZZ, [[1; p] for p in pts1])
    A2 = matrix(ZZ, [[1; p] for p in pts2[s]])
    # every affine relation among pts1 holds among the images, and the ranks agree
    return rank(A1) == rank(A2) && iszero(kernel(A1, side = :left) * A2) ? s : nothing
end
