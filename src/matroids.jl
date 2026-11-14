# (E6,P6)-matroids: the dictionary with lines on a cubic surface, the
# tritangent criterion (Theorem 2.4) and W(E6)-orbits (Theorem 2.5).
#
# Everything here is plain combinatorics on index sets in 1:27. It needs only
# `Gosset_lines` and `W_perms` from src/gosset.jl.

# --- Lines, cross polytopes, tritangent planes ------------------------------

"""
`meets[i,j]` is true when the lines `i != j` meet, i.e. <x_i,x_j> = -2/3.
In the R^8 coordinates this is the vanishing of the dot product
(skew lines have dot product 4, and each line has norm 8).
"""
meets = [i != j && iszero(sum(Gosset_lines[i, :] .* Gosset_lines[j, :]))
        for i in 1:27, j in 1:27]

"The 10 lines meeting `y`: the vertices of the cross-polytope facet opposite `y`."
lozenge(y) = findall(meets[:, y])

"The 16 lines skew to `y`: its neighbours N(y) in 2_21."
neighbours(y) = [x for x in 1:27 if x != y && !meets[x, y]]

"The antipode x* of `x` in lozenge(y): the third line of the tritangent plane through x, y."
antipode(x, y) = only(z for z in lozenge(y) if meets[z, x])

"The 5 antipodal pairs {x, x*} of lozenge(y), i.e. the tritangent planes through `y`."
tritangent_pairs = [unique(sort([x, antipode(x, y)]) for x in lozenge(y)) for y in 1:27]

"The 45 tritangent planes, as sorted triples of lines."
tritangent_planes = unique(sort([p; y]) for y in 1:27 for p in tritangent_pairs[y])

# --- Matroid criteria (Theorem 2.4) -----------------------------------------

"""
    is_E6P6_matroid(Q)

Criterion (iii): for no line `y` does `Q` contain exactly one antipodal pair
of lozenge(y).
"""
is_E6P6_matroid(Q) = !isempty(Q) &&
    all(count(p -> issubset(p, Q), tritangent_pairs[y]) != 1 for y in 1:27)

"""
    has_tritangent_exchange(Q)

Criterion (iv): whenever `a, b` in `Q` meet, there is a second meeting pair
`c, d` in `Q` with a + b = c + d.
"""
function has_tritangent_exchange(Q)
    pt(i) = Gosset_lines[i, :]
    meeting = [(a, b) for a in Q for b in Q if a < b && meets[a, b]]
    return !isempty(Q) && all(any(p != q && pt(p[1]) + pt(p[2]) == pt(q[1]) + pt(q[2])
                                  for q in meeting) for p in meeting)
end

"The degree of `q` in the graph Gamma_Q joining meeting lines of `Q`."
valency(Q, q) = count(x -> meets[x, q], Q)

# --- W(E6)-orbits (Theorem 2.5) ---------------------------------------------

"The image of the set `Q` under the permutation `g`."
apply_perm(Q, g) = sort!(g[Q])

"The W(E6)-orbit of `Q`."
W_orbit(Q) = Set(apply_perm(Q, g) for g in W_perms)

"The lexicographically smallest element of the W(E6)-orbit of `Q`."
minimal_rep(Q) = minimum(apply_perm(Q, g) for g in W_perms)

"The order of the stabilizer of `Q` in W(E6)."
stabilizer_order(Q) = count(g -> apply_perm(Q, g) == Q, W_perms)

"The number of sets in the W(E6)-orbit of `Q`."
orbit_size(Q) = div(length(W_perms), stabilizer_order(Q))

"The stabilizer of `Q` in W(E6), as a permutation group on 27 points."
E6P6_automorphism_group(Q) =
    permutation_group(27, [x for x in W_in_S27 if on_sets(Q, x) == Q])

# --- Data -------------------------------------------------------------------

"Read a list of index sets, one per line."
read_sets(file) = [parse.(Int, split(line)) for line in readlines(file)]

datadir = joinpath(@__DIR__, "..", "data")

"Representatives of the 230 W(E6)-orbits of (E6,P6)-matroids (lex-minimal, sorted by size)."
E6_mats = read_sets(joinpath(datadir, "matroid_orbit_reps.dat"))

"The index set of the lines with the given names."
line_indices(names...) = sort([findfirst(==(s), line_labels) for s in names])

# Two special matroids of the abstract: Q1 is the complement of a triad, and
# Qnr = Q1 + c14 is the matroid that is not realizable over any field.
Q1  = line_indices("a1","a2","a3","a4","a5","a6","b1","b2","b3","b4","b5","b6",
            "c12","c13","c23","c45","c46","c56")
Qnr = sort([Q1; line_indices("c14")])
