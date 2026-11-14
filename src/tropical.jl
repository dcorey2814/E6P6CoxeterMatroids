# The tropical E6-Grassmannian (Section 4.1).
#
# The fan itself is computed outside Julia with gfan (scripts/gfan_traverse.sh);
# this file prepares gfan's input, reads its output and checks its properties.
# Max convention throughout: in_w(f) collects the terms of maximal w-weight.

# --- Input for gfan ---------------------------------------------------------

# gfan needs variable names none of which is a prefix of another.
var_names = ["x" * lpad(i, 2, "0") for i in 1:27]
R_gfan, x_gfan = polynomial_ring(QQ, var_names)
phi_gfan = hom(R, R_gfan, x_gfan)

# W(E6) acts on the coordinates by signed permutations. Here are the six
# simple reflections: the i-th entry is the image of the i-th coordinate.
s1_vec = [
    -a2, a1, a3, a4, a5, a6,
    b2, -b1, b3, b4, b5, b6,
    c12, -c23, -c24, -c25, -c26, c13, c14, c15, c16, c34, c35, c36, c45, c46, c56
]

s2_vec = [
    -c23, c13, -c12, a4, a5, a6,
    b1, b2, b3, -c56, c46, -c45,
    a3, -a2, c14, c15, c16, a1, c24, c25, c26, c34, c35, c36, b6, -b5, b4
]

s3_vec = [
    a1, -a3, a2, a4, a5, a6,
    b1, b3, -b2, b4, b5, b6,
    -c13, c12, c14, c15, c16, c23, -c34, -c35, -c36, c24, c25, c26, c45, c46, c56
]

s4_vec = [
    a1, a2, -a4, a3, a5, a6,
    b1, b2, b4, -b3, b5, b6,
    c12, c14, -c13, c15, c16, c24, -c23, c25, c26, c34, c45, c46, -c35, -c36, c56
]

s5_vec = [
    a1, a2, a3, -a5, a4, a6,
    b1, b2, b3, b5, -b4, b6,
    c12, c13, c15, -c14, c16, c23, c25, -c24, c26, c35, -c34, c36, c45, c56, -c46
]

s6_vec = [
    a1, a2, a3, a4, -a6, a5,
    b1, b2, b3, b4, b6, -b5,
    c12, c13, c14, c16, -c15, c23, c24, c26, -c25, c34, c36, -c35, c46, -c45, c56
]

signed_symmetries = [s1_vec, s2_vec, s3_vec, s4_vec, s5_vec, s6_vec]

"The signs and the underlying permutation of a signed permutation of the coordinates."
s_signs(s) = [haskey(position_abc, z) ? 1 : -1 for z in s]
s_nosign(s) = [position_abc[e * z] for (e, z) in zip(s_signs(s), s)]

"The indicator vector of the index set `Q_ind`, used as a weight vector."
function matroid_index_to_weight_vector(Q_ind)
    wQ = zeros(Int, 27)
    for i in Q_ind
        wQ[i] = 1
    end
    return wQ
end

# The matroid whose indicator vector lies in a maximal cone of the tropical
# variety; it is the starting point of the traversal.
Q_start = [1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 19, 23, 24, 25, 26]

gfan_list(v) = "(" * join(v, ", ") * ")"

"The input of `gfan_initialforms --ideal --pair`: ring, generators of I, weight vector."
gfan_ideal_input(w) = "Q[" * join(var_names, ", ") * "]\n\n{\n" *
    join(["  $f" for f in gens(phi_gfan(I))], ",\n") * "\n}\n\n" * gfan_list(w) * "\n"

"The symmetries in gfan's format (0-based permutations, then signs), to be appended to the starting cone."
gfan_symmetry_input() =
    "\n\n{\n" * join([gfan_list(s_nosign(s) .- 1) for s in signed_symmetries], ",\n") * "\n}\n" *
    "\n{\n"   * join([gfan_list(s_signs(s)) for s in signed_symmetries], ",\n") * "\n}"

# --- Output of gfan ---------------------------------------------------------

"Split a gfan/polymake fan file into its sections: name => lines (comments removed)."
function gfan_sections(path)
    sections = Dict{String, Vector{String}}()
    name = ""
    for line in readlines(path)
        line = strip(first(split(line, "#")))
        if isempty(line)
            name = ""
        elseif name == ""
            name = line
            sections[name] = String[]
        else
            push!(sections[name], line)
        end
    end
    return sections
end

ints(line) = parse.(Int, split(replace(line, r"[{}]" => " ")))

"""
    parse_tropical_file(path)

Read the output of `gfan_tropicaltraverse --symmetry`. Returns the sections of
the file, the rays (modulo lineality), and one representative of each
W(E6)-orbit of cones and of maximal cones, as lists of ray indices (1-based).
"""
function parse_tropical_file(path)
    sections = gfan_sections(path)
    rays = ints.(sections["RAYS"])
    cones = [ints(l) .+ 1 for l in sections["CONES_ORBITS"]]
    max_cones = [ints(l) .+ 1 for l in sections["MAXIMAL_CONES_ORBITS"]]
    return sections, rays, cones, max_cones
end

filename = joinpath(datadir, "gfan", "cayley_plane_trop.dat")
trop_OP2, rays_OP2, cones_OP2, max_cones_OP2 = parse_tropical_file(filename)

"A weight vector in the relative interior of the cone `c`: the sum of its rays."
vec_rep(c) = sum(rays_OP2[c])

"Basis of the orthogonal complement of the lineality space."
orth_lineality = ints.(trop_OP2["ORTH_LINEALITY_SPACE"])

"Do `u` and `v` span the same ray modulo the lineality space?"
function same_ray(u, v)
    pu, pv = ([sum(o .* z) for o in orth_lineality] for z in (u, v))
    k = findfirst(!iszero, pu)
    return !isnothing(k) && pu[k] * pv[k] > 0 && pu[k] * pv == pv[k] * pu
end

"The W(E6)-orbits of the rays, as lists of ray indices; W permutes the coordinates."
function ray_orbits()
    index = Dict(r => i for (i, r) in enumerate(rays_OP2))
    perms = [ints(l) .+ 1 for l in trop_OP2["SYMMETRY_GENERATORS"]]
    orbits, seen = Vector{Int}[], falses(length(rays_OP2))
    for i in eachindex(rays_OP2)
        seen[i] && continue
        orbit = [i]
        seen[i] = true
        for j in orbit, g in perms       # `orbit` grows while we loop over it
            k = index[rays_OP2[j][g]]
            seen[k] || (seen[k] = true; push!(orbit, k))
        end
        push!(orbits, orbit)
    end
    return orbits
end

# --- Lineality space --------------------------------------------------------

"Equations of the lineality space of Trop(V(I)), read off a reduced Gröbner basis."
function lineality_space_H_rep(I)
    Grob = groebner_basis(I, complete_reduction = true)
    W_perp = Vector{Vector{Int64}}()
    for f in Grob
        exponents_f = collect(exponents(f))
        for i in 1:length(exponents_f) - 1, j in i+1:length(exponents_f)
            push!(W_perp, exponents_f[i] - exponents_f[j])
        end
    end
    n, W_perp = rref(matrix(QQ, W_perp))
    return W_perp[1:n,:]
end

"A basis of the lineality space of Trop(V(I)), as the rows of a matrix."
function lineality_space_V_rep(I)
    H_rep = lineality_space_H_rep(I)
    V_rep = transpose(kernel(H_rep, side = :right))
    n, V_rep = rref(V_rep)
    return V_rep
end

# --- Matroid subdivisions of 2_21 -------------------------------------------

"""
    induced_subdivision(w)

The maximal cells of the regular subdivision of 2_21 induced by `w` (upper
hull, to match the max convention), as index sets.
"""
induced_subdivision(w) = [collect(c) for c in maximal_cells(subdivision_of_points(Gosset_lines, -w))]

"Is every cell of the subdivision induced by `w` an (E6,P6)-matroid polytope?"
is_matroid_subdivision(w) = all(is_E6P6_matroid, induced_subdivision(w))
