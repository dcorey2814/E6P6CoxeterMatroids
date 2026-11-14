# The Gosset polytope 2_21, the 27 lines, and the Weyl group W(E6).
# (Section 2 of the abstract.)
#
# The 27 vertices of 2_21 are always listed in the order
#     a1,...,a6, b1,...,b6, c12, c13, ..., c56,
# and a subset of them is stored as a sorted vector of indices in 1:27.

# --- The 27 lines -----------------------------------------------------------

"Names of the 27 lines, in the order used everywhere in this repository."
line_labels = vcat(["a$i" for i in 1:6], ["b$i" for i in 1:6],
                   ["c$i$j" for i in 1:5 for j in i+1:6])

# Coordinates of the vertices of 2_21 in R^8:
#   a_i = 2e_i + 2e_7,   b_i = 2e_i + 2e_8,   c_ij = (e_1 + ... + e_8) - 2e_i - 2e_j.

function ai(i)
    v = zeros(QQ, 8)
    v[i] = v[7] = 2
    return v
end

function bi(i)
    v = zeros(QQ, 8)
    v[i] = v[8] = 2
    return v
end

function cij(i,j)
    @assert i != j
    v = fill(QQ(1), 8)
    v[i] = v[j] = -1
    return v
end

"The 27 x 8 matrix whose rows are the vertices of the Gosset polytope 2_21."
function Gosset_lines_creation()
    Gosset_lines = zeros(QQ, 27, 8)
    r = 1
    for i in 1:6
        Gosset_lines[r, :] = ai(i)
        r += 1
    end

    for i in 1:6
        Gosset_lines[r, :] = bi(i)
        r += 1
    end

    for i in 1:5, j in i+1:6
        Gosset_lines[r, :] = cij(i,j)
        r += 1
    end
    return Gosset_lines
end

Gosset_lines = Gosset_lines_creation()

# --- Roots and the edge criterion (Gelfand--Serganova) ----------------------

# The 36 positive roots of E6 in the same R^8 coordinates.
positive_roots_in_R8 = []
for i in 1:5, j in i+1:6
    v = zeros(QQ, 8)
    v[i], v[j] = (1, -1)
    push!(positive_roots_in_R8, v)
end
push!(positive_roots_in_R8, [0,0,0,0,0,0,QQ(1),QQ(-1)])

for I in subsets(8,4)
    if length(intersect(I,[7,8])) != 1 || !(1 in I)
        continue
    end
    v = zeros(QQ, 8)
    J = setdiff(1:8, I)
    v[I] = fill(1//2, 4)
    v[J] = fill(-1//2, 4)
    push!(positive_roots_in_R8, v)
end

primitive_vector(v) = [QQ(b) for b in pm.common.primitive(v)]
primitive_roots_in_R8 = primitive_vector.(positive_roots_in_R8)

"The edge directions of the polytope `P`."
function edge_directions(P)
    output = []
    for b in faces(P,1)
        s, t = vertices(b)
    push!(output, t-s)
    end
    return output
end

"""
    is_matroidal(P_index)

Decide whether the lines indexed by `P_index` form an (E6,P6)-matroid, straight
from the definition: every edge of their convex hull is parallel to a root.
The much faster test `is_E6P6_matroid` (src/matroids.jl) is equivalent.
"""
function is_matroidal(P_index)
    P_verts_QQ = Gosset_lines[P_index, :]
    P = convex_hull(P_verts_QQ)
    ed = [primitive_vector(v) for v in edge_directions(P)]

    for v in ed
        if v in primitive_roots_in_R8 || -v in primitive_roots_in_R8
            continue
        else
            return false
        end
    end
    return true
end

# --- The Weyl group W(E6) acting on the 27 lines ----------------------------

E6 = root_system(:E, 6)
W = weyl_group(E6)
(w1, w2, w3, w4, w5, w6) = [fundamental_weight(E6, i) for i in 1:6]

# The 27 lines as weights of the minuscule representation (the W-orbit of w1).
lines_to_weights = Dict(
    "a1" => w1 - w2,
    "a2" => -w1 - w2 + w3,
    "a3" => -w2 - w3 + w4,
    "a4" => -w4 + w5,
    "a5" => -w5 + w6,
    "a6" => -w6,
    "b1" => w1,
    "b2" => -w1 + w3,
    "b3" => -w3 + w4,
    "b4" => w2 - w4 + w5,
    "b5" => w2 - w5 + w6,
    "b6" => w2 - w6,
    "c12" => w2 - w3,
    "c13" => -w1 + w2 + w3 - w4,
    "c14" => -w1 + w4 - w5,
    "c15" => -w1 + w5 - w6,
    "c16" => -w1 + w6,
    "c23" => w1 + w2 - w4,
    "c24" => w1 - w3 + w4 - w5,
    "c25" => w1 - w3 + w5 - w6,
    "c26" => w1 - w3 + w6,
    "c34" => w3 - w5,
    "c35" => w3 - w4 + w5 - w6,
    "c36" => w3 - w4 + w6,
    "c45" => -w2 + w4 - w6,
    "c46" => -w2 + w4 - w5 + w6,
    "c56" => -w2 + w5
)

weights_to_lines = Dict(lines_to_weights[a] => a for a in keys(lines_to_weights))

# Gosset 2_21 polytope via weights w1, ..., w6, in the order of `line_labels`.
Gosset_weights = [lines_to_weights[a] for a in line_labels]

# The action of W(E6) on the 27 lines, as permutations in S_27.
pos = Dict(Gosset_weights[i] => i for i in eachindex(Gosset_weights))

W_on_27(i, x) = pos[Gosset_weights[i] * x]

"All 51840 elements of W(E6) as plain vectors: `g[i]` is the image of line `i`."
W_perms = [[W_on_27(i, x) for i in 1:27] for x in W];

"The same 51840 elements as permutations in S_27."
W_in_S27 = perm.(W_perms);
