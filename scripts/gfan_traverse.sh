#!/bin/sh
# Compute the tropical E6-Grassmannian with gfan (Theorem 4.1).
# Tested with gfan 0.6.2; the traversal takes about five minutes.
#
#   sh scripts/gfan_traverse.sh
#
# Input:  data/gfan/cayley_plane_init.dat, data/gfan/symmetries.dat
#         (written by scripts/gfan_input.jl)
# Output: data/gfan/cayley_plane_cone.dat  the starting cone and the symmetries
#         data/gfan/cayley_plane_trop.dat  the tropical variety up to symmetry
set -e
cd "$(dirname "$0")/../data/gfan"

# 1. The starting cone: the initial ideal in_w(I) for the starting weight w,
#    together with a compatible reduced Groebner basis of I.
gfan_initialforms --ideal --pair < cayley_plane_init.dat > cayley_plane_cone.dat

# 2. Append the generators of W(E6) and their signs.
cat symmetries.dat >> cayley_plane_cone.dat

# 3. Traverse the tropical variety up to symmetry.
gfan_tropicaltraverse --symmetry --symsigns --nocones < cayley_plane_cone.dat > cayley_plane_trop.dat
