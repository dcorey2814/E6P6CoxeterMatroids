# Matroidal strata of E6/P6 (Sections 3.2 and 3.3).
#
# For a matroid Q and a chart y in Q, the stratum M_Q is cut out in the torus
# with coordinates Q ∩ N(y) by equations f_v = 0 (v in E_y) and conditions
# f_v != 0 (v in lozenge(y) ∩ Q); see equation (4) of the abstract.

# --- Combinatorics of a chart: exchange graph and leaf orders ---------------
# (index sets only; no polynomials)

"""
    chart_matchings(Q, y)

For each `v` in lozenge(y), the monomials of f_v: the tritangent pairs {u,u*}
inside `Q` with u + u* = v + y. They are pairwise disjoint and lie in N(y).
"""
chart_matchings(Q, y) = Dict(v =>
    [p for p in tritangent_pairs[antipode(v, y)] if p != sort([v, y]) && issubset(p, Q)]
    for v in lozenge(y))

"The equations of the chart: the matchings of the `v` in lozenge(y) outside `Q` with f_v != 0."
function exchange_graph(Q, y)
    m = chart_matchings(Q, y)
    return [m[v] for v in lozenge(y) if !(v in Q) && !isempty(m[v])]
end

"""
    leaf_order(matchings)

An ordering of the matchings in which each one contains a vertex of degree one
in the union of itself and its successors, or `nothing` if there is none.
Along a leaf order the equations can be solved one variable at a time.
"""
function leaf_order(matchings)
    remaining = collect(eachindex(matchings))
    order = Int[]
    while !isempty(remaining)
        deg(u) = count(u in p for i in remaining for p in matchings[i])
        k = findfirst(i -> any(deg(u) == 1 for p in matchings[i] for u in p), remaining)
        isnothing(k) && return nothing
        push!(order, popat!(remaining, k))
    end
    return order
end

"A chart `y` in `Q` admitting a leaf order, with as few equations as possible; `nothing` if none."
function good_chart(Q)
    charts = [y for y in Q if !isnothing(leaf_order(exchange_graph(Q, y)))]
    return isempty(charts) ? nothing : argmin(y -> length(exchange_graph(Q, y)), charts)
end

"""
    stratum_dimension(Q)

The dimension |Q ∩ N(y)| - |E_y| of the stratum of a realizable matroid `Q`,
computed in a chart with a leaf order. Only the orbit of Q1 has no such chart;
there the three equations are the 2x2 minors of a 2x3 matrix, of codimension 2.
"""
function stratum_dimension(Q)
    y = good_chart(Q)
    isnothing(y) && return length(intersect(Q, neighbours(first(Q)))) - 2
    return length(intersect(Q, neighbours(y))) - length(exchange_graph(Q, y))
end

# --- Equations of a stratum in a chart --------------------------------------

factors(f) = collect(keys(factor(f).fac))

"""
    E6P6_realization_space(Q, y)

The stratum of `Q` (a list of coordinates) in the chart `y`, as `[IQ, SQ]`:
`IQ` is the ideal of the equations f_v, v in lozenge(y) outside `Q`, and `SQ`
lists the polynomials that must be invertible (the coordinates of Q ∩ N(y)
and the factors of the f_v with v in lozenge(y) ∩ Q).
"""
function E6P6_realization_space(Q, y = nothing)
    if isnothing(y)
        y = choose_best_chart(Q)
    end
    elim_vars = isolate_vars(y)
    Q = sort(Q)
    Qc = setdiff(abc, Q)
    hom_vars = []
    for a in abc
        if a in Q
            push!(hom_vars, a)
        else
            push!(hom_vars, R(0))
        end
    end
    phi = hom(R, R, hom_vars)
    IQ = ideal(R, [phi(elim_vars[x]) for x in intersect(Qc, collect(keys(elim_vars)))])
    SQ1 = [phi(elim_vars[x]) for x in intersect(Q, collect(keys(elim_vars)))]
    SQ2 = [x for x in setdiff(Q, collect(keys(elim_vars))) if x != y]
    SQ = unique!(vcat([factors(f) for f in vcat(SQ1, SQ2)]...))
    return [IQ, SQ]
end

"Saturate the ideal `J` with respect to each polynomial in `U`."
function saturate_all(J, U)
    for u in U
        J = saturation(J, ideal(R, [u]))
    end
    return J
end

"`Q` (index set) is realizable over the algebraic closure of QQ: its stratum is nonempty."
stratum_is_nonempty(Q) = !isone(saturate_all(E6P6_realization_space(abc[Q])...))

"""
    realizable_over_every_prime(Q)

Compute a Gröbner basis over ZZ of the stratum of `Q` (index set) in a chart,
where the units are inverted with the unused coordinate p_y as extra variable.
If it is not the unit ideal and all leading coefficients are ±1, the stratum
has a point over the algebraic closure of every prime field.
"""
function realizable_over_every_prime(Q)
    y = something(good_chart(Q), first(Q))
    IQ, SQ = E6P6_realization_space(abc[Q], abc[y])
    RZ, _ = polynomial_ring(ZZ, line_labels)
    to_ZZ(f) = map_coefficients(ZZ, f; parent = RZ)
    K = ideal(RZ, elem_type(RZ)[to_ZZ.(gens(IQ)); to_ZZ(abc[y] * prod(SQ; init = R(1))) - 1])
    G = collect(groebner_basis(K; ordering = degrevlex(RZ)))
    return all(is_unit(leading_coefficient(f; ordering = degrevlex(RZ))) for f in G) &&
           !any(is_unit, G)
end

# --- Closures of strata (Remark 3.7) ----------------------------------------

"The ideal in P^26 of the closure of the stratum of `Q` (a list of coordinates)."
E6P6_stratum(Q) = saturate_all(I + ideal(R, setdiff(abc, Q)), Q)

"Does the closure of the stratum of `Q1` meets the stratum of `Q2`?"
function test_Q2_meets_Q1bar(Q1, Q2)
    I1, I2 = E6P6_stratum.([Q1, Q2])
    return !isone(saturate_all(I1 + I2, Q2))
end

"The frontier condition for `Q2` ⊂ `Q1`: the closure of M_Q1 either misses M_Q2 or contains it."
function frontier_condition_holds(Q1, Q2)
    I1, I2 = E6P6_stratum.([Q1, Q2])
    return isone(saturate_all(I1 + I2, Q2)) || is_subset(I1, I2)
end
