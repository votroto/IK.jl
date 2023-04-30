using DynamicPolynomials
using MultivariatePolynomials

function lift_mono_tree!(var_map, mon)
    unfolded = replicate(variables(mon), exponents(mon))
    mapped = map(a -> get(var_map, a, a), unfolded)
    _reduce(jump_quadratic_product, mapped; init=1)
end

function lifting_vars_tree!(eqs, var_map)
    lifts = Dict()
    monos = unique(sort(mapreduce(monomials, vcat, eqs)))

    for mon in monos
        if !(monos in keys(lifts))
            lifts[mon] = lift_mono_tree!(var_map, mon)
        end
    end

    lifts
end

"""Creates the lifted pose constraint by allowing for cancellations and pre-
computing lifting variables in a tree-like pattern."""
function lift_tree(d, r, α, M, c, s)
    ids = eachindex(d)
    @polyvar C[ids] S[ids]

    E = build_pose_constraint_poly(d, r, α, C, S, M)
    var_map = Dict([C; S] .=> [c; s])
    lvars = lifting_vars_tree!(E, var_map)
    lifted = lift_poly.(Ref(lvars), E)

    lifted
end