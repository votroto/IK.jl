using DynamicPolynomials
using MultivariatePolynomials

function lifting_vars_shor!(m, eqs, var_map)
    eq_ms = unique(sort(mapreduce(monomials, vcat, eqs)))
    eq_vs = unique(sort(mapreduce(variables, vcat, eqs)))
    half_deg = div(maximum(degree, eq_ms), 2, RoundUp)

    ms = reverse(monomials(reverse(eq_vs), 0:half_deg))
    fms = filter(x -> any(m -> divides(x, m), eq_ms), ms)
    is = eachindex(fms)
    @variable m -1 <= W[is, is] <= 1 Symmetric

    for i in eachindex(fms)
        for j in eachindex(fms)
            str_i = uppercase(string(fms[i]))
            str_j = uppercase(string(fms[j]))
            set_name(W[i, j], "W[$str_i,$str_j]")
        end
    end

    vv = Dict()
    for i in eachindex(fms)
        for j in eachindex(fms)
            if i < j
                continue
            end
            if any(m -> divides(fms[i] * fms[j], m), eq_ms)
                w = get!(vv, fms[i] * fms[j], W[i, j])
                @constraint m w == W[i, 1] * W[j, 1]
            else
                JuMP.delete(m, W[i, j])
            end
        end
    end

    for v in eq_vs
        i = findfirst(x -> x == v, fms)
        @constraint m W[i, 1] == var_map[v]
    end

    @constraint m W[1, 1] == 1
    vv
end

"""Creates the lifted pose constraint a la moment relaxations."""
function lift_shor(d, r, α, M, c, s, model)
    ids = eachindex(d)
    @polyvar C[ids] S[ids]

    E = build_pose_constraint_poly(d, r, α, C, S, M)
    var_map = Dict([C; S] .=> [c; s])
    lvars = lifting_vars_shor!(model, E, var_map)
    lifted = lift_poly.(Ref(lvars), E)

    lifted
end