using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using MultivariatePolynomials

include("utils.jl")
include("jump_extensions.jl")
include("denavit_hartenberg.jl")
include("ik_modelling.jl")

function _default_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2, "Presolve" => 2, "Threads" => 4)
end

function _scip_optimizer()
        optimizer_with_attributes(SCIP.Optimizer)
end

function xlift_poly(vv, poly)
        function inner(mo)
                vv[mo]
        end
        pp = map_monomials(inner, poly)

        println()
        pp
end

function xlifting_vars!(m, eqs, var_map)
        mo = unique(sort([mon for q in eqs for mon in monomials(q)]))
        _p = sum(mo)

        ms = reverse(monomials(variables(_p), 0:div(maxdegree(_p), 2, RoundUp)))
        fms = (filter(x -> any(m -> divides(x, m), mo), ms))
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
                        if any(m -> divides(fms[i] * fms[j], m), mo)
                                w = get!(vv, fms[i] * fms[j], W[i, j])
                                @constraint m w == W[i, 1] * W[j, 1]
                        else
                                JuMP.delete(m, W[i, j])
                        end

                end
        end
        for v in variables(_p)
                i = findfirst(x -> x == v, fms)
                @constraint m W[i, 1] == var_map[v]

        end
        @constraint m W[1, 1] == 1
        #@constraint m V[1] == 1
        W, fms, vv
end

function solve_inverse_kinematics_loc(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer(), init=θ)

        m = Model(optimizer)

        ids = eachindex(d)
        @variable(m, c[ids])
        @variable(m, s[ids])
        constrain_trig_var.(c, s, θl, θh, init)

        E, pc, ps = build_eqs(d, r, α, M)
        E = mapcoefficients.(round_zero, E)

        var_map = Dict([pc; ps] .=> [c; s])
        lvars = xlifting_vars!(m, E, var_map)
        lE = xlift_poly.(Ref(lvars), E)

        @constraint m lE .== 0
        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        @objective m Min sum(w .* lin_angdiff_approx.(c, s, θ))

        optimize!(m)

        extract_solution(c, s, m)
end



