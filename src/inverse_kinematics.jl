using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using MultivariatePolynomials

include("utils.jl")
include("jump_extensions.jl")
include("denavit_hartenberg.jl")

function _default_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2, "Presolve" => 2, "Threads"=>4)#, "MIPFocus"=>3)#, "Cuts"=>0)
end

function _scip_optimizer()
        optimizer_with_attributes(SCIP.Optimizer)
end

function build_eqs(d, r, α, M)
        T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
        iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

        ids = eachindex(d)
        fwd, rev = _split_manipulator(ids)

        @polyvar c[ids] s[ids]

        prod(T, fwd) .- M * prod(iT, rev), c, s
end

function map_monomials(f, poly)
        sum(coefficients(poly) .* map(f, monomials(poly)), init=0)
end

function lift_poly(lift_vars, poly)
        map_monomials(e -> lift_vars[e], poly)
end

#=
d offset
r radius
α twist
M desired pose
θ initial angle
w angle weights
=#

function lift_monon!(m, var_map, mon)
        function fefe(a::Number, b)
                #println(a)
                b
        end
        function fefe(a, b)
                eee = a * b
                #println(a|>typeof, b|>typeof)
                nvar = @variable(m, lower_bound = -1, upper_bound = 1, base_name = string(eee))
                @constraint m eee == nvar
                nvar
        end

        function fefe(a::VariableRef, b::VariableRef)
                eee = a * b
                la = JuMP.lower_bound(a)
                ua = JuMP.upper_bound(a)
                lb = JuMP.lower_bound(b)
                ub = JuMP.upper_bound(b)

                stst = JuMP.start_value(a) * JuMP.start_value(b)

                ni = min(la * lb, ua * lb, la * ub, ua * ub)
                na = max(la * lb, ua * lb, la * ub, ua * ub)

                name = string(a) * string(b)
                nvar = @variable(m, lower_bound = ni, upper_bound = na, base_name = name, start = stst)
                @constraint m eee == nvar
                nvar
        end
        pre = map(a -> get(var_map, a, a), replicate(variables(mon), exponents(mon)))
        preduce(fefe, pre; init=1)
end

function lifting_vars!(m, eqs, var_map)
        nnn = Dict()
        mo = [mon for q in eqs for mon in monomials(q)]

        for mon in mo
                v = lift_monon!(m, var_map, mon)
                nnn[mon] = v
        end

        nnn
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
        @variable m -1<=W[is, is]<=1 Symmetric


        for i in eachindex(fms)
                #set_name(V[i], "y_{" * uppercase(string(fms[i])) * "}")
                for j in eachindex(fms)
                        set_name(W[i, j], "W[" * uppercase(string(fms[i]) * "," * string(fms[j])) * "]")
                end
        end

        vv = Dict()
        for i in eachindex(fms)
                for j in eachindex(fms)
                        if i < j
                                continue
                        end
                        if any(m -> divides(fms[i] * fms[j], m), mo)
                                #@constraint m W[i, j] == V[i] * V[j]
                                w = get!(vv, fms[i] * fms[j], W[i, j])
                                @constraint m w ==  W[i,1] * W[j,1]
#				if W[i, j] == w
#println(W[i,j])
#end
                        else

                                JuMP.delete(m, W[i, j])
                        end

                end
        end
@show vv
        #@constraint m W[:, 1] .== V

        for v in variables(_p)
                i = findfirst(x -> x == v, fms)
                @constraint m W[i,1] == var_map[v]

        end
        @constraint m W[1, 1] == 1
#@constraint m V[1] == 1
        W, fms, vv
end

function solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer(), init=θ)

        m = Model(optimizer)

        slim = sin_min_max.(θl, θh)
        clim = cos_min_max.(θl, θh)

        ids = eachindex(d)
        @variable(m, clim[i][1] <= c[i in ids] <= clim[i][2])
        @variable(m, slim[i][1] <= s[i in ids] <= slim[i][2])

        set_start_value.(c, cos.(init))
        set_start_value.(s, sin.(init))

        E, pc, ps = build_eqs(d, r, α, M)
        E = mapcoefficients.(c -> (abs(c) > (1e-12)) ? c : 0.0, E)

        var_map = Dict([pc; ps] .=> [c; s])
        ls = lifting_vars!(m, E, var_map)
        @constraint m lift_poly.(Ref(ls), E) .== 0


        #W, fms, vv = xlifting_vars!(m, E, var_map)
        #Xpols = xlift_poly.(Ref(vv), E)
        #@constraint m Xpols .== 0


        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        #@objective m Min sum(2 .* w .* (1 .- c .* cos.(θ) .- s .* sin.(θ)))

#println(m)
        optimize!(m)

        stat = termination_status(m)

        sol = has_values(m) ? atan.(value.(s), value.(c)) : missing
        obj = stat == OPTIMAL ? objective_value(m) : missing

        sol, obj, stat, solve_time(m)
end



