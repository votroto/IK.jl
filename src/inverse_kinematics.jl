using Base.Iterators: drop, peel, take
using DynamicPolynomials
using Gurobi
using JuMP
using MultivariatePolynomials
using SCIP

include("utils.jl")
include("jump_extensions.jl")
include("denavit_hartenberg.jl")
include("modelling.jl")

function map_monomials(f, poly)
        sum(coefficients(poly) .* map(f, monomials(poly)), init=0)
end

function lift_poly(lift_vars, poly)
        map_monomials(e -> lift_vars[e], poly)
end

function lift_monon!(m, var_map, mon)
        function inner(a::Number, b)
                b
        end
        function inner(a, b)
                eee = a * b
                nvar = @variable(m, lower_bound = -1, upper_bound = 1, base_name = string(eee))
                @constraint m eee == nvar
                nvar
        end

        function inner(a::VariableRef, b::VariableRef)
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
        _reduce(inner, pre; init=1)
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

function chain_constraint_standard(d, r, α, M, c, s, m)
        ids = eachindex(d)
        @polyvar pc[ids] ps[ids]
        
        fwd, rev = build_eqs(d, r, α, pc, ps)

        chain_poly_dirty = prod(fwd) .- M * prod(rev)
        chain_poly_clean = mapcoefficients.(round_zero, chain_poly_dirty)

        var_map = Dict([pc; ps] .=> [c; s])
        lvars = lifting_vars!(m, chain_poly_clean, var_map)
        chain_jump = lift_poly.(Ref(lvars), chain_poly_clean)
        
        chain_jump
end

function solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer(), init=θ)

        m = Model(optimizer)
        
        ids = eachindex(d)
        @variable(m, c[ids])
        @variable(m, s[ids])
        constrain_trig_vars.(c, s, θl, θh, init)

        E = chain_constraint_standard(d, r, α, M, c, s, m)

        @constraint m E .== 0
        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        @objective m Min sum(w .* lin_angdiff_approx.(c, s, θ))

        optimize!(m)

        extract_solution(c, s, m)
end