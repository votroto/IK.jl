using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using MultivariatePolynomials

include("utils.jl")
include("jump_extensions.jl")
include("denavit_hartenberg.jl")
include("local_kinematics.jl")

function solve_inverse_kinematics_poly(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer(), init=θ)
        m = Model(optimizer)

        ids = eachindex(d)
        @variable(m, c[ids])
        @variable(m, s[ids])
        constrain_trig_var.(c, s, θl, θh, init)

        E, pc, ps = build_eqs(d, r, α, M)
        E = mapcoefficients.(round_zero, E)
        E = map(e -> e([pc; ps] => [c; s]), E[1:3, :])

        @constraint m E .== 0
        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        @objective m Min sum(w .* lin_angdiff_approx.(c, s, θ))

        optimize!(m)

        extract_solution(c, s, m)
end