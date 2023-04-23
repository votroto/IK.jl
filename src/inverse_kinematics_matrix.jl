using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using MultivariatePolynomials

include("modelling.jl")
include("utils.jl")
include("jump_extensions.jl")
include("denavit_hartenberg.jl")

function chain_constraint_mat(d, r, α, M, c, s)
        F, R = build_eqs(d, r, α, c, s)

        LHS = _reduce(jump_quadratic_product, F)
        RHS = _reduce(jump_quadratic_product, R)

        view(LHS .- M * RHS, 1:3, :)
end

function solve_inverse_kinematics_matrix(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer(), init=θ)

        m = Model(optimizer)

        ids = eachindex(d)
        @variable(m, c[ids])
        @variable(m, s[ids])
        constrain_trig_vars.(c, s, θl, θh, init)

        E = chain_constraint_mat(d, r, α, M, c, s)
        
        @constraint m E .== 0
        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        @objective m Min sum(w .* lin_angdiff_approx.(c, s, θ))
        optimize!(m)

        extract_solution(c, s, m)
end




