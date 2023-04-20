using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using MultivariatePolynomials

include("utils.jl")
include("jump_extensions.jl")
include("denavit_hartenberg.jl")

function build_mat_eqs(d, r, α, c, s)
        T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
        iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

        ids = eachindex(d)
        fwd, rev = _split_manipulator(ids)

        F = map(T, fwd)
        R = map(iT, rev)

        LHS = _reduce(jump_quadratic_product, F)
        RHS = _reduce(jump_quadratic_product, R)

        LHS, RHS
end

function solve_inverse_kinematics_matrix(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer(), init=θ)

        m = Model(optimizer)

        ids = eachindex(d)
        @variable(m, c[ids])
        @variable(m, s[ids])
        constrain_trig_var.(c, s, θl, θh, init)

        LHS, RHS = build_mat_eqs(d, r, α, c, s)
        
        @constraint m LHS .== M * RHS
        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        @objective m Min sum(w .* lin_angdiff_approx.(c, s, θ))

        optimize!(m)

        extract_solution(c, s, m)
end




