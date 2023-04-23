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
include("local_kinematics.jl")

function chain_constraint_poly(d, r, α, M, c, s)
        ids = eachindex(d)
        @polyvar pc[ids] ps[ids]
        
        fwd, rev = build_eqs(d, r, α, pc, ps)

        chain_poly_dirty = prod(fwd) .- M * prod(rev)
        chain_poly_clean = mapcoefficients.(round_zero, chain_poly_dirty)
        chain_jump = map(e -> e([pc; ps] => [c; s]), chain_poly_clean)

        chain_jump
end

function solve_inverse_kinematics_poly(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer(), init=θ)

        m = Model(optimizer)

        ids = eachindex(d)
        @variable(m, c[ids])
        @variable(m, s[ids])
        constrain_trig_vars.(c, s, θl, θh, init)

        E = chain_constraint_poly(d, r, α, M, c, s)

        @constraint m E .== 0
        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        @objective m Min sum(2 .* w .* (1 .- c .* cos.(θ) .- s .* sin.(θ)))
        optimize!(m)

        extract_solution(c, s, m)
end