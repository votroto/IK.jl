using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using MultivariatePolynomials

include("../utils.jl")
include("../jump_extensions.jl")
include("../denavit_hartenberg.jl")
include("../modelling.jl")

function lift_poly(d, r, α, M, c, s)
        ids = eachindex(d)
        @polyvar C[ids] S[ids]
        
        fwd, rev = build_eqs(d, r, α, C, S)

        chain_poly_dirty = prod(fwd) .- M * prod(rev)
        chain_poly_clean = mapcoefficients.(round_zero, chain_poly_dirty)
        chain_jump = map(e -> e([C; S] => [c; s]), chain_poly_clean)

        chain_jump
end