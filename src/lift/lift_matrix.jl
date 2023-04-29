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

function lift_matrix(d, r, α, M, c, s)
        F, R = build_eqs(d, r, α, c, s)

        LHS = _reduce(jump_quadratic_product, F)
        RHS = _reduce(jump_quadratic_product, R)

        view(LHS .- M * RHS, 1:3, :)
end