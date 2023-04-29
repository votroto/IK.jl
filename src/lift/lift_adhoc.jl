using Base.Iterators: drop, peel, take
using DynamicPolynomials
using Gurobi
using JuMP
using MultivariatePolynomials
using SCIP

include("../utils.jl")
include("../jump_extensions.jl")
include("../denavit_hartenberg.jl")
include("../modelling.jl")

function lift_mono_adhoc!(var_map, mon)
        unfolded = replicate(variables(mon), exponents(mon))
        mapped = map(a -> get(var_map, a, a), unfolded)
        _reduce(jump_quadratic_product, mapped; init=1)
end

function lifting_vars_adhoc!(eqs, var_map)
        lifts = Dict()
        monos = unique(sort(mapreduce(monomials, vcat, eqs)))

        for mon in monos
                if !(monos in keys(lifts))
                        lifts[mon] = lift_mono_adhoc!(var_map, mon)
                end
        end

        lifts
end

function lift_adhoc(d, r, α, M, c, s)
        ids = eachindex(d)
        @polyvar pc[ids] ps[ids]

        fwd, rev = build_eqs(d, r, α, pc, ps)

        chain_poly_dirty = prod(fwd) .- M * prod(rev)
        chain_poly_clean = mapcoefficients.(round_zero, chain_poly_dirty)

        var_map = Dict([pc; ps] .=> [c; s])
        lvars = lifting_vars_adhoc!(chain_poly_clean, var_map)
        chain_jump = lift_poly.(Ref(lvars), chain_poly_clean)

        chain_jump
end