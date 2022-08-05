using Base.Iterators: peel, drop, take
using Gurobi
using JuMP

include("denavit_hartenberg.jl")

function _default_optimizer()
    optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2)
end

#=
d offset
r radius
α twist
M desired pose
θ initial angle
w angle weights
=#

function _split_at(xs, n)
    take(xs, n), drop(xs, n)
end

function _split_manipulator(ids)
    f, s = _split_at(ids, div(length(ids), 2, RoundUp))
    f, reverse(collect(s))
end

function add_dh_chain!(m, T, ids)
    h, t = peel(T.(ids))

    X = [@variable(m, [1:4, 1:4]) for i in drop(ids, 2)]
    es = prod.(zip([h, X...], t))

    for (e,x) in zip(es, X)
        @constraint m e .== x
    end

    last(es)
end

function solve_inverse_kinematics(d, r, α, M, θ, w; optimizer=_default_optimizer())
    T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
    iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

    ids = eachindex(d)
    fwd, rev = _split_manipulator(ids)

    m = Model(optimizer)

    @variable(m, -1 <= c[i in ids] <= 1)
    @variable(m, -1 <= s[i in ids] <= 1)

    F = add_dh_chain!(m, T, fwd)
    R = add_dh_chain!(m, iT, rev)

    @constraint m F .== M * R

    @constraint m c .^ 2 + s .^ 2 .== 1
    @constraint m (c .+ 1) .* tan.(θl ./ 2) .- s .<= 0
    @constraint m (c .+ 1) .* tan.(θh ./ 2) .- s .>= 0

    @objective m Min sum(2 .* w .* (1 .- c .* cos.(θ) .- s .* sin.(θ)))

    optimize!(m)

    stat = termination_status(m)
    sol = has_values(m) ? atan.(value.(s), value.(c)) : missing
    obj = stat == OPTIMAL ? objective_value(m) : missing

    sol, obj, stat, solve_time(m)
end