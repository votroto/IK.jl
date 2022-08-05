using Base.Iterators: peel, drop, take
using Gurobi
using JuMP

include("denavit_hartenberg.jl")

function _default_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2)
end

function _split_manipulator(ids)
	mid = div(length(ids), 2)
        f, s = take(ids, mid), drop(ids, mid)
        f, reverse(collect(s))
end

function add_dh_chain!(m, T, ids)
        h, t = peel(T.(ids))

	X = [@variable(m, [1:4, 1:4]) for i in drop(ids, 2)]
        es = prod.(zip([h, X...], t))

        for (e, x) in zip(es, X)
                @constraint m e .== x
        end

        (length(es) == 0) ? h : last(es)
end

#=
d offset
r radius
α twist
M desired pose
θ initial angle
w angle weights
=#

function solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer())

        T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
        iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

        ids = eachindex(d)
        fwd, rev = _split_manipulator(ids)

        m = Model(optimizer)

        @variable(m, -1 <= c[ids] <= 1)
        @variable(m, -1 <= s[ids] <= 1)

        F = add_dh_chain!(m, T, fwd)
        R = add_dh_chain!(m, iT, rev)

        @constraint m F .== M * R

        @constraint m c .^ 2 + s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) - s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) - s .>= 0

        @objective m Min sum(2 * w .* (1 .- c .* cos.(θ) - s .* sin.(θ)))

        optimize!(m)

        stat = termination_status(m)
        sol = has_values(m) ? atan.(value.(s), value.(c)) : missing
        obj = stat == OPTIMAL ? objective_value(m) : missing

        sol, obj, stat, solve_time(m)
end
