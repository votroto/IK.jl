using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using MultivariatePolynomials

include("jump_extensions.jl")
include("denavit_hartenberg.jl")

function _matrix_default_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2, "Presolve" => 2)#"BranchDir" => 1, "Heuristics" => 0, "Cuts" => 0)
end

function _split_manipulator(ids)
        mid = div(length(ids), 2)
        f, s = take(ids, mid), drop(ids, mid)
        f, reverse(collect(s))
end

function _reduce(f, xs)
        mid = div(length(xs), 2)
        h, t = take(xs, mid), drop(xs, mid)
        f(reduce(f, h), reduce(f, t))
end

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


function matrix_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
        optimizer=_matrix_default_optimizer(), init=θ)

        m = Model(optimizer)

        slim = sin_min_max.(θl, θh)
        clim = cos_min_max.(θl, θh)

        ids = eachindex(d)
        @variable(m, clim[i][1] <= c[i in ids] <= clim[i][2])
        @variable(m, slim[i][1] <= s[i in ids] <= slim[i][2])

        set_start_value.(c, cos.(init))
        set_start_value.(s, sin.(init))

        LHS, RHS = build_mat_eqs(d, r, α, c, s)
        @constraint m LHS .== M * RHS

        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        @objective m Min sum(2 .* w .* (1 .- c .* cos.(θ) .- s .* sin.(θ)))

        optimize!(m)

        stat = termination_status(m)

        sol = has_values(m) ? atan.(value.(s), value.(c)) : missing
        obj = stat == OPTIMAL ? objective_value(m) : missing

        println(θl)
        println(sol)
        println(θh)

        sol, obj, stat, solve_time(m)
end




