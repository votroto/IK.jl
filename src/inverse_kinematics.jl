using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP

include("denavit_hartenberg.jl")

function _default_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2, "Presolve"=>2)
end


function _scip_optimizer()
        optimizer_with_attributes(SCIP.Optimizer)
end

function _split_manipulator(ids)
        mid = div(length(ids), 2)
        f, s = take(ids, mid), drop(ids, mid)
        f, reverse(collect(s))
end

function build_eqs(d, r, α, M)
        T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
        iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

        ids = eachindex(d)
        fwd, rev = _split_manipulator(ids)

        @polyvar c[ids] s[ids]

        prod(T, fwd) .- M * prod(iT, rev), c, s
end


function lift_mono!(m, var_map, powers, mmm)
        temp = 1
        for (v, p) in powers
                for j in 1:p
                        eee = temp * var_map[v]
                        nvar = get!(mmm, eee, @variable(m))
                        @constraint m eee == nvar
                        temp = nvar
                end
        end
        temp
end


function lift_poly!(m, var_map, poly)
        poly = mapcoefficients(c->round(c, digits=12), poly)
        mmm = Dict()
        es = powers.(monomials(poly))
        cs = coefficients(poly)
        ns = map(e -> lift_mono!(m, var_map, e,mmm), es)

        sum(cs .* ns, init=0)
end

function cos_min_max(l, h)
        l, h = sort([l, h])
        rl, rh = mod(l, 2 * pi), mod(h, 2 * pi)
        rh = (rh > rl) ? rh : rh + 2 * pi

        inflcos = filter(x -> rl <= x <= rh, [pi * i / 2 for i in 0:4])
        infl = map(cos, [inflcos; rl; rh])
        minimum(infl), maximum(infl)
end

function sin_min_max(l, h)
        cos_min_max(l - pi / 2, h - pi / 2)
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

        m = Model(optimizer)

        slim = sin_min_max.(θl, θh)
        clim = cos_min_max.(θl, θh)

        ids = eachindex(d)
        @variable(m, clim[i][1] <= c[i in ids] <= clim[i][2])
        @variable(m, slim[i][1] <= s[i in ids] <= slim[i][2])

        E, pc, ps = build_eqs(d, r, α, M)
        var_map = Dict([pc; ps] .=> [c; s])

        @constraint m lift_poly!.(m, Ref(var_map), E) .== 0

        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        @objective m Min sum(2 .* w .* (1 .- c .* cos.(θ) .- s .* sin.(θ)))

        optimize!(m)

        stat = termination_status(m)

        sol = has_values(m) ? atan.(value.(s), value.(c)) : missing
        obj = stat == OPTIMAL ? objective_value(m) : missing

        sol, obj, stat, solve_time(m)
end