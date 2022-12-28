using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using MultivariatePolynomials

include("jump_extensions.jl")
include("denavit_hartenberg.jl")

function _default_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2, "Presolve" => 2)
end

function _matrix_default_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2, "Presolve" => 2)#"BranchDir" => 1, "Heuristics" => 0, "Cuts" => 0)
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


function build_eqs_simple(d, r, α, c, s)
        T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])

        ids = eachindex(d)



        prod(T, ids), c, s
end


function map_monomials(f, poly)
        sum(coefficients(poly) .* map(f, monomials(poly)), init=0)
end

function lift_poly(lift_vars, poly)
        map_monomials(e -> lift_vars[e], poly)
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

function replicate(xs, ns)
        nxs = similar(xs, sum(ns))
        ni = 1
        for i in eachindex(xs)
                for j in 1:ns[i]
                        nxs[ni] = xs[i]
                        ni += 1
                end
        end
        nxs
end

function preduce(f, xs; init=1)
        mid = div(length(xs), 2)
        h, t = take(xs, mid), drop(xs, mid)
        f(reduce(f, h; init), reduce(f, t; init))
end


function _reduce(f, xs)
        mid = div(length(xs), 2)
        h, t = take(xs, mid), drop(xs, mid)
        f(reduce(f, h), reduce(f, t))
end

function lift_monon!(m, var_map, mon)
        function fefe(a::Number, b)
                #println(a)
                b
        end
        function fefe(a, b)
                eee = a * b
                #println(a|>typeof, b|>typeof)
                nvar = @variable(m, lower_bound = -1, upper_bound = 1, base_name = string(eee))
                @constraint m eee == nvar
                nvar
        end

        function fefe(a::VariableRef, b::VariableRef)

                eee = a * b
                la = JuMP.lower_bound(a)
                ua = JuMP.upper_bound(a)
                lb = JuMP.lower_bound(b)
                ub = JuMP.upper_bound(b)

                stst = JuMP.start_value(a) * JuMP.start_value(b)

                ni = min(la * lb, ua * lb, la * ub, ua * ub)
                na = max(la * lb, ua * lb, la * ub, ua * ub)
                #println(ni, " ",na)
                nvar = @variable(m, lower_bound = ni, upper_bound = na, base_name = string(eee), start = stst)
                @constraint m eee == nvar
                nvar
        end
        pre = map(a -> get(var_map, a, a), replicate(variables(mon), exponents(mon)))
        preduce(fefe, pre; init=1)
end

function lifting_vars!(m, eqs, var_map)
        nnn = Dict()
        mo = [mon for q in eqs for mon in monomials(q)]

        for mon in mo
                v = lift_monon!(m, var_map, mon)
                nnn[mon] = v
        end

        nnn
end

function solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer(), init=θ)

        m = Model(optimizer)

        slim = sin_min_max.(θl, θh)
        clim = cos_min_max.(θl, θh)

        ids = eachindex(d)
        @variable(m, clim[i][1] <= c[i in ids] <= clim[i][2])
        @variable(m, slim[i][1] <= s[i in ids] <= slim[i][2])

        set_start_value.(c, cos.(init))
        set_start_value.(s, sin.(init))

        E, pc, ps = build_eqs(d, r, α, M)
        E = mapcoefficients.(c -> (abs(c) > (1e-12)) ? c : 0.0, E)

        var_map = Dict([pc; ps] .=> [c; s])
        ls = lifting_vars!(m, E, var_map)

        @constraint m lift_poly.(Ref(ls), E) .== 0

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






function matrix_kinematic_bounds(d, r, α, θl, θh, M;
        optimizer=_matrix_default_optimizer(), init=θ)

        m = Model(optimizer)

        slim = sin_min_max.(θl, θh)
        clim = cos_min_max.(θl, θh)

        ids = eachindex(d)
        @variable(m, clim[j][1] <= c[j in ids] <= clim[j][2])
        @variable(m, slim[j][1] <= s[j in ids] <= slim[j][2])

        set_start_value.(c, cos.(init))
        set_start_value.(s, sin.(init))

        LHS, RHS = build_mat_eqs(d, r, α, c, s)
        @constraint m LHS .== M * RHS

        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

	g = nothing
	angmax = zeros(length(d))
	for i in eachindex(d)
		g = @constraint(m, s[i] >= 0)
		@objective m Max sum(2 * (1 - c[i] * cos(0) - s[i] * sin(0)))
		optimize!(m)
		angmax[i] = atan(value(s[i]), value(c[i]))
		JuMP.delete(m, g)
	end
angmin = zeros(length(d))
	for i in eachindex(d)
		g = @constraint(m, s[i] <= 0)
		@objective m Max sum(2 * (1 - c[i] * cos(0) - s[i] * sin(0)))
		optimize!(m)
		angmin[i] = atan(value(s[i]), value(c[i]))
		JuMP.delete(m, g)
	end
angmin, angmax
end











function find_bounds(d, r, α, θl, θh; optimizer=_default_optimizer())
        m = Model(optimizer)

        slim = sin_min_max.(θl, θh)
        clim = cos_min_max.(θl, θh)

        ids = eachindex(d)
        @variable(m, clim[i][1] <= c[i in ids] <= clim[i][2], start = 0)
        @variable(m, slim[i][1] <= s[i in ids] <= slim[i][2], start = 0)
        @variable(m, M[1:4, 1:4], start = 0)

        @polyvar pc[ids] ps[ids]
        E, pc, ps = build_eqs_simple(d, r, α, pc, ps)

        E = mapcoefficients.(c -> (abs(c) > (1e-12)) ? c : 0.0, E)
        var_map = Dict([pc; ps] .=> [c; s])
        ls = lifting_vars!(m, E, var_map)

        M = lift_poly.(Ref(ls), E)

        @constraint m c .^ 2 .+ s .^ 2 .== 1
        @constraint m (c .+ 1) .* tan.(θl / 2) .- s .<= 0
        @constraint m (c .+ 1) .* tan.(θh / 2) .- s .>= 0

        Lmin = fill(-1e6, 4, 4)
        Lmax = fill(1e6, 4, 4)
        for i in 1:3, j in 4:4
                @objective m Min M[i, j]
                optimize!(m)

                Lmin[i, j] = value(M[i, j])

                @objective m Max M[i, j]
                optimize!(m)

                Lmax[i, j] = value(M[i, j])
        end


        Lmin, Lmax
end

