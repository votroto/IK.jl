using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using Ipopt
using MultivariatePolynomials

include("denavit_hartenberg.jl")

function _default_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2, "Presolve" => 2)
end


function _weak_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2, "OptimalityTol" => 1e-2, "MIPGap" => 1e-1)
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
        reduce(f, [reduce(f, h; init), reduce(f, t; init)]; init)
end

function lift_mono!(m, var_map, mon, mmm)
        ps = powers(mon)

        temp = 1
        nvar = 1
        for (v, p) in ps
                for j in 1:p
                        eee = temp * var_map[v]
                        nvar = get!(mmm, eee, @variable(m))#, lower_bound=-6, upper_bound=6))
                        @constraint m eee == nvar
                        temp = nvar
                end
        end
        nvar
end


function lift_monon!(m, var_map, mon, mmm)
        function fefe(a, b)
                a, b
		na = get(mmm, a, get(var_map, a, a))
                nb =  get(mmm, b, get(var_map, b, b))
		eee = na*nb
		if eee in keys(mmm)
			nvar = mmm[eee]
else
nvar = mmm[eee]= @variable(m,lower_bound=-1,upper_bound=1,base_name=string(eee))

                @constraint m eee == nvar
end
                nvar
        end

        pre = replicate(variables(mon), exponents(mon))

        preduce(fefe, pre; init=1)
end


function lift_poly(lift_vars, poly)
        es = monomials(poly)
        cs = coefficients(poly)
        ns = map(e -> lift_vars[e], es)

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

function lifting_vars!(m, eqs, var_map)
        mmm = Dict()
        nnn = Dict()
        mo = reverse(sortmonovec(unique([mon for q in eqs for mon in monomials(q)][:]))[2])

        for mon in mo
                v = lift_monon!(m, var_map, mon, mmm)
                nnn[mon] = v
        end
        nnn
end

function solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
        optimizer=_scip_optimizer(), init=zeros(length(d)))

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
