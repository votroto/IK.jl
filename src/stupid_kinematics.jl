using Base.Iterators: peel, drop, take
using DynamicPolynomials
using Gurobi
using SCIP
using JuMP
using Ipopt

include("denavit_hartenberg.jl")

function _default_optimizer()
        optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2)
end


function _split_manipulator(ids)
        mid = div(length(ids), 2)
        f, s = take(ids, mid), drop(ids, mid)
        f, reverse(collect(s))
end

function build_eqs(d, r, α, M,m)
        T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
        iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

        ids = eachindex(d)
        fwd, rev = _split_manipulator(ids)

        @polyvar c[ids] s[ids]

        e = prod(T, fwd) .- M * prod(iT, rev)
	e = mapcoefficients.(c -> (abs(c) > (1e-12)) ? c : 0.0, e)

	@show maxdeg = maximum(maxdegree.(e))
	@show monovec = monomials([c;s],0:div(maxdeg,2,RoundUp))

	@show length(monovec)
println("aa")
	@variable m mv[1:length(monovec)]
println("aa")
	@variable m my[1:length(monovec),1:length(monovec)]
println("aa")
	@constraint m my .== mv * mv'
println("aa")
	monomat = monovec * monovec'
println("aa")
	mdict = Dict(monomat[:] .=> my[:])
println("aa")
	[dot(coefficients(ee), get.(Ref(mdict), monomials(ee), missing)) for ee in e], get.(Ref(mdict), c, missing),get.(Ref(mdict), s, missing)

end

function stupid_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
        optimizer=_default_optimizer())

        m = Model(optimizer)

        e, c, s = build_eqs(d, r, α, M, m)

        @constraint m e .== 0

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
