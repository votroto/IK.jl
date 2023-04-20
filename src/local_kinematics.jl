using Base.Iterators: peel, drop, take
using ADNLPModels
using Ipopt
using NLPModelsIpopt

include("denavit_hartenberg.jl")

function build_obj(x)
        return sum(2 .* w .* (1 .- cos.(x) .* cos.(θ) .- sin.(x) .* sin.(θ)))
end

function build_constr(x, d, r, α, M, th)
        T(i) = dh_t(x[i], d[i], α[i], r[i])

        ids = eachindex(d)
        pp = prod(T, ids)

        return [x; norm(pp - M)]
end

function infeasibility(x, d, r, α, M)
        T(i) = dh_t(x[i], d[i], α[i], r[i])

        ids = eachindex(d)
        pp = prod(T, ids)

        return maximum(abs.(pp .- M))
end

function local_inverse_kinematics(d, r, α, θl, θh, M, θ, w; ang=θ)


        obj(x) = build_obj(x)
	con(x) = build_constr(x, d, r, α, M, θ)
        nlp = ADNLPModel(obj, ang,
                con,
                [θl; 0], [θh;0])

        stats = ipopt(nlp; tol=1e-3, print_level=0, max_iter=200)

	infeas = infeasibility(stats.solution, d, r, α, M)
        println("inf", infeasibility(ang, d, r, α, M), " ->", infeasibility(stats.solution, d, r, α, M))
	println("obj",obj(ang), "->", stats.objective)
	println()
        stats.solution, stats.objective, infeas
end