using Base.Iterators: peel, drop, take
using ADNLPModels
using Ipopt
using NLPModelsIpopt

include("denavit_hartenberg.jl")
include("modelling.jl")

function build_obj(x, w, θ)
        return sum(w .* lin_angdiff_approx.(cos.(x), sin.(x), θ))
end

function build_constr(x, d, r, α, M)
        T(i) = dh_t(x[i], d[i], α[i], r[i])

        ids = eachindex(d)
        pp = prod(T, ids)

        return [x; norm(pp - M)]
end

function local_inverse_kinematics(d, r, α, θl, θh, M, θ, w; ang=θ)
        obj(x) = build_obj(x, w, θ)
	con(x) = build_constr(x, d, r, α, M)

        nlp = ADNLPModel(obj, ang, con, [θl; 0], [θh;0])
        stats = ipopt(nlp; tol=1e-3, print_level=0, max_iter=200)

        stats.solution, stats.objective
end