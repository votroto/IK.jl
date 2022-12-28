using Base.Iterators: peel, drop, take
using HomotopyContinuation
using ADNLPModels
using Ipopt
using NLPModelsIpopt

include("denavit_hartenberg.jl")

function build_obj(x, d, r, α, M, th)
        T(i) = dh_t(x[i], d[i], α[i], r[i])

        ids = eachindex(d)
        pp = prod(T, ids)

        return norm(pp - M)^2# + 1e-6 * sum((x.-th).^2)
end


function infeasibility(x, d, r, α, M)
        T(i) = dh_t(x[i], d[i], α[i], r[i])

        ids = eachindex(d)
        pp = prod(T, ids)

        return maximum(abs.(pp .- M))
end

function local_inverse_kinematics(d, r, α, θl, θh, M, θ, w)
        obj(x) = build_obj(x, d, r, α, M, θ)
        nlp = ADNLPModel(obj, θ,
                identity,
                θl, θh)

        stats = ipopt(nlp)

        infeas = infeasibility(stats.solution, d, r, α, M)
        print(infeas)
        stats.solution, stats.objective, infeas
end



function local_kinematic_bounds(d, r, α, θl, θh, M, θ)

        obj(x) = build_obj(x, d, r, α, M, θ)

        function inner(i, t)
                nlp = ADNLPModel(x -> t * x[i], θ,
                        x -> [x; obj(x)],
                        [θl; -1e-3], [θh; 1e-3])

                stats = ipopt(nlp)
                t*stats.objective
        end
	angmax = [inner(i, -1) for i in eachindex(d)]
	angmin = [inner(i, 1) for i in eachindex(d)]

	angmin, angmax
end