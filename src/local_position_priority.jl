include("./denavit_hartenberg.jl")

using Base.Iterators: peel, drop, take
using ADNLPModels
using Ipopt
using NLPModelsIpopt

function build_K(d, r, α, x)
    T(i) = dh_t(x[i], d[i], α[i], r[i])
    prod(map(T, eachindex(x)))
end

function build_obj(d, r, α, M, x, w, θ)
    K = build_K(d, r, α, x)
    Ka = K[1:3, 1:3]
    Ma = M[1:3, 1:3]

    return norm(Ma - Ka)^2 + (θ - x)' * diagm(w) * (θ - x)
end

function build_constr(d, r, α, M, x)
    K = build_K(d, r, α, x)
    Kx = K[1:3, 4]
    Mx = M[1:3, 4]

    return [norm(Mx - Kx)^2; x]
end

function local_position_priority(d, r, α, θl, θh, M, θ, w; init=θ)
    obj(x) = build_obj(d, r, α, M, x, w, θ)
    con(x) = build_constr(d, r, α, M, x)

    nlp = ADNLPModel(obj, init, con, [0; θl], [0; θh])
    stats = ipopt(nlp; tol=1e-3, print_level=0, max_iter=200)

    stats.solution, stats.objective
end