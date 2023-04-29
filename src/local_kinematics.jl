using Base.Iterators: peel, drop, take
using ADNLPModels
using Ipopt
using NLPModelsIpopt

function build_obj(x, w, θ)
    return sum(w .* lin_abs_angdiff_proxy.(cos.(x), sin.(x), θ))
end

function build_constr(x, d, r, α, M)
    Fs, Rs = build_eqs(d, r, α, cos.(x), sin.(x))
    E = prod(Fs) - M * prod(Rs)

    return [x; norm(E)]
end

function local_inverse_kinematics(d, r, α, θl, θh, M, θ, w; ang=θ)
    obj(x) = build_obj(x, w, θ)
    con(x) = build_constr(x, d, r, α, M)

    nlp = ADNLPModel(obj, ang, con, [θl; 0], [θh; 0])
    stats = ipopt(nlp; tol=1e-3, print_level=0, max_iter=200)

    stats.solution, stats.objective
end