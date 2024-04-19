using Base.Iterators: peel, drop, take
using ADNLPModels
using Ipopt
using NLPModelsIpopt

function _build_objective(x, w, θ)
    return sum(lin_abs_angdiff_proxy.(cos.(x), sin.(x), θ, w))
end

function _build_constraint(x, d, r, α, M)
    Fs, Rs = build_pose_constraint(d, r, α, cos.(x), sin.(x))
    E = vec(prod(Fs) - M * prod(Rs))

    return [x; E]
end

"""
    local_inverse_kinematics(d, r, α, θl, θh, M, θ, w; init=θ)

Computes a local inverse kinematics solution, starting from `init`.
"""
function local_inverse_kinematics(d, r, α, θl, θh, M, θ, w; init=θ)
    obj(x) = _build_objective(x, w, θ)
    con(x) = _build_constraint(x, d, r, α, M)
    z = vec(zero(M))

    nlp = ADNLPModel(obj, init, con, Float64[θl; z], Float64[θh; z])
    stats = ipopt(nlp; tol=1e-3, print_level=0, max_iter=200)

    stats.solution, stats.objective
end


function _build_constraintgb(x, cs, fs)
    E = [f(cs => [cos.(x); sin.(x)]) for f in fs]

    return [x; E]
end

"""
    local_inverse_kinematics(d, r, α, θl, θh, M, θ, w; init=θ)

Computes a local inverse kinematics solution, starting from `init`.
"""
function local_inverse_kinematicsgb(cs, fs, θl, θh, θ, w; init=θ)
    obj(x) = _build_objective(x, w, θ)
    con(x) = _build_constraintgb(x, cs, fs)
    z = fill(0., length(fs))

    nlp = ADNLPModel(obj, init, con, Float64[θl; -z], Float64[θh; z])
    stats = ipopt(nlp; tol=1e-3, max_iter=200)

    stats.solution, stats.objective
end