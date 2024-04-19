include("./benchmark_ik_method.jl")

include("./rand_parameters.jl")
include("./kuka_parameters.jl")
include("./icub_arm_parameters.jl")


include("../src/utils.jl")
include("../src/modelling.jl")
include("../src/quaternion.jl")
include("../src/trutman.jl")
include("../src/denavit_hartenberg.jl")
include("../src/forward_kinematics.jl")
include("../src/lift.jl")
include("../src/inverse_kinematics.jl")
include("../src/local_kinematics.jl")

using Groebner
using DynamicPolynomials

function rat_feas_pose(d, r, α, θl, θh; tol=1e-2)
    x = θl .+ rand(length(θl)) .* (θh .- θl)

    prod(dh_matrix_rat.(x, d, α, r; tol))
end

function rat_pose_constraint_half(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_matrix_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_matrix_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

    fwd, rev = _split_manipulator(eachindex(d))

    res = prod(T, fwd) - M*prod(iT, rev)
    view(res, 1:3, :)[:]
end

function simple_mainp(joint_count)
    r = rand(1:99, joint_count)
    d = rand(1:99, joint_count)
    α = rand([-π / 2, 0, π / 2], joint_count)

    θl = vec(fill(-3, joint_count))
    θh = vec(fill(+3, joint_count))

	w = normalize(ones(joint_count), 1)
	θ = zeros(joint_count)

    d, r, α, θl, θh, w, θ
end

d, r, α, θl, θh, w, θ = params_kuka_iiwa() # params_icub_v2(8)
M = rationalize_transformation(random_feasible_pose(d, r, α, θl, θh); tol=1e-3)

@polyvar c[eachindex(d)] s[eachindex(d)]
constrs = rat_pose_constraint_half(M, d, r, α, c, s; tol=1e-2)
@show constrs = [vec(constrs); c .^ 2 .+ s .^ 2 .- 1]

gbA = groebner(constrs)
gb12 = filter(x->maxdegree(x) <= 2, gbA)

gb_inverse_kinematics(gb12, c, s, θl, θh, θ, w)
