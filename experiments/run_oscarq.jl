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

using Oscar
using LinearAlgebra
using DynamicPolynomials

function rat_feas_pose_quaternion(d, r, α, θl, θh; tol=1e-2)
    x = θl .+ rand(length(θl)) .* (θh .- θl)

    prod(dh_quaternion_rat.(x, d, α, r; tol))
end

function rat_pose_constraint_half_quaternion(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_quaternion_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_quaternion_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

    fwd, rev = _split_manipulator(eachindex(d))

    res = prod(T, fwd) - prod(iT, rev)*M
    vec(res)
end

function simple_mainp(joint_count)
    r = rand(1:99, joint_count)
    d = rand(1:99, joint_count)
    α = rand([-π / 2, π / 2], joint_count)

    θl = vec(fill(-3, joint_count))
    θh = vec(fill(+3, joint_count))

	w = normalize(ones(joint_count), 1)
	θ = zeros(joint_count)

    d, r, α, θl, θh, w, θ
end

d, r, α, θl, θh, w, θ = params_kuka_iiwa() #simple_mainp(7)
M = rat_feas_pose_quaternion(d, r, α, θl, θh)

QQ, (C1, C2, C3, C4, C5, C6, C7, S1, S2, S3, S4, S5, S6, S7) = Oscar.polynomial_ring(Oscar.QQ, ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "S1", "S2", "S3", "S4", "S5", "S6", "S7"])
c, s = [C1, C2, C3, C4, C5, C6, C7], [S1, S2, S3, S4, S5, S6, S7]

@polyvar zc[1:7] zs[1:7]

zconstrs = rat_pose_constraint_half_quaternion(M, d, r, α, zc, zs; tol=1e-2)

constrs = [z([zc;zs]=>[c;s]) for z in zconstrs]
@show constrs = [vec(constrs); c .^ 2 .+ s .^ 2 .- 1]

gbA = groebner_basis(ideal(constrs); complete_reduction=true)

using Serialization

serialize("oscarqgba.bin", gbA)
