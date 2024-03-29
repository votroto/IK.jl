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

function rat_pose_constraint_half(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_matrix_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_matrix_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

    fwd, rev = _split_manipulator(eachindex(d))

    res = prod(T, fwd) - prod(iT, rev)*M
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

d, r, α, θl, θh, w, θ = simple_mainp(7)
M = rationalize_transformation(random_feasible_pose(d, r, α, θl, θh))

QQ, (C1, C2, C3, C4, C5, C6, C7, S1, S2, S3, S4, S5, S6, S7) = Oscar.polynomial_ring(Oscar.QQ, ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "S1", "S2", "S3", "S4", "S5", "S6", "S7"])
c, s = [C1, C2, C3, C4, C5, C6, C7], [S1, S2, S3, S4, S5, S6, S7]

constrs = rat_pose_constraint_half(M, d, r, α, c, s; tol=1e-2)
@show constrs = [vec(constrs); c .^ 2 .+ s .^ 2 .- 1]

gbA = groebner_basis(ideal(constrs); complete_reduction=true)
gb12 = filter(x->total_degree(x) <= 2, gbA.gens |> collect)

gbB = groebner_basis(ideal(gb12); complete_reduction=true)
@show all(iszero, gbA.gens .- gbB.gens)
