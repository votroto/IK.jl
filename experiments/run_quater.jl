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

using DynamicPolynomials

function con_mat(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_matrix_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_matrix_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

    res = T(3)*T(4)*T(5)-iT(2)*iT(1)*M*iT(7)*iT(6)
    #res = prod(T, eachindex(d))-M

    res
end

function rat_feas_poses(d, r, α, θl, θh)
    x = θl .+ rand(length(θl)) .* (θh .- θl)

    M = prod(dh_matrix.(x, d, α, r))
    Q = prod(dh_quaternion.(x, d, α, r))

    M, Q
end

function con_quat(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_quaternion_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_quaternion_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

    res = T(3)*T(4)*T(5)-iT(2)*iT(1)*M*iT(7)*iT(6)
#    res = prod(T, eachindex(d))-M

    res

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


d, r, α, θl, θh, w, θ = simple_mainp(7)
M, Q = rat_feas_poses(d, r, α, θl, θh)
