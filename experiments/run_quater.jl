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

function trutman_pose_constraint(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_matrix_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_matrix_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

    #res = T(3)*T(4)*T(5)-iT(2)*iT(1)*M*iT(7)*iT(6)
    res = prod(T, eachindex(d))-M
    view(res, 1:3, :)[:]
end

function rat_feas_pose_quaternion(d, r, α, θl, θh; tol=1e-2)
    x = θl .+ rand(length(θl)) .* (θh .- θl)

    prod(dh_quaternion_rat.(x, d, α, r; tol))
end

function trutman_pose_constraint_quaternion(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_quaternion_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_quaternion_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

    res = T(3)*T(4)*T(5)-iT(2)*iT(1)*M*iT(7)*iT(6)
#    res = prod(T, eachindex(d))-M

    vec(res)
    
end

#r = [0; 0; 0; 0; 0; 0; 0]
#d = [340; 0; 400; 0; 400; 0; 126]
#α = [-π / 2, π / 2, -π / 2, π / 2, -π / 2, π / 2, 0]
#θl = [-2.9671; -2.0944; -2.9671; -2.0944; -2.9671; -2.0944; -3.0543]
#θh = [2.9671; 2.0944; 2.9671; 2.0944; 2.9671; 2.0944; 3.0543]



#=
r = [47,22,82,19,85,15,89]
d = [74,26,11,86,81,16,35]
α = [-0.5073481922772854, 2.498091544796509, 0.8097835725701669, 2.651635327336065, 2.498091544796509, -2.882187579277969, 0.2213144423477913]

M = [
    41//105	 -88//105	8//21	-11807//100
    8//21	  11//21	16//21	389//100
    -88//105 -16//105	11//21	-24931//100
    0 0 0 1
]
=#
#
#M = Rational{BigInt}[19569//43481 25200//43481 29540//43481 -4792//17; 33600//43481 -27569//43481 1260//43481 14030//29; 19460//43481 22260//43481 -31881//43481 9110//33; 0 0 0 1]
#M = rat_feas_pose_quaternion(d, r, α, θl, θh; tol=1e-2)




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
M = rat_feas_pose_quaternion(d, r, α, θl, θh) # rationalize_transformation(random_feasible_pose(d, r, α, θl, θh))

#M = rat_feas_pose(d, r, α, θl, θh; tol=1e-2)

@polyvar c[eachindex(d)] s[eachindex(d)]
constrs = trutman_pose_constraint(M, d, r, α, c, s; tol=1e-2)
constrs = [constrs; c .^ 2 .+ s .^ 2 .- 1]

gbA = groebner(constrs)
gb12 = filter(x->maxdegree(x) <= 2, gbA)

lx, lo = local_inverse_kinematics(d, r, α, θl, θh, M, θ, w; init=θ)
gb_inverse_kinematics(gb12, c, s, θl, θh, θ, w; init=lx)

#then find the degree two polys
#compare to the original by recomputing GB and assert equal
#try with quats