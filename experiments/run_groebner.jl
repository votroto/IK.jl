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

    res = prod(T, fwd) - prod(iT, rev)*R
    view(res, 1:3, :)[:]
end

function rat_pose_constraint_full(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_matrix_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_matrix_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

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

d, r, α, θl, θh, w, θ = params_kuka_iiwa() # simple_mainp(8)
M = rationalize_transformation(random_feasible_pose(d, r, α, θl, θh))

#M = rat_feas_pose(d, r, α, θl, θh; tol=1e-2)

@polyvar c[eachindex(d)] s[eachindex(d)]
#constrs = trutman_pose_constraint(M, d, r, α, c, s; tol=1e-2)
constrs = rat_pose_constraint_full(M, d, r, α, c, s; tol=1e-2)
@show constrs = [vec(constrs); c .^ 2 .+ s .^ 2 .- 1]

gbA = groebner(constrs)
gb12 = filter(x->maxdegree(x) <= 2, gbA)

gbB = groebner(gb12)

allb = [gb12[: .!= i] for i in 1:45]

#=
filename_date = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

open("constrs$filename_date.txt", "w") do f
    println.(f,constrs)
end

open("gbA$filename_date.txt", "w") do f
    println.(f,gbA)
end

open("gb12$filename_date.txt", "w") do f
    println.(f,gb12)
end

open("gbB$filename_date.txt", "w") do f
    println.(f,gbB)
end
=#

#=
lx, lo = local_inverse_kinematics(d, r, α, θl, θh, M, θ, w; init=θ)
gb_inverse_kinematics(gb12, c, s, θl, θh, θ, w; init=lx)
=#

#then find the degree two polys
#compare to the original by recomputing GB and assert equal
#try with quats