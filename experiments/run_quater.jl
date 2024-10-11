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

<<<<<<< HEAD
using HomotopyContinuation, LinearAlgebra

function fpose(d, r, α, θl, θh)
    x = θl .+ rand(length(θl)) .* (θh .- θl)
    Q = prod(dh_quaternion.(x, d, α, r))

    x, Q, [vec(Q.r); _quaternion_to_translation(Q)] 
end


function fapose(d, r, α, θl, θh)
=======
using DynamicPolynomials

function con_mat(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_matrix_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_matrix_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

    res = T(3)*T(4)*T(5)-iT(2)*iT(1)*M*iT(7)*iT(6)
    #res = prod(T, eachindex(d))-M

    res
end

function rat_feas_poses(d, r, α, θl, θh)
>>>>>>> 2be260c6dcd975c62e05a764f862797e9a437f0e
    x = θl .+ rand(length(θl)) .* (θh .- θl)
    Q = prod(dh_quaternion.(x, d, α, r))
    M = prod(dh_matrix.(x, d, α, r))

<<<<<<< HEAD
    M,Q
end

=======
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


>>>>>>> 2be260c6dcd975c62e05a764f862797e9a437f0e
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
<<<<<<< HEAD
M, Q = feas_pose(d, r, α, θl, θh)

function pose_constraint(M, d, r, α, c, s)
    fwd, rev = build_pose_constraint_quaternion(d, r, α, c, s)
    Mm = DualQuaternion(Quaternion(M[1], M[2], M[3], M[4]), Quaternion(M[5], M[6], M[7], M[8]))
    Q = prod(fwd) - Mm * prod(rev)

    #T(i) = dh_quaternion(c[i], s[i], d[i], α[i], r[i])
    #Q = prod(T, eachindex(d))     
    vec(Q) 
end


@var c[1:7] s[1:7] M[1:8]
d, r, α, θl, θh, w, θ = simple_mainp(7)

y, My, my = fpose(d, r, α, θl, θh)
cy = cos.(y / 2)
sy = sin.(y / 2)

x, Mx, mx = fpose(d, r, α, θl, θh)
cx = cos.(x / 2)
sx = sin.(x / 2)

constrs = pose_constraint(M, d, r, α, c, s)
constrs = [vec(constrs); c .^ 2 .+ s .^ 2 .- 1]
sys = System(constrs; variables=[c; s], parameters=vec(M))

@show norm([constr(c=>cx, s=>sx, vec(M)=>vec(Mx)) for constr in constrs])
@show norm([constr(c=>cy, s=>sy, vec(M)=>vec(My)) for constr in constrs])

#@show rcold = solve(sys; target_parameters=vec(Mx))
#@show rself = solve(sys, solutions(rcold); start_parameters=vec(Mx), target_parameters=vec(My))
@show rwarm = solve(sys, [[cx; sx]]; start_parameters=vec(Mx), target_parameters=vec(My))

nothing
=======
M, Q = rat_feas_poses(d, r, α, θl, θh)
>>>>>>> 2be260c6dcd975c62e05a764f862797e9a437f0e
