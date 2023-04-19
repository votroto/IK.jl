params = "replication_parameters"

include("$(params).jl")
include("../src/forward_kinematics.jl")
include("../src/inverse_kinematics.jl")
include("../src/inverse_kinematics_matrix.jl")
include("../src/inverse_kinematics_poly.jl")
include("../src/local_kinematics.jl")

r, d, α, θl, θh = random_manipulator(5)
w = normalize(ones(5), 1)
θ = zeros(5)
M = random_feasible_pose(d, r, α, θl, θh)
rand_pose = random_feasible_pose(d, r, α, θl, θh)
init, lobj, = local_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w)
sol, obj, stat, tim = solve_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w; init)

p1 = solve_forward_kinematics(sol[1:3], d[1:3], r[1:3], α[1:3])
p2 = solve_forward_kinematics(sol[4:5], d[4:5], r[4:5], α[4:5])

p2 = rand_pose / p1
#=
@show r, d, α, θl, θh, w, θ = get_icub_v2(7)
open("$params.v2.R7.$(rand(UInt32)).infeas.txt","a") do f
#open("$params.$(rand(UInt32)).infeas.txt","a") do f
    for i in 1:300
        rand_pose = random_feasible_pose(d, r, α,  θl, θh)
	init,lobj, = local_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w)
        sol, obj, stat, tim = solve_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w; init)
	act_pose = solve_forward_kinematics(sol, d, r, α)
	println("$tim\t$obj")
	println(f, "$tim\t$obj")
	#@show pe, re = pose_error(rand_pose, act_pose)
        #println("$tim\t$pe\t$re\t$obj\t$lobj")
	#println(f, "$tim\t$pe\t$re\t$obj\t$lobj")
    end
end
=#