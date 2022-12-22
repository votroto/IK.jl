params = "kuka_parameters"

include("$(params).jl")
include("../src/forward_kinematics.jl")
include("../src/inverse_kinematics.jl")

open("$params.scip3.txt","a") do f
    for i in 1:100
        rand_pose = random_feasible_pose(d, r, α, θh, θl)
        sol, obj, stat, tim = solve_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w)
	act_pose = solve_forward_kinematics(sol, d, r, α)
	er = sum(abs.(rand_pose .- act_pose))/16
        println("$tim\t$er")
        println(f, "$tim\t$er")
    end
end

#=
open("$params.2.infeas.txt","a") do f
    for i in 1:100
	pis = ones(length(d)) * π
        rand_pose = random_feasible_pose(d, r, α, -pis, pis)
        sol, obj, stat, tim = solve_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w)
	act_pose = solve_forward_kinematics(sol, d, r, α)
	er = sum(abs.(rand_pose .- act_pose))/16
        println("$tim\t$er")
        println(f, "$tim\t$er")
    end
end
=#