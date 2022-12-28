params = "kuka_parameters"

include("$(params).jl")
include("../src/forward_kinematics.jl")
include("../src/inverse_kinematics.jl")
include("../src/inverse_kinematics_matrix.jl")
include("../src/inverse_kinematics_poly.jl")
include("../src/local_kinematics.jl")

function bench(f)

        rand_pose = random_feasible_pose(d, r, α, θh, θl)

        ii, iiobj,infeas = local_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w)
	if infeas >= 1e-4
println(infeas)
		return
	end
        sol, obj, stat, tim = matrix_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w; init=ii)
	sol, obj, stat, tim = matrix_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w; )
        act_pose = solve_forward_kinematics(sol, d, r, α)
        er = sum(abs.(rand_pose .- act_pose)) / 16
        println("$tim\t$er")
        println(f, "$tim\t$er")
end


#lmin,lmax = find_bounds(d[1:2], r[1:2], α[1:2], θl[1:2], θh[1:2])

#bench(devnull)
rand_pose = random_feasible_pose(d, r, α, θh, θl)

#
sol, obj, = solve_inverse_kinematics_poly(d, r, α, θl, θh, rand_pose, θ, w)
sol, obj, = solve_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w)
#sol1, obj1, = solve_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w; )
#=
open("$params.test.txt", "a") do f
        for i in 1:1
                bench(f)
        end
end
=#
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