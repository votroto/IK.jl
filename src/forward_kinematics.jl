using LinearAlgebra

include("denavit_hartenberg.jl")

function solve_forward_kinematics(x, d, r, α)
    prod(dh_t.(x, d, α, r))
end

function random_feasible_pose(d, r, α, θh, θl)
    x = θl .+ rand(length(θl)) .* (θh .- θl)
    prod(dh_t.(x, d, α, r))
end


function pose_error(A, B)
	pos_e = norm(A[1:3, 4] - B[1:3, 4])
	rot_e = acos(0.5 * (tr(A[1:3,1:3]\B[1:3,1:3])-1))

	pos_e, rot_e
end
