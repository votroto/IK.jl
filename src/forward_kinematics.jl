using LinearAlgebra

function solve_forward_kinematics(x, d, r, α)
    prod(dh_matrix.(x, d, α, r))
end

function random_feasible_pose(d, r, α, θl, θh)
    x = θl .+ rand(length(θl)) .* (θh .- θl)
    prod(dh_matrix.(x, d, α, r))
end

function random_feasible_pose_hq(d, r, α, θl, θh)
    x = θl .+ rand(length(θl)) .* (θh .- θl)
    
    prod(dh_matrix_rat.(x, d, α, r)), prod(dh_quaternion.(x, d, α, r))
end

function pose_error(A, B)
    pos_e = norm(A[1:3, 4] - B[1:3, 4])
    rot_e = acos(clamp(0.5 * (tr(A[1:3, 1:3] \ B[1:3, 1:3]) - 1), -1.0, 1.0))

    pos_e, rot_e
end
