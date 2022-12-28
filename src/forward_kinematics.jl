include("denavit_hartenberg.jl")

function solve_forward_kinematics(x, d, r, α)
    prod(dh_t.(x, d, α, r))
end

function random_feasible_pose(d, r, α, θh, θl)
    x = θl .+ rand(length(θl)) .* (θh .- θl)
    prod(dh_t.(x, d, α, r))
end
