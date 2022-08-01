include("denavit_hartenberg.jl")

function solve_forward_kinematics(x, d, r, α)
    prod(dh_t.(x, d, α, r))
end