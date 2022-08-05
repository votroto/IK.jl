include("eight_parameters.jl")
include("../src/forward_kinematics.jl")
include("../src/inverse_kinematics.jl")

solve_inverse_kinematics(d, r, α, random_feasible_pose(d, r, α, θh, θl), θ, w)