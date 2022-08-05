include("../src/inverse_kinematics.jl")
include("kuka_parameters.jl")

solve_inverse_kinematics(d, r, α, θl, θh, M_actual, θ, w)