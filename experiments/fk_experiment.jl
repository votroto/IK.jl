include("../src/forward_kinematics.jl")
include("kuka_parameters.jl")

sol = [
    -0.408833022287725
    -0.973777295554142
    -0.823238993601953
    -1.44436562369012
    -0.0076269233501136
    -1.246589257526
    1.16830383119891
]

solve_forward_kinematics(sol, d, r, α)