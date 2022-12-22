include("kuka_parameters.jl")
include("../src/forward_kinematics.jl")
include("../src/ahomo_kinematics.jl")

rand_pose = random_feasible_pose(d, r, α, θh, θl)
homo_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w)