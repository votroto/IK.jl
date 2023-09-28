using Revise

using IK

include("../experiments/icub_arm_parameters.jl")
include("../src/local_position_priority.jl")

d, r, α, θl, θh, w, θ = params_icub_v2(7)
pose_gen = IK.random_feasible_pose
warm_start = local_inverse_kinematics
desired = pose_gen(d, r, α, θl, θh)
local_x, local_obj = local_position_priority(d, r, α, θl, θh, desired, θ, w)
