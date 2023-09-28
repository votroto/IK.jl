using Revise

using IK

include("../experiments/icub_arm_parameters.jl")

function stats_sample(d, r, α, θl, θh, w, θ; pose_gen, warm_start)
    desired = pose_gen(d, r, α, θl, θh)
    local_x, local_obj = warm_start(d, r, α, θl, θh, desired, θ, w)

    x, obj, ret, tim = solve_inverse_kinematics(d, r, α, θl, θh, desired, θ, w; init=local_x)
    actual = solve_forward_kinematics(x, d, r, α)
    loc_err, rot_err = IK.pose_error(desired, actual)

    loc_err, rot_err, obj, local_obj, tim, ret
end

params = params_icub_v2(7)
pose_gen = IK.random_feasible_pose
warm_start = local_inverse_kinematics
result = stats_sample(params...; pose_gen, warm_start)