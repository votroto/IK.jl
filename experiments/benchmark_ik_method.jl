include("../src/forward_kinematics.jl")

no_warm_start(d, r, α, θl, θh, desired, θ, w) = θ, 0.0

function stats_sample(method, d, r, α, θl, θh, w, θ; warm_start)
    desired = random_feasible_pose(d, r, α, θl, θh)
    local_x, local_obj = warm_start(d, r, α, θl, θh, desired, θ, w)
    
    x, obj, ret, tim = method(d, r, α, θl, θh, desired, θ, w; init=local_x)
    actual = solve_forward_kinematics(x, d, r, α)
    loc_err, rot_err = pose_error(desired, actual)

    loc_err, rot_err, obj, local_obj, tim, ret
end