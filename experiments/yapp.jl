
function yapp_sample(method, d, r, α, θl, θh, w, θ; pose_gen, warm_start)
    desired = pose_gen(d, r, α, θl, θh)
    local_x, local_obj = warm_start(d, r, α, θl, θh, desired, θ, w)

    x, obj, ret, tim = method(d, r, α, θl, θh, desired, θ, w; init=local_x)
    actual = solve_forward_kinematics(local_x, d, r, α)
    loc_err, rot_err = pose_error(desired, actual)

    desired[1:3,4]..., loc_err, obj, ret, tim, local_obj
end