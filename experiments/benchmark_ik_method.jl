include("../src/forward_kinematics.jl")

mean(xs) = sum(xs) / length(xs)
no_warm_start(d, r, α, θl, θh, desired, θ, w) = θ

function pose_error_sample(method, d, r, α, θl, θh, w, θ; warm_start)
    desired = random_feasible_pose(d, r, α, θl, θh)
    local_x = warm_start(d, r, α, θl, θh, desired, θ, w)
    
    x = method(d, r, α, θl, θh, desired, θ, w; init=local_x)
    M_actual = solve_forward_kinematics(x, d, r, α)
    
    pose_error(desired, M_actual)
end

#function average_pose_error(method, params, samples; warm_start=no_warm_start)
#    ss = [pose_error_sample(method, params; warm_start) for i in 1:samples]
#    mean(ss)
#end

