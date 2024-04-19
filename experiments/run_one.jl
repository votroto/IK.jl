include("./benchmark_ik_method.jl")

include("./rand_parameters.jl")
include("./kuka_parameters.jl")
include("./icub_arm_parameters.jl")


include("../src/utils.jl")
include("../src/modelling.jl")
include("../src/quaternion.jl")
include("../src/trutman.jl")
include("../src/denavit_hartenberg.jl")
include("../src/forward_kinematics.jl")
include("../src/lift.jl")
include("../src/inverse_kinematics.jl")
include("../src/local_kinematics.jl")

d, r, α, θl, θh, w, θ = params_icub_v2(8)

function _random_feasible_pose_hq(d, r, α, θl, θh)
    x = θl .+ rand(length(θl)) .* (θh .- θl)

    prod(dh_matrix.(x, d, α, r)), prod(dh_quaternion.(x, d, α, r))
end
using Random



H, Q = _random_feasible_pose_hq(d, r, α, θl, θh)

sol, obj = local_inverse_kinematics(d, r, α, θl, θh, H, θ, w)
#xh, oh, sh, th = solve_inverse_kinematics(d, r, α, θl, θh, H, θ, w)
xh, oh, sh, th = solve_inverse_kinematics(d, r, α, θl, θh, H, θ, w; init=sol)

#=
xh, oh, sh, th = solve_inverse_kinematics(d, r, α, θl, θh, H, θ, w; init=sol)
xq, oq, sq, tq = solve_inverse_kinematics(d, r, α, θl/2, θh/2, Q, θ/2, w; init=xh)
=#
