include("./benchmark_ik_method.jl")

include("./kuka_parameters.jl")
include("./rand_parameters.jl")
include("./icub_arm_parameters.jl")

include("../src/utils.jl")
include("../src/modelling.jl")
include("../src/jump_extensions.jl")
include("../src/denavit_hartenberg.jl")
include("../src/forward_kinematics.jl")
include("../src/inverse_kinematics.jl")
include("../src/local_kinematics.jl")
include("../src/lift/lift_matrix.jl")

using Random


pose_gen = random_feasible_pose
#params = params_random_4rad(7)
params = params_icub_v2(10)

warm_start = local_inverse_kinematics

d, r, α, θl, θh, w, θ = params
desired = pose_gen(d, r, α, θl, θh)
local_x, local_obj = warm_start(d, r, α, θl, θh, desired, θ, w)

#@show xs, objs, rets, tims = solve_inverse_kinematics(d, r, α, θl, θh, desired, θ, w; init=local_x)
@show x, obj, ret, tim = gb_inverse_kinematics(d, r, α, θl, θh, desired, θ, w; init=local_x)
@show x2, obj2, ret2, tim2 = gb_inverse_kinematics2(d, r, α, θl, θh, desired, θ, w; init=local_x)

#actual = solve_forward_kinematics(x, d, r, α)
#actuals = solve_forward_kinematics(xs, d, r, α)
#
#loc_err, rot_err = pose_error(desired, actual)
