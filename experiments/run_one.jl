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
include("../src/lift/lift_tree.jl")

include("../src/quaternion.jl")

d, r, α, θl, θh, w, θ = params_kuka_iiwa()

desiredh, desiredq = random_feasible_pose_hq(d, r, α, θl, θh)


xh, objh, reth, timh = solve_inverse_kinematics(d, r, α, θl, θh, desiredh, θ, w; lift_method=lift_tree)

xq, objq, retq, timq = solve_inverse_kinematics(d, r, α, θl, θh, desiredq, θ, w; lift_method=lift_tree_q)