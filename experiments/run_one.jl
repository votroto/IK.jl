include("./benchmark_ik_method.jl")

include("./rand_parameters.jl")
include("./kuka_parameters.jl")
include("./icub_arm_parameters.jl")


include("../src/utils.jl")
include("../src/modelling.jl")
include("../src/quaternion.jl")
include("../src/denavit_hartenberg.jl")
include("../src/forward_kinematics.jl")
include("../src/lift/lift_adhoc.jl")
include("../src/lift/lift_matrix.jl")
include("../src/inverse_kinematics.jl")
include("../src/local_kinematics.jl")

include("../src/quaternion.jl")

d, r, α, θl, θh, w, θ = params_kuka_iiwa()

desiredh, desiredq = random_feasible_pose_hq(d, r, α, θl, θh)

θi, obji = local_inverse_kinematics(d, r, α, θl, θh, desiredh, θ, w)

#xh, objh, reth, timh = solve_inverse_kinematics(d, r, α, θl, θh, desiredh, θi, w; lift_method=lift_matrix)
xq, objq, retq, timq = solve_inverse_kinematics(d, r, α, θl/2, θh/2, desiredq, θi, w; lift_method=lift_matrix_q)
