include("./benchmark_ik_method.jl")

include("./rand_parameters.jl")

include("../src/utils.jl")
include("../src/modelling.jl")
include("../src/denavit_hartenberg.jl")
include("../src/forward_kinematics.jl")
include("../src/lift/lift_adhoc.jl")
include("../src/inverse_kinematics.jl")
include("../src/local_kinematics.jl")



d, r, α, θl, θh, w, θ = params_random_4rad(7)

desiredh = random_feasible_pose(d, r, α, θl, θh)

xh, objh, reth, timh = solve_inverse_kinematics(d, r, α, θl, θh, desiredh, θ, w)
