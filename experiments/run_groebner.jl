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
include("../src/trutman.jl")

using Groebner
using DynamicPolynomials

d, r, α, θl, θh, w, θ = params_random_orth(6)
Mflt, _ = random_feasible_pose_hq(d, r, α, θl, θh)
M = rationalize_transformation(Mflt)

@polyvar c[1:7] s[1:7]
F, R = build_pose_constraint(d, r, α, c, s)
constrs = view(prod(F) - M * prod(R), 1:3, :)

