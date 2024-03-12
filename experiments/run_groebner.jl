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

using Groebner
using DynamicPolynomials

d, r, α, θl, θh, w, θ = params_kuka_iiwa()
Mflt, _ = random_feasible_pose_hq(d, r, α, θl, θh)
M = rationalize_transformation(Mflt)



function trutman_pose_constraint(M, d, r, α, c, s)
    T(i) = dh_matrix_rat(c[i], s[i], d[i], α[i], r[i])
    iT(i) = dh_matrix_rat_inverse(c[i], s[i], d[i], α[i], r[i])

    res = T(3)*T(4)*T(5)-iT(1)*iT(2)*M*iT(7)*iT(6)
    view(res, 1:3, :)[:]
end

@polyvar c[eachindex(d)] s[eachindex(d)]
constrs = trutman_pose_constraint(M, d, r, α, c, s)
constrs = [constrs; c .^ 2 .+ s .^ 2 .- 1]



