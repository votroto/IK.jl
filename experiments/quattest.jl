include("../src/quaternion.jl")
include("../src/denavit_hartenberg.jl")
include("../src/trutman.jl")
include("../src/forward_kinematics.jl")
include("kuka_parameters.jl")

n = 3
H, Q = random_feasible_pose_hq(rand(n), rand(n), randn(n), -3*rand(n), 3*rand(n))

#=
function rtq(r::Matrix)
    q = [tr(r) + 1, r[3, 2] - r[2, 3], r[1, 3] - r[3, 1], r[2, 1] - r[1, 2]] / (2 * sqrt(tr(r) + 1))
    #normalize!(q)
    Quaternion(q[1], q[2:4]...)
end
=#

rh = H[1:3,1:3]