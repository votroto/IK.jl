include("kuka_parameters.jl")
include("../src/forward_kinematics.jl")
include("../src/inverse_kinematics.jl")
include("../src/local_kinematics.jl")

rand_pose = random_feasible_pose(d, r, α, θh, θl)


loc, locobj = local_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w)
glo, gloobj, = solve_inverse_kinematics(d, r, α, θl, θh, rand_pose, θ, w; init=loc)

println(loc)
println(locobj)
println(glo)
println(gloobj)