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

params = params_icub_v2(7)
params = params_kuka_iiwa()

warm_start = local_inverse_kinematics

d, r, α, θl, θh, w, θ = params
θ = θl .+ rand(length(d)) .* (θh .- θl)

desired = random_feasible_pose(d, r, α, θl, θh)

#@time lx, lo = zlocal_inverse_kinematics(d, r, α, θl, θh, desired, θ, w; init=θ)
#zx, zobj, zret, ztim = zsolve_inverse_kinematics(d, r, α, θl, θh, desired, θ, w)#; init=lx)
#ax, aobj, aret, atim = asolve_inverse_kinematics(d, r, α, θl, θh, desired, θ, w)#; init=lx)

zz, = solve_inverse_kinematics(d, r, α, θl, θh, desired, θ, w; optimizer=_feasible_optimizer())
zx, zobj, zret, ztim = solve_inverse_kinematics(d, r, α, θl, θh, desired, θ, w;optimizer=_scip())


function ferr(x,d,r,α)
    actual = solve_forward_kinematics(x, d, r, α)
    @show zloc_err, zrot_err = pose_error(desired, actual)
end

