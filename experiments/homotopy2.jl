
using HomotopyContinuation, LinearAlgebra

include("../src/utils.jl")
include("../src/quaternion.jl")
include("../src/modelling.jl")
include("../src/trutman.jl")
include("../src/denavit_hartenberg.jl")
include("../src/forward_kinematics.jl")
include("../src/lift.jl")
include("../src/inverse_kinematics.jl")
include("../src/local_kinematics.jl")

function simple_mainp(joint_count)
    r = rand(1:99, joint_count)
    d = rand(1:99, joint_count)
    α = rand([π / 2, π / 2], joint_count)

    θl = vec(fill(-3, joint_count))
    θh = vec(fill(+3, joint_count))

    w = normalize(ones(joint_count), 1)
    θ = zeros(joint_count)

    d, r, α, θl, θh, w, θ
end

function fpose(d, r, α, θl, θh)
    x = θl .+ rand(length(θl)) .* (θh .- θl)
    x, prod(dh_matrix.(x, d, α, r))
end

function pose_constraint(M, d, r, α, c, s)
    fwd, rev = build_pose_constraint(d, r, α, c, s)
    res = prod(fwd) - M * prod(rev)
    view(res, 1:3, :)[:]
end


@var c[1:6] s[1:6] M[1:4,1:4]
d, r, α, θl, θh, w, θ = simple_mainp(6)

x, Mx = fpose(d, r, α, θl, θh)
cx = cos.(x)
sx = sin.(x)

y, My = fpose(d, r, α, θl, θh)
cy = cos.(y)
sy = sin.(y)


constrs = pose_constraint(M, d, r, α, c,s)
constrs = [vec(constrs); c .^ 2 .+ s .^ 2 .- 1]
sys = System(constrs; variables=[c; s], parameters=vec(M))

@show norm([constr(c=>cx, s=>sx, vec(M)=>vec(Mx)) for constr in constrs])
@show norm([constr(c=>cy, s=>sy, vec(M)=>vec(My)) for constr in constrs])

@show rcold = solve(sys; target_parameters=vec(Mx))
@show rself = solve(sys, solutions(rcold); start_parameters=vec(Mx), target_parameters=vec(My))
@show rwarm = solve(sys, [[cx; sx]]; start_parameters=vec(Mx), target_parameters=vec(My))

nothing