
using HomotopyContinuation, LinearAlgebra
#=
include("../src/utils.jl")
include("../src/quaternion.jl")
include("../src/modelling.jl")
include("../src/trutman.jl")
include("../src/denavit_hartenberg.jl")
include("../src/forward_kinematics.jl")
include("../src/lift.jl")
include("../src/inverse_kinematics.jl")
include("../src/local_kinematics.jl")
=#
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

function pose_constraint_half(M, d, r, α, c, s; tol=1e-4)
    T(i) = dh_matrix(c[i], s[i], d[i], α[i], r[i])
    iT(i) = dh_matrix_inverse(c[i], s[i], d[i], α[i], r[i])

    fwd, rev = _split_manipulator(eachindex(d))

    res = prod(T, fwd) - M * prod(iT, rev)
    view(res, 1:3, :)[:]
end


@var c[1:7] s[1:7] p[1:3] Mv[1:4,1:4]

d, r, α, θl, θh, w, θ = simple_mainp(7) # params_icub_v2(8)

x, M = fpose(d, r, α, θl, θh)
x2, M2 = fpose(d, r, α, θl, θh)


cv = cos.(x)
sv = sin.(x)

constrs = pose_constraint_half(Mv, d, r, α, c, s)
constrs = [vec(constrs); c .^ 2 .+ s .^ 2 .- 1]

cF = System(constrs; variables=[c; s], parameters=vec(Mv))

@show rnext = solve(cF, [[cv; sv]]; start_parameters=vec(M), target_parameters=vec(M))

subs.(constrs, c=>cv, s=>sv, vec(Mv)=>vec(M))




