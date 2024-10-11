
using DynamicPolynomials, LinearAlgebra

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
    view(res, 1:3, :)
end


@polyvar c[1:6] s[1:6] M[1:4,1:4]
d, r, α, θl, θh, w, θ = simple_mainp(6)

x, Mx = fpose(d, r, α, θl, θh)
cx = cos.(x)
sx = sin.(x)

y, My = fpose(d, r, α, θl, θh)
cy = cos.(y)
sy = sin.(y)


constrsa = [DynamicPolynomials.map_coefficients(f-> round(f;digits=12), p) for p in pose_constraint(M, d, r, α, c,s)]
constrsb = c .^ 2 .+ s .^ 2 .- 1

open("/tmp/bertini.input", "w") do f
    println(f, "CONFIG")
    println(f, "USERHOMOTOPY: 1;")
    println(f, "END;")

    println(f, "INPUT")
    println(f, "variable c1, c2, c3, c4, c5, c6, s1, s2, s3, s4, s5, s6;")
    println(f, "function k1, k2, k3, k4, k5, k6, ik11, ik12, ik13, ik14, ik21, ik22, ik23, ik24, ik31, ik32, ik33, ik34;")
    println(f, "pathvariable t;")
    println(f, "parameter M11, M12, M13, M14, M21, M22, M23, M24, M31, M32, M33, M34;")

    
    println(f, "M11 = $(Mx[1,1])*(t) + $(Mx[1,1])*(t-1);")
    println(f, "M12 = $(Mx[1,2])*(t) + $(Mx[1,2])*(t-1);")
    println(f, "M13 = $(Mx[1,3])*(t) + $(Mx[1,3])*(t-1);")
    println(f, "M14 = $(Mx[1,4])*(t) + $(Mx[1,4])*(t-1);")
    println(f, "M21 = $(Mx[2,1])*(t) + $(Mx[2,1])*(t-1);")
    println(f, "M22 = $(Mx[2,2])*(t) + $(Mx[2,2])*(t-1);")
    println(f, "M23 = $(Mx[2,3])*(t) + $(Mx[2,3])*(t-1);")
    println(f, "M24 = $(Mx[2,4])*(t) + $(Mx[2,4])*(t-1);")
    println(f, "M31 = $(Mx[3,1])*(t) + $(Mx[3,1])*(t-1);")
    println(f, "M32 = $(Mx[3,2])*(t) + $(Mx[3,2])*(t-1);")
    println(f, "M33 = $(Mx[3,3])*(t) + $(Mx[3,3])*(t-1);")
    println(f, "M34 = $(Mx[3,4])*(t) + $(Mx[3,4])*(t-1);")


    println(f, "k1 = $(constrsb[1]);")
    println(f, "k2 = $(constrsb[2]);")
    println(f, "k3 = $(constrsb[3]);")
    println(f, "k4 = $(constrsb[4]);")
    println(f, "k5 = $(constrsb[5]);")
    println(f, "k6 = $(constrsb[6]);")
    
    println(f, "ik11 = $(constrsa[1, 1]);")
    println(f, "ik12 = $(constrsa[1, 2]);")
    println(f, "ik13 = $(constrsa[1, 3]);")
    println(f, "ik14 = $(constrsa[1, 4]);")
    println(f, "ik21 = $(constrsa[2, 1]);")
    println(f, "ik22 = $(constrsa[2, 2]);")
    println(f, "ik23 = $(constrsa[2, 3]);")
    println(f, "ik24 = $(constrsa[2, 4]);")
    println(f, "ik31 = $(constrsa[3, 1]);")
    println(f, "ik32 = $(constrsa[3, 2]);")
    println(f, "ik33 = $(constrsa[3, 3]);")
    println(f, "ik34 = $(constrsa[3, 4]);")

    println(f, "END;")
end

open("/tmp/bertini.start", "w") do f
    println(f, "1")
    println(f)
    println(f, cx[1], " 0.0;")
    println(f, cx[2], " 0.0;")
    println(f, cx[3], " 0.0;")
    println(f, cx[4], " 0.0;")
    println(f, cx[5], " 0.0;")
    println(f, cx[6], " 0.0;")
    println(f, sx[1], " 0.0;")
    println(f, sx[2], " 0.0;")
    println(f, sx[3], " 0.0;")
    println(f, sx[4], " 0.0;")
    println(f, sx[5], " 0.0;")
    println(f, sx[6], " 0.0;")
end