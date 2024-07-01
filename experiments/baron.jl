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

function rat_feas_pose(d, r, α, θl, θh; tol=1e-2)
    x = θl .+ rand(length(θl)) .* (θh .- θl)

    prod(dh_matrix_rat.(x, d, α, r; tol))
end

function rat_pose_constraint_half(M, d, r, α, c, s; tol=1e-2)
    T(i) = dh_matrix_rat(c[i], s[i], d[i], α[i], r[i]; tol)
    iT(i) = dh_matrix_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol)

    fwd, rev = _split_manipulator(eachindex(d))

    res = prod(T, fwd) - M*prod(iT, rev)
    view(res, 1:3, :)[:]
end

function simple_mainp(joint_count)
    r = rand(1:99, joint_count)
    d = rand(1:99, joint_count)
    α = rand_between(-π / 2, π / 2, joint_count)

    θl = vec(fill(-3, joint_count))
    θh = vec(fill(+3, joint_count))

	w = normalize(ones(joint_count), 1)
	θ = zeros(joint_count)

    d, r, α, θl, θh, w, θ
end




function zzz(i, f)

    d, r, α, θl, θh, w, θ = simple_mainp(7) # params_icub_v2(8)
    M = rationalize_transformation(random_feasible_pose(d, r, α, θl, θh); tol=1e-3)
    
    @polyvar c[eachindex(d)] s[eachindex(d)]
    constrs = rat_pose_constraint_half(M, d, r, α, c, s; tol=1e-2)
    constrs = [vec(constrs); c .^ 2 .+ s .^ 2 .- 1]
    
    
    
    ids = 1:7
    
    println(f, "set IDS;")
    println(f, "var c{IDS} <= 1, >= -1;")
    println(f, "var s{IDS} <= 1, >= -1;")

    for i in eachindex(constrs)
        println(f, "s.t. eq$i: $(constrs[i]) = 0;")
    end

    for i in ids
        println(f, "s.t. ub$i: $(lin_angdiff_proxy(c[i], s[i], θl[i])) >= 0;")
        println(f, "s.t. lb$i: $(lin_angdiff_proxy(c[i], s[i], θh[i])) <= 0;")
    end
    
    println(f, "minimize move: $(sum(lin_abs_angdiff_proxy.(c, s, θ, w)));")

    println(f, "data;")
    println(f, "set IDS := 1 2 3 4 5 6 7;")
    
end


for i in 1:100
   

    open("baron/mod_$i.txt", "w") do f
        zzz(i, f)
    end
end