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

function rat_feas_pose(d, r, α, θl, θh)
    x = θl .+ rand(length(θl)) .* (θh .- θl)

    prod(dh_matrix_rat.(x, d, α, r))
end

r = [0; 0; 0; 0; 0; 0; 0]
d = [340; 0; 400; 0; 400; 0; 126]
α = [-π / 2, π / 2, -π / 2, π / 2, -π / 2, π / 2, 0]

M = [
      7910762//10000679  -4904166//10000679 3658221//10000679 -287711//1000000
      4705866//50003395  34419087//50003395 7193110//10000679  -11941//40000
    -30227313//50003395 -26728166//50003395 5906790//10000679  819169//1000000
            0//1                0//1              0//1              1//1
]

function trutman_pose_constraint(M, d, r, α, c, s)
    T(i) = dh_matrix_rat(c[i], s[i], d[i], α[i], r[i])
    iT(i) = dh_matrix_rat_inverse(c[i], s[i], d[i], α[i], r[i]; tol=1e-5)

    res = T(3)*T(4)*T(5)-iT(2)*iT(1)*M*iT(7)*iT(6)
    view(res, 1:3, :)[:]
end

@polyvar c[eachindex(d)] s[eachindex(d)]
constrs = trutman_pose_constraint(M, d, r, α, c, s)
constrs = [constrs; c .^ 2 .+ s .^ 2 .- 1]

#@polyvar c1 c2 c3 c4 c5 c6 c7 s1 s2 s3 s4 s5 s6 s7
#cg = Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, BigInt}[-39553810 * c6 * c7 * c2 * c1 - 24520830 * c6 * s7 * c2 * c1 - 4705866 * c6 * c7 * c2 * s1 + 34419087 * c6 * s7 * c2 * s1 - 18291105 * s6 * c2 * c1 - 35965550 * s6 * c2 * s1 + 50003395 * c3 * c4 * c5 - 30227313 * c6 * c7 * s2 + 26728166 * c6 * s7 * s2 + 29533950 * s6 * s2 - 50003395 * s3 * s5 -4705866 * c6 * c7 * c1 + 34419087 * c6 * s7 * c1 + 50003395 * s3 * c4 * c5 + 39553810 * c6 * c7 * s1 + 24520830 * c6 * s7 * s1 - 35965550 * s6 * c1 + 50003395 * c3 * s5 + 18291105 * s6 * s1 -39553810 * c6 * c7 * s2 * c1 - 24520830 * c6 * s7 * s2 * c1 - 4705866 * c6 * c7 * s2 * s1 + 34419087 * c6 * s7 * s2 * s1 - 18291105 * s6 * s2 * c1 + 30227313 * c6 * c7 * c2 - 26728166 * c6 * s7 * c2 - 35965550 * s6 * s2 * s1 - 29533950 * s6 * c2 - 50003395 * s4 * c5 -39553810 * s6 * c7 * c2 * c1 - 24520830 * s6 * s7 * c2 * c1 - 4705866 * s6 * c7 * c2 * s1 + 34419087 * s6 * s7 * c2 * s1 + 18291105 * c6 * c2 * c1 + 35965550 * c6 * c2 * s1 - 30227313 * s6 * c7 * s2 + 26728166 * s6 * s7 * s2 - 50003395 * c3 * s4 - 29533950 * c6 * s2 -4705866 * s6 * c7 * c1 + 34419087 * s6 * s7 * c1 + 39553810 * s6 * c7 * s1 + 24520830 * s6 * s7 * s1 + 35965550 * c6 * c1 - 18291105 * c6 * s1 - 50003395 * s3 * s4 -39553810 * s6 * c7 * s2 * c1 - 24520830 * s6 * s7 * s2 * c1 - 4705866 * s6 * c7 * s2 * s1 + 34419087 * s6 * s7 * s2 * s1 + 18291105 * c6 * s2 * c1 + 30227313 * s6 * c7 * c2 - 26728166 * s6 * s7 * c2 + 35965550 * c6 * s2 * s1 + 29533950 * c6 * c2 - 50003395 * c4 24520830 * c7 * c2 * c1 - 39553810 * s7 * c2 * c1 - 34419087 * c7 * c2 * s1 - 4705866 * s7 * c2 * s1 - 50003395 * c3 * c4 * s5 - 50003395 * s3 * c5 - 26728166 * c7 * s2 - 30227313 * s7 * s2 -50003395 * s3 * c4 * s5 - 34419087 * c7 * c1 - 4705866 * s7 * c1 + 50003395 * c3 * c5 - 24520830 * c7 * s1 + 39553810 * s7 * s1 24520830 * c7 * s2 * c1 - 39553810 * s7 * s2 * c1 - 34419087 * c7 * s2 * s1 - 4705866 * s7 * s2 * s1 + 26728166 * c7 * c2 + 30227313 * s7 * c2 + 50003395 * s4 * s5 463813151355769 * c2 * c1 + 909317312698475 * c2 * s1 + 4000271600000000 * c3 * s4 - 4136294153784249 * s2 4000271600000000 * s3 * s4 + 909317312698475 * c1 - 463813151355769 * s1 463813151355769 * s2 * c1 + 909317312698475 * s2 * s1 + 4136294153784249 * c2 + 4000271600000000 * c4 + 4000271600000000 c1 ^ 2 + s1 ^ 2 - 1 c2 ^ 2 + s2 ^ 2 - 1 c3 ^ 2 + s3 ^ 2 - 1 c4 ^ 2 + s4 ^ 2 - 1 c5 ^ 2 + s5 ^ 2 - 1 c6 ^ 2 + s6 ^ 2 - 1 c7 ^ 2 + s7 ^ 2 - 1];

fafa = groebner(constrs)

#then find the degree two polys
#compare to the original by recomputing GB and assert equal
#try with quats