using DynamicPolynomials
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

joint_count = 7
@polyvar c[1:7] s[1:7]

θl = vec(fill(-3, joint_count))
θh = vec(fill(+3, joint_count))

w = normalize(ones(joint_count), 1)
θ = zeros(joint_count)

function zzz(i, f)
    file_path = "simple/con$i.out"
    file_content = open(file_path, "r") do file
        read(file, String)
    end

    readd = eval(Meta.parse(strip(file_content)[1:end-1]))
    divided = [DynamicPolynomials.map_coefficients(g->Float64(g), p/leading_coefficient(p)) for p in readd]

    ids = 1:7
    
    println(f, "set IDS;")
    println(f, "var c{IDS} <= 1, >= -1;")
    println(f, "var s{IDS} <= 1, >= -1;")

    for i in eachindex(divided)
        println(f, "s.t. eq$i: $(divided[i]) = 0;")
    end

    for i in ids
        println(f, "s.t. ub$i: $(lin_angdiff_proxy(c[i], s[i], θl[i])) >= 0;")
        println(f, "s.t. lb$i: $(lin_angdiff_proxy(c[i], s[i], θh[i])) <= 0;")
    end
    
    println(f, "minimize move: $(sum(lin_abs_angdiff_proxy.(c, s, θ, w)));")

    println(f, "data;")
    println(f, "set IDS := 1 2 3 4 5 6 7;")
    
end


info_line = "SCIP2"
for i in 1:100
    open("ampl/mod_$i.txt", "w") do f
        zzz(i, f)
    end
end