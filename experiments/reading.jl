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


@polyvar c[1:7] s[1:7]

function zzz(i)
    file_path = "simple/con$i.out"
    file_content = open(file_path, "r") do file
        read(file, String)
    end

    readd = eval(Meta.parse(strip(file_content)[1:end-1]))
    divided = [DynamicPolynomials.map_coefficients(g->Float64(g), p/leading_coefficient(p)) for p in readd]

    @show file_path
    
    local_x, local_obj = local_inverse_kinematicsgb([c;s], divided, fill(-3.0, 7), fill(3.0,7), zeros(7), ones(7)/7)
    x, obj, ret, tim = gb_inverse_kinematics_an(divided, c, s, fill(-3.0, 7), fill(3.0,7), zeros(7), ones(7)/7, init=local_x)
    
    @show local_x
    @show x

    @show [d(c=>cos.(local_x), s=>sin.(local_x)) for d in divided]
    @show [d(c=>cos.(x), s=>sin.(x)) for d in divided]

    loc_err, rot_err = NaN, NaN
    result = loc_err, rot_err, obj, local_obj, tim, ret
end

info_line = "SCIP2"
println("$info_line")
open("/tmp/DATA_$info_line.txt", "w") do f
    println(f, "# $info_line")
    println(f, "# loc_err rot_err obj local_obj tim ret")
    for i in 1:100
        result = zzz(i)
        println(f, join(result, " "))
    end
end