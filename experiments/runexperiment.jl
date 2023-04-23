include("./benchmark_ik_method.jl")

include("./kuka_parameters.jl")
include("./rand_parameters.jl")
include("./icub_arm_parameters.jl")

include("../src/inverse_kinematics.jl")
include("../src/local_kinematics.jl")

using Dates

experiment = ARGS[1]
sample_count = parse(Int, ARGS[2])

warm_start = local_inverse_kinematics
ik_method = solve_inverse_kinematics

if experiment == "kuka"
    params = params_kuka_iiwa()
elseif experiment == "rand_6rad"
    params = params_random_6rad(7)
elseif experiment == "rand_4rad_cold"
    warm_start = no_warm_start
    params = params_random_4rad(7)
elseif experiment == "rand_4rad"
    params = params_random_4rad(7)
elseif experiment == "rand_orth"
    params = params_random_orth(7)
elseif experiment == "icub_v2_7"
    params = params_icub_v2(7)
elseif experiment == "icub_v2_8"
    params = params_icub_v2(8)
elseif experiment == "icub_v2_9"
    params = params_icub_v2(9)
elseif experiment == "icub_v2_10"
    params = params_icub_v2(10)
end

filename_date = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")

println("$(experiment) $sample_count $(filename_date)")
open("DATA_bench_$(experiment)_$(filename_date).txt", "w") do f

    println(f, "# $(experiment) $sample_count $(filename_date)")
    println(f, "# loc_err rot_err obj local_obj tim ret")
    for i in 1:sample_count
        result = stats_sample(ik_method, params...; warm_start)
        println(f, join(result, " "))
    end
end
