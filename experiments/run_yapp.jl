include("./benchmark_ik_method.jl")
include("./yapp.jl")

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

using Dates

function parse_args(arg_samples)
    sample_count = parse(Int, arg_samples)

    pose_gen = random_bounded_pose(-0.8, 0.8, -0.8, 0.8, 0.0, 1.0)
    params = params_kuka_iiwa()
    warm_start = local_inverse_kinematics

    return pose_gen, params, warm_start, sample_count
end

function run_experiment(pose_gen, params, warm_start, samples)
    filename_date = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    info_line = "kuka_yapp_$(filename_date)"
    ik_method = solve_inverse_kinematics

    println("$info_line")
    open("DATA_$info_line.txt", "w") do f
        println(f, "# $info_line")
        println(f, "# x, y, z, loc_err, obj, ret")
        for sample_i in 1:samples
            @show sample_i
            result = yapp_sample(ik_method, params...; pose_gen, warm_start)
            println(f, join(result, " "))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_experiment(parse_args(ARGS...)...)
end