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

using Dates

function parse_args(arg_pose, arg_param, arg_warm, arg_samples)
    sample_count = parse(Int, arg_samples)

    opts_warm = Dict(
        "cold" => no_warm_start,
        "warm" => local_inverse_kinematics
    )

    opts_pose = Dict(
        "uniform" => uniform_pose_gen,
        "feasible" => random_feasible_pose
    )

    opts_param = Dict(
        "kuka" => params_kuka_iiwa(),
        "rand_6rad" => params_random_6rad(7),
        "rand_4rad" => params_random_4rad(7),
        "rand_orth" => params_random_orth(7),
        "icub_v2_7" => params_icub_v2(7),
        "icub_v2_8" => params_icub_v2(8),
        "icub_v2_9" => params_icub_v2(9),
        "icub_v2_10" => params_icub_v2(10)
    )

    pose_gen = opts_pose[arg_pose]
    params = opts_param[arg_param]
    warm_start = opts_warm[arg_warm]

    return pose_gen, params, warm_start, sample_count
end

function run_experiment(_pose, _params, _warm, _samples)
    parsed = parse_args(_pose, _params, _warm, _samples)
    pose_gen, params, warm_start, samples = parsed

    filename_date = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    info_line = "$(_pose)_$(_params)_$(_warm)_$(filename_date)"
    ik_method = solve_inverse_kinematics

    println("$info_line")
    open("DATA_$info_line.txt", "w") do f
        println(f, "# $info_line")
        println(f, "# loc_err rot_err obj local_obj tim ret")
        for _ in 1:samples
            result = stats_sample(ik_method, params...; pose_gen, warm_start)
            println(f, join(result, " "))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_experiment(ARGS...)
end