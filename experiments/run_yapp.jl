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

function parse_args(arg_pose, arg_param, arg_samples)
    sample_count = parse(Int, arg_samples)

    opts_pose = Dict(
        "uniform" => uniform_pose_gen,
        "feasible" => random_feasible_pose,
        "kuka_bounded" => random_bounded_pose(-0.8, 0.8, -0.8, 0.8, 0.0, 1.0),
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

    return pose_gen, params, sample_count
end

function run_experiment(_pose, _params, _samples)
    parsed = parse_args(_pose, _params, _samples)
    pose_gen, params, samples = parsed

    warm_start = local_inverse_kinematics
    filename_date = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    info_line = "$(_pose)_$(_params)_$(filename_date)"
    ik_method = solve_inverse_kinematics

    println("$info_line")
    open("CLOUD_$info_line.txt", "w") do f
        println(f, "# $info_line")
        println(f, "# x, y, z, loc_err, obj, ret, tim, local_obj")
        for sample_i in 1:samples
            @show sample_i
            result = yapp_sample(ik_method, params...; pose_gen, warm_start)
            println(f, join(result, " "))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_experiment(ARGS...)
end