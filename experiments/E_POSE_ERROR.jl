include("./benchmark_ik_method.jl")

include("./kuka_parameters.jl")

include("../src/inverse_kinematics_matrix.jl")
include("../src/local_kinematics.jl")

ik_method = first ∘ solve_inverse_kinematics_matrix
params = params_kuka_iiwa()
warm_start = first ∘ local_inverse_kinematics

pose_error_sample(ik_method, params...; warm_start)