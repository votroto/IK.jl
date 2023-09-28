module IK

export solve_forward_kinematics
export solve_inverse_kinematics
export local_inverse_kinematics

include("utils.jl")
include("modelling.jl")
include("jump_extensions.jl")
include("denavit_hartenberg.jl")
include("forward_kinematics.jl")
include("lift/lift_matrix.jl")
include("local_kinematics.jl")
include("inverse_kinematics.jl")

include("solve_position_priority.jl")

end
