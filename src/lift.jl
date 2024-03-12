"""Creates the lifted pose constraint by multiplying matrices and adding lifting 
variables as needed."""
function lift(d, r, α, M::AbstractArray, c, s)
    F, R = build_pose_constraint(d, r, α, c, s)

    LHS = prod(F)
    RHS = prod(R)

    view(LHS - M * RHS, 1:3, :)
end

function lift(d, r, α, M::DualQuaternion, c, s)
    F, R = build_pose_constraint_q(d, r, α, c, s)

    vec(prod(F) - M * prod(R))
end