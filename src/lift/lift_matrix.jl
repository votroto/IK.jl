"""Creates the lifted pose constraint by multiplying matrices and adding lifting 
variables as needed."""
function lift_matrix(d, r, α, M, c, s)
    F, R = build_pose_constraint(d, r, α, c, s)

    LHS = prod(F)
    RHS = prod(R)

    view(LHS - M * RHS, 1:3, :)
end

"""Creates the lifted pose constraint by multiplying matrices and adding lifting 
variables as needed."""
function lift_matrixq(d, r, α, M, c, s)
    F, R = build_pose_constraintq(d, r, α, c, s)

    LHS = prod(F)
    RHS = prod(R)

    view(LHS - M * RHS, 1:3, :)
end