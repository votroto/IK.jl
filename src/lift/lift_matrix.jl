"""Creates the lifted pose constraint by multiplying matrices and adding lifting 
variables as needed."""
function lift_matrix(d, r, α, M, c, s)
    F, R = build_eqs(d, r, α, c, s)

    LHS = _reduce(jump_quadratic_product, F)
    RHS = _reduce(jump_quadratic_product, R)

    view(LHS - M * RHS, 1:3, :)
end