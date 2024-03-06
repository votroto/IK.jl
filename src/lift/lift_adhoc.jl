using DynamicPolynomials
using MultivariatePolynomials

"""Creates the lifted pose constraint by adding lifting variables as needed."""
function lift_adhoc(d, r, α, M, c, s)
    fwd, rev = build_pose_constraint(d, r, α, c, s)
    prod(fwd) - M * prod(rev)

    #ids = eachindex(d)
    #@polyvar C[ids] S[ids]
#
    #E = build_pose_constraint_poly(d, r, α, C, S, M)[1:3,:]
    #@show E
    #lifted = map(e -> e([C; S] => [c; s]), E)
    #println.(lifted)
    #lifted
end