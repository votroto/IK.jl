using DynamicPolynomials
using MultivariatePolynomials

"""Creates the lifted pose constraint by adding lifting variables as needed."""
function lift_adhoc(d, r, α, M, c, s)
    ids = eachindex(d)
    @polyvar C[ids] S[ids]

    E = build_eqs_poly(d, r, α, C, S, M)
    lifted = map(e -> e([C; S] => [c; s]), E)

    lifted
end