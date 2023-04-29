using DynamicPolynomials
using MultivariatePolynomials

"""Creates the lifted pose constraint by adding lifting variables as needed."""
function lift_adhoc(d, r, α, M, c, s)
    ids = eachindex(d)
    @polyvar C[ids] S[ids]

    fwd, rev = build_eqs(d, r, α, C, S)

    chain_poly_dirty = prod(fwd) .- M * prod(rev)
    chain_poly_clean = mapcoefficients.(round_zero, chain_poly_dirty)
    chain_jump = map(e -> e([C; S] => [c; s]), chain_poly_clean)

    chain_jump
end