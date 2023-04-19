using Base.Iterators: peel, drop, take
using HomotopyContinuation

include("denavit_hartenberg.jl")

function _split_manipulator(ids)
        mid = div(length(ids), 2)
        f, s = take(ids, mid), drop(ids, mid)
        f, reverse(collect(s))
end

function build_eqs(c, s, d, r, α, M)
        T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
        iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

        ids = eachindex(d)
        fwd, rev = _split_manipulator(ids)

        prod(T, fwd) .- M * prod(iT, rev)
end

function homo_inverse_kinematics(d, r, α, θl, θh, M, θ, w)
        ids = eachindex(d)
        @var c[ids] s[ids]

        C = [c .^ 2 .+ s .^ 2 .- 1]
        E = build_eqs(c, s, d, r, α, M)
	sys = System([E[1:4,1:3][:]; C...])
	solve(sys)


end