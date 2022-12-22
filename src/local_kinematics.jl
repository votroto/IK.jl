using Base.Iterators: peel, drop, take
using HomotopyContinuation

include("denavit_hartenberg.jl")

function _split_manipulator(ids)
        mid = div(length(ids), 2)
        f, s = take(ids, mid), drop(ids, mid)
        f, reverse(collect(s))
end

function build_obj(x, θ, d, r, α, M, w)
        T(i) = dh_t(x[i], d[i], α[i], r[i])
        iT(i) = dh_inv_t(x[i], d[i], α[i], r[i])

        ids = eachindex(d)
        fwd, rev = _split_manipulator(ids)

        pp = (prod(T, fwd).-M*prod(iT, rev))[1:3, 1:4][:]
        return sum(pp .^ 2)
end

function local_inverse_kinematics(d, r, α, θl, θh, M, θ, w)
        lb = [θl...]
        ub = [θh...]

        obj(x) = build_obj(x, θ, d, r, α, M, w)
        nlp = ADNLPModel(obj, θ,
                x -> x,
                lb, ub)
        stats = ipopt(nlp)
	stats.solution, stats.objective
end