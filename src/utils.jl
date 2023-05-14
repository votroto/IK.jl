using Base.Iterators: take, drop

"""Rounds tiny values to zero"""
round_zero(x; atol=1e-12) = abs(x) <= atol ? zero(x) : x

"""Performs a tree-like fold over a magma `f`"""
function _reduce(f, xs; init=1)
    n = length(xs)

    if n == 0
        init
    elseif n == 1
        first(xs)
    else
        mid = div(length(xs), 2)
        h, t = take(xs, mid), drop(xs, mid)
        f(_reduce(f, h; init), _reduce(f, t; init))
    end
end

"""Finds the minimum and maximum of cos(x), for l <= x <= h"""
function cos_min_max(l, h)
    if h - l >= 2pi
        return -1.0, 1.0
    end

    ql = 0.5pi * ceil(l / (0.5pi))
    crit = [ql:0.5pi:h; l; h]
    minimum(cos, crit), maximum(cos, crit)
end

function sin_min_max(l, h)
    cos_min_max(l - pi / 2, h - pi / 2)
end

function replicate(xs, ns)
    mapreduce(fill, vcat, xs, ns)
end