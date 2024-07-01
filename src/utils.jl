using Base.Iterators: take, drop

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

function rand_between(l, h)
    l + rand() * (h - l)
end

function rand_between(l, h, n)
    ntuple(_ -> rand_between(l, h), n)
end