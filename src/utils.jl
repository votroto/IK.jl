using Base.Iterators: take, drop

# Thanks @ivirshup!
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

round_zero(x; atol=1e-12) = abs(x) <= atol ? 0 : x

function _reduce(f, xs; init=one(eltype(xs)))
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

function cos_min_max(l, h)
        if h - l >= 2pi
                return -1.0, 1.0
        end

        _l = round(l / (0.5pi)) - 1
        _h = round(h / (0.5pi)) + 1

        crit = [0.5pi * i for i in _l:_h]
        filter!(x -> l <= x <= h, crit)
        vals = [crit; l; h]
        minimum(cos, vals), maximum(cos, vals)
end

function sin_min_max(l, h)
        cos_min_max(l - pi / 2, h - pi / 2)
end

function replicate(xs, ns)
        mapreduce(fill, vcat, xs, ns)
end