using Base.Iterators:take,drop

function _split_manipulator(ids)
        mid = div(length(ids), 2)
        f, s = take(ids, mid), drop(ids, mid)
        f, reverse(collect(s))
end

function _reduce(f, xs)
        mid = div(length(xs), 2)
        h, t = take(xs, mid), drop(xs, mid)
        f(reduce(f, h), reduce(f, t))
end


function cos_min_max(l, h)
        l, h = sort([l, h])
        rl, rh = mod(l, 2 * pi), mod(h, 2 * pi)
        rh = (rh > rl) ? rh : rh + 2 * pi

        inflcos = filter(x -> rl <= x <= rh, [pi * i / 2 for i in 0:4])
        infl = map(cos, [inflcos; rl; rh])
        minimum(infl), maximum(infl)
end

function sin_min_max(l, h)
        cos_min_max(l - pi / 2, h - pi / 2)
end


function replicate(xs, ns)
        nxs = similar(xs, sum(ns))
        ni = 1
        for i in eachindex(xs)
                for j in 1:ns[i]
                        nxs[ni] = xs[i]
                        ni += 1
                end
        end
        nxs
end

function preduce(f, xs; init=1)
        mid = div(length(xs), 2)
        h, t = take(xs, mid), drop(xs, mid)
        f(reduce(f, h; init), reduce(f, t; init))
end

