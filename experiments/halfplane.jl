using Symbolics

function halfplane(lower, inner, upper, x, y)
    lc, ls = cos(lower), sin(lower)
    ic, is = cos(inner), sin(inner)
    uc, us = cos(upper), sin(upper)


    (uc - lc)*(x - lc) + (us - ls)*(y - ls) >~ 0
end

@variables x y

lo = -pi/2
in = 0
up = pi/2

halfplane(lo, in, up, x, y)