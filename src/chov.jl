using DynamicPolynomials
using LinearAlgebra

@polyvar x[1:3]

mms = monomials(x, 1:3)
ccs1 = rand(length(mms))
ccs2 = rand(length(mms))

p1 = dot(mms, ccs1)
p2 = dot(mms, ccs2)

eqs = [p1, p2]

ums = filter(x-> maxdegree(x) > 2, unique(reduce(vcat, monomials.(eqs))))
@polyvar _lift[eachindex(ums)]

d = Dict(ums .=> _lift)

ceqs = map_monomials.(m -> get(d, m, m), eqs)