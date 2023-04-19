using JuMP
#using DynamicPolynomials
#using SumOfSquares
using MosekTools
#using SemialgebraicSets

m = Model(Mosek.Optimizer)
@variable m c
@variable m f
@variable m g
e = 1 - 2c
d = 0.5 - f
h = 1 - 2g
Q = [
        1 0.5 c d
        0.5 e f g
        c f h 0.5
        d g 0.5 1
]

@constraint(m, Q >= 0, PSDCone())


@objective m Min c
optimize!(m)
lc = objective_value(m)

@objective m Max c
optimize!(m)
hc = objective_value(m)


@objective m Min f
optimize!(m)
lf = objective_value(m)

@objective m Max f
optimize!(m)
hf = objective_value(m)


@objective m Min g
optimize!(m)
lg = objective_value(m)

@objective m Max g
optimize!(m)
hg = objective_value(m)

function ub(c, f)

        m = Model(Mosek.Optimizer)
        @variable m g
        e = 1 - 2c
        d = 0.5 - f
        h = 1 - 2g
        Q = [
                1 0.5 c d
                0.5 e f g
                c f h 0.5
                d g 0.5 1
        ]

        @constraint(m, Q >= 0, PSDCone())


        @objective m Max g
        set_silent(m)
        optimize!(m)
        if termination_status(m) == MOI.OPTIMAL
                objective_value(m)
        else
                NaN
        end
end


function lb(c, f)

        m = Model(Mosek.Optimizer)
        @variable m g
        e = 1 - 2c
        d = 0.5 - f
        h = 1 - 2g
        Q = [
                1 0.5 c d
                0.5 e f g
                c f h 0.5
                d g 0.5 1
        ]

        @constraint(m, Q >= 0, PSDCone())


        @objective m Min g
        set_silent(m)
        optimize!(m)
        if termination_status(m) == MOI.OPTIMAL
                objective_value(m)
        else
                NaN
        end
end



function ubf(c,g)

        m = Model(Mosek.Optimizer)
        @variable m f
        e = 1 - 2c
        d = 0.5 - f
        h = 1 - 2g
        Q = [
                1 0.5 c d
                0.5 e f g
                c f h 0.5
                d g 0.5 1
        ]

        @constraint(m, Q >= 0, PSDCone())


        @objective m Max f
        set_silent(m)
        optimize!(m)
        if termination_status(m) == MOI.OPTIMAL
                objective_value(m)
        else
                NaN
        end
end




function lbf(c,g)

        m = Model(Mosek.Optimizer)
        @variable m f
        e = 1 - 2c
        d = 0.5 - f
        h = 1 - 2g
        Q = [
                1 0.5 c d
                0.5 e f g
                c f h 0.5
                d g 0.5 1
        ]

        @constraint(m, Q >= 0, PSDCone())


        @objective m Min f
        set_silent(m)
        optimize!(m)
        if termination_status(m) == MOI.OPTIMAL
                objective_value(m)
        else
                NaN
        end
end


function lbc(f,g)

        m = Model(Mosek.Optimizer)
        @variable m c
        e = 1 - 2c
        d = 0.5 - f
        h = 1 - 2g
        Q = [
                1 0.5 c d
                0.5 e f g
                c f h 0.5
                d g 0.5 1
        ]

        @constraint(m, Q >= 0, PSDCone())


        @objective m Min c
        set_silent(m)
        optimize!(m)
        if termination_status(m) == MOI.OPTIMAL
                objective_value(m)
        else
                NaN
        end
end


open("chul.dat", "w") do f
	la = NaN
        for cc in lc:0.05:hc
                for ff in lf:0.05:hf
                        u = ub(cc, ff)
                        println(f, "$cc $ff $u")
                end
                println(f)
        end


        for cc in lc:0.05:hc
                for ff in lf:0.05:hf
                        u = lb(cc, ff)
                        println(f, "$cc $ff $u")
                end
                println(f)
        end


        for cc in lc:0.05:hc
                for gg in lg:0.05:hg
                        u = ubf(cc, gg)
                        println(f, "$cc $u $gg")
                end
                println(f)
        end



        for cc in lc:0.05:hc
                for gg in lg:0.05:hg
                        u = lbf(cc, gg)
                        println(f, "$cc $u $gg")
                end
                println(f)
        end


        for ff in lf:0.05:hf
                for gg in lg:0.05:hg
                        u = lbc(ff, gg)
                        println(f, "$u $ff $gg")
                end
                println(f)
        end
end


#=
include("replication_parameters.jl")

include("../src/forward_kinematics.jl")
include("../src/utils.jl")




function map_monomials(f, poly)
        sum(coefficients(poly) .* map(f, monomials(poly)), init=0)
end


#=
d offset
r radius
α twist
M desired pose
θ initial angle
w angle weights
=#

function lmo!(m, ine, mon)
        @show dm = degree(mon)
        @show mon
        function fefe(a, b)
                L, = @polyvar L[rand(UInt32)]
                @show ce = L[1] - a * b
                @show ce = -a * b
                push!(m, ce)
                push!(ine, -L[1]^2 - 1)
                L[1]

        end
        if dm < 3
                return mon
        end
        if dm == 3
                a, b, c = variables(mon)
                ab = fefe(a, b)
                return ab * c
        end
        if dm == 4
                a, b, c, d = variables(mon)
                ab = fefe(a, b)
                cd = fefe(c, d)
                return ab * cd
        end
        pre = replicate(variables(mon), exponents(mon))
        reduce(fefe, pre; init=1)
end

function lva!(m, ine, eqs)
        nnn = Dict()
        mo = unique(sort([mon for q in eqs for mon in monomials(q)]))

        for mon in mo
                v = lmo!(m, ine, mon)
                nnn[mon] = v
        end

        nnn
end



function build_eqs(d, r, α, M)
        T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
        iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

        ids = eachindex(d)
        fwd, rev = _split_manipulator(ids)

        @polyvar c[ids] s[ids]

        prod(T, fwd) .- M * prod(iT, rev), c, s
end

function sos(d, r, α, θl, θh, M, θ, w; de, aa=false)

        T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
        iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

        ids = eachindex(d)
        fwd, rev = _split_manipulator(ids)
        lids = length(ids)
        @polyvar c[1:lids] s[1:lids]



        circ = c .^ 2 .+ s .^ 2 .- 1
        bndl = -((c .+ 1) .* tan.(θl / 2) .- s)
        bndh = (c .+ 1) .* tan.(θh / 2) .- s
        cc = c .^ 2 .- 1
        ss = s .^ 2 .- 1
        ineqs = [bndl; bndh; cc; ss]


        E = prod(T, fwd) .- M * prod(iT, rev)
        G = mapcoefficients.(c -> (abs(c) > (1e-12)) ? c : 0.0, E)

        if aa
                nn = lva!(circ, ineqs, G)
                @show H = vcat(map_monomials.(x -> nn[x], G[1:3, :]), G[4, :]')
        else
                H = G
        end
        @show eqs = [H[:]; circ]

        @show cons = basicsemialgebraicset(algebraicset(eqs), ineqs)

        @show obj = sum(2 .* w .* (1 .- c .* cos.(θ) .- s .* sin.(θ)))

        m = SOSModel(Mosek.Optimizer)
        @variable m a
        @constraint m mo obj >= a domain = cons maxdegree = de
        @objective m Max a

        optimize!(m)

        ν3 = moment_matrix(mo)
        extractatoms(ν3, 1e-3)
end


njoints = 5
r, d, α, θl, θh = random_manipulator(njoints)
w = normalize(ones(njoints), 1)
θ = zeros(njoints)
M = random_feasible_pose(d, r, α, θl, θh)
sos(d, r, α, θl, θh, M, θ, w; de=2)
=#