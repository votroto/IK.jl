include("./denavit_hartenberg.jl")

using Base.Iterators: peel, drop, take
using ADNLPModels
using Ipopt
using NLPModelsIpopt

norm2diff(x, y) = (x - y) * (x - y)

function build_K(d, r, α, c, s)
    T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
    prod(T, eachindex(c))
end

function solve_position_priority(d, r, α, θl, θh, M, θ, w;
    lift_method=lift_matrix, optimizer=_default_optimizer(), init=θ)

    ids = eachindex(d)
    m = Model(optimizer)

    @variable(m, c[ids])
    @variable(m, s[ids])
    constrain_trig_vars.(c, s, θl, θh, init)

   @show  K = build_K(d, r, α, c, s)
    map_coefficients_inplace!.(round_zero, K[1:3,1:4])

    Ma = M[1:3, 1:3]
    Ka = K[1:3, 1:3]
    Mx = M[1:3, 4]
    Kx = K[1:3, 4]

    @constraint m Mx .== Kx
    @constraint m c .^ 2 .+ s .^ 2 .== 1
    @constraint m lin_angdiff_proxy.(c, s, θl) .>= 0
    @constraint m lin_angdiff_proxy.(c, s, θh) .<= 0

    @objective m Min mapreduce(norm2diff, +, Ma, Ka)
    optimize!(m)

@show 'a'
    _extract_solution(c, s, m)
end