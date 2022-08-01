using LinearAlgebra
using Gurobi
using JuMP

function _default_optimizer()
    optimizer_with_attributes(Gurobi.Optimizer, "Nonconvex" => 2)
end

"""
d offset
r radius
α twist
M desired pose
θ initial angle
w angle weights
"""
function solve_inverse_kinematics(d, r, α, M, θ, w; optimizer=_default_optimizer())
    T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
    iT(i) = dh_inv_t(c[i], s[i], d[i], α[i], r[i])

    m = Model(optimizer)

    @variable m -1 <= c[1:7] <= 1
    @variable m -1 <= s[1:7] <= 1
    @variable m x[1:4, 1:4, 1:3]

    @constraint m x[4, 1:3, :] .== 0
    @constraint m x[4, 4, :] .== 1
    @constraint m T(3) * T(4) .== x[:, :, 1]
    @constraint m x[:, :, 1] * T(5) .== x[:, :, 2]
    @constraint m x[:, :, 2] * T(6) .== x[:, :, 3]
    @constraint m x[:, :, 3] * T(7) .== iT(2) * iT(1) * M

    @constraint m c .^ 2 + s .^ 2 .== 1
    @constraint m (c .+ 1) .* tan.(θl ./ 2) .- s .<= 0
    @constraint m (c .+ 1) .* tan.(θh ./ 2) .- s .>= 0

    @objective m Min sum(2 .* w .* (1 .- c .* cos.(θ) .- s .* sin.(θ)))

    optimize!(m)

    status = termination_status(m)

    sol = if status == MOI.OPTIMAL
        atan.(value.(s), value.(c))
    else
        undefined
    end

    sol, status
end
