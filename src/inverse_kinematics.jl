using PolyJuMP
using JuMP
using Gurobi

function _default_optimizer()
    grb = Gurobi.Optimizer()
    #MOI.set(grb, MOI.Silent(), true)
    #MOI.set(grb, MOI.RawOptimizerAttribute("Nonconvex"), 2)
    #MOI.set(grb, MOI.RawOptimizerAttribute("Presolve"), 2)
    #MOI.set(grb, MOI.RawOptimizerAttribute("Threads"), 4)
    () -> PolyJuMP.QCQP.Optimizer(grb)
end

function _extract_solution(c, s, m)
    stat = termination_status(m)

    vs = has_values(m) ? clamp.(value.(s), -1, 1) : fill(NaN, length(c))
    vc = has_values(m) ? clamp.(value.(c), -1, 1) : fill(NaN, length(c))

    sol = atan.(vs, vc)
    obj = stat == OPTIMAL ? objective_value(m) : NaN

    sol, obj, stat, solve_time(m)
end

"""
    solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
    lift_method=lift_matrix, optimizer=_default_optimizer(), init=θ)

Computes the global inverse kinematics solution using a chosen lift_method, 
starting from `init`.
"""
function solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w; 
    optimizer=_default_optimizer(), init=θ)

    ids = eachindex(d)
    m = Model(optimizer)

    @variable(m, c[ids])
    @variable(m, s[ids])
    constrain_trig_vars.(c, s, θl, θh, init)

    @constraint m lift(d, r, α, M, c, s) .== 0
    @constraint m c .^ 2 .+ s .^ 2 .== 1
    @constraint m lin_angdiff_proxy.(c, s, θl) .>= 0
    @constraint m lin_angdiff_proxy.(c, s, θh) .<= 0

    @objective m Min sum(lin_abs_angdiff_proxy.(c, s, θ, w))
    optimize!(m)

    #_extract_solution(c, s, m)
end
