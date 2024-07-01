using PolyJuMP
using JuMP
using SCIP
using Gurobi


function _default_optimizer()
    grb = Gurobi.Optimizer()

    MOI.set(grb, MOI.RawOptimizerAttribute("Threads"), 1)
    #MOI.set(grb, MOI.RawOptimizerAttribute("SolutionLimit"), 1)
    #MOI.set(grb, MOI.RawOptimizerAttribute("Presolve"), 0)
    #MOI.set(grb, MOI.RawOptimizerAttribute("Heuristics"), 0)
    #MOI.set(grb, MOI.RawOptimizerAttribute("NLPHeur"), 0)
    #MOI.set(grb, MOI.RawOptimizerAttribute("Cuts"), 0)
    
    () -> PolyJuMP.QCQP.Optimizer(grb)
end


function _gb_optimizer()
    optimizer_with_attributes(SCIP.Optimizer,
        "parallel/maxnthreads" => 4,
        "limits/time" => 100
    )
    #optimizer_with_attributes(Gurobi.Optimizer,
    #    "Threads" => 4,
    #)
end

recover_angle(::AbstractArray, s, c) = atan(s, c)
recover_angle(::DualQuaternion, s, c) = atan(s, c) * 2

function _extract_solution(M, c, s, m)
    stat = termination_status(m)

    vs = has_values(m) ? clamp.(value.(s), -1, 1) : fill(NaN, length(c))
    vc = has_values(m) ? clamp.(value.(c), -1, 1) : fill(NaN, length(c))

    sol = recover_angle.(Ref(M), vs, vc)
    obj = stat == OPTIMAL ? objective_value(m) : NaN

    sol, obj, stat, solve_time(m)
end



function _extract_solution(x, m)
    stat = termination_status(m)

    sol = has_values(m) ? value.(x) : fill(NaN, length(c))    
    obj = stat == OPTIMAL ? objective_value(m) : NaN

    sol, obj, stat, solve_time(m)
end


"""
    solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
    lift_method=lift_matrix, optimizer=_default_optimizer(), init=θ)

Computes the global inverse kinematics solution using a chosen lift_method, 
starting from `init`.
"""

column(x::VariableRef) = Gurobi.c_column(backend(owner_model(x)), index(x))

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

    #@objective m Min sum(lin_abs_angdiff_proxy.(c, s, θ, w))
    optimize!(m)

    _extract_solution(M, c, s, m)
end

function gb_inverse_kinematics(eqs, C, S, θl, θh, θ, w; init=θ)
    optimizer = optimizer_with_attributes(
        Gurobi.Optimizer,
        "FuncNonlinear" => 1
    )

    ids = eachindex(θ)
    m = direct_model(optimizer)

    @variable(m, θl[i] <= x[i in ids] <= θh[i])
    @variable(m, 0 <= y[i in ids] <= 1)

    @variable(m, -1 <= c[ids] <= 1)
    @variable(m, -1 <= s[ids] <= 1)

    
    for i in ids
        GRBaddgenconstrCos(backend(m), C_NULL, column(x[i]), column(c[i]), C_NULL)
        GRBaddgenconstrSin(backend(m), C_NULL, column(x[i]), column(s[i]), C_NULL)

        GRBaddgenconstrPWL(backend(m), C_NULL, column(x[i]), column(y[i]), 5, collect(θ[i]-2pi:pi:θ[i]+2pi), [0.,1.,0.,1.,0.])
    end

#    constrain_trig_vars.(c, s, θl, θh, init)

    jump_eqs = [e([C; S] => [c; s]) for e in eqs]
    @constraint m jump_eqs .== 0
    #@constraint m lin_angdiff_proxy.(c, s, θl) .>= 0
    #@constraint m lin_angdiff_proxy.(c, s, θh) .<= 0

    @objective m Min sum(y) #lin_abs_angdiff_proxy.(c, s, θ, w))
    optimize!(m)

    _extract_solution(x, m)
end





function gb_inverse_kinematics_an(eqs, C, S, θl, θh, θ, w; init=θ)
    optimizer=_gb_optimizer()

    ids = eachindex(C)
    m = Model(optimizer)

    @variable(m, c[ids])
    @variable(m, s[ids])
    constrain_trig_vars.(c, s, θl, θh, init)

    @constraint m [e([C;S]=>[c;s]) for e in eqs] .== 0
    @constraint m c .^ 2 .+ s .^ 2 .== 1
    @constraint m lin_angdiff_proxy.(c, s, θl) .>= 0
    @constraint m lin_angdiff_proxy.(c, s, θh) .<= 0

    @objective m Min sum(lin_abs_angdiff_proxy.(c, s, θ, w))
    optimize!(m)

    _extract_solution([], c, s, m)
end


function gb_inverse_kinematicsold(eqs, C, S, θl, θh, θ, w; init=θ)
    optimizer = optimizer_with_attributes(
        Gurobi.Optimizer
    )

    ids = eachindex(θ)
    m = Model(optimizer)


    @variable(m, c[ids])
    @variable(m, s[ids])


    constrain_trig_vars.(c, s, θl, θh, init)

    jump_eqs = [e([C; S] => [c; s]) for e in eqs]
    @constraint m jump_eqs .== 0
    @constraint m lin_angdiff_proxy.(c, s, θl) .>= 0
    @constraint m lin_angdiff_proxy.(c, s, θh) .<= 0

    @objective m Min sum(lin_abs_angdiff_proxy.(c, s, θ, w))
    optimize!(m)

    _extract_solution([], c, s, m)
end