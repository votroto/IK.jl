using JuMP
using Gurobi

function _default_optimizer()
    optimizer_with_attributes(
        Gurobi.Optimizer,
        #MOI.Silent() => true,
        "SolutionLimit" => 1,
        "FuncNonlinear" => 1,
        "Nonconvex" => 2,
        "Presolve" => 2,
        "Threads" => 4
    )
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
    lift_method=lift_matrix, optimizer=_default_optimizer(), init=θ)

    ids = eachindex(d)
    m = Model(optimizer)

    @variable(m, c[ids])
    @variable(m, s[ids])
    constrain_trig_vars.(c, s, θl, θh, init)

    @constraint m lift_method(d, r, α, M, c, s) .== 0
    @constraint m c .^ 2 .+ s .^ 2 .== 1
    @constraint m lin_angdiff_proxy.(c, s, θl) .>= 0
    @constraint m lin_angdiff_proxy.(c, s, θh) .<= 0

    @objective m Min sum(w .* lin_abs_angdiff_proxy.(c, s, θ))
    optimize!(m)

    _extract_solution(c, s, m)
end

column(x::VariableRef) = Gurobi.c_column(backend(owner_model(x)), index(x))

function gb_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
    lift_method=lift_matrix, init=θ)

    optimizer = _default_optimizer()

    ids = eachindex(θ)
    m = direct_model(optimizer)

    xstrt = init
    ystrt = abs.(mod2pi.(init .- θ .+ 7*pi) .- pi) 
    cstrt = cos.(init)
    sstrt = sin.(init)
    @variable(m, θl[i] <= x[i in ids] <= θh[i], start = xstrt[i])
    @variable(m, 0 <= y[i in ids] <= pi, start = ystrt[i])
    @variable(m, -1 <= c[i in ids] <= 1, start = cstrt[i])
    @variable(m, -1 <= s[i in ids] <= 1, start = sstrt[i])
    
    constrain_trig_vars.(c, s, θl, θh, init)

    for i in ids
        GRBaddgenconstrCos(backend(m), C_NULL, column(x[i]), column(c[i]), C_NULL)
        GRBaddgenconstrSin(backend(m), C_NULL, column(x[i]), column(s[i]), C_NULL)
        
        GRBaddgenconstrPWL(backend(m), C_NULL, column(x[i]), column(y[i]), 5, collect(θ[i]-2pi:pi:θ[i]+2pi), [0.,pi,0.,pi,0.])
    end
    
    @constraint m lift_method(d, r, α, M, c, s) .== 0
    @objective m Min sum(w .* y) 
    
    optimize!(m)

    _extract_solution(c, s, m)    
end

