using JuMP
using Gurobi

const GRB_ENV_REF = Ref{Gurobi.Env}()

function __init__()
    global GRB_ENV_REF
    GRB_ENV_REF[] = Gurobi.Env()
    return
end

function _default_optimizer()
    global GRB_ENV_REF
    optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV_REF[]),
        #MOI.Silent() => true,
        "BestObjStop" => 0.0,
        "FuncNonlinear" => 1,
        "Nonconvex" => 2,
        "Presolve" => 2,
        "Threads" => 4
    )
end

function _extract_solution(x, m)
    stat = termination_status(m)

    sol = has_values(m) ? value.(x) : fill(NaN, length(x))
    obj = stat == OPTIMAL ? objective_value(m) : NaN

    sol, obj, stat
end

col(x::VariableRef) = Gurobi.c_column(backend(owner_model(x)), index(x))

"""
    solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
    optimizer=_default_optimizer(), init=θ)

Computes the global inverse kinematics solution using a chosen lift_method,
starting from `init`.
"""
function solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
    optimizer=_default_optimizer(), init=θ)

    ids = eachindex(θ)
    m = direct_model(optimizer)

    @variable(m, θl[i] <= x[i in ids] <= θh[i], start = init[i])
    @variable(m, -1 <= c[i in ids] <= 1, start = cos(init[i]))
    @variable(m, -1 <= s[i in ids] <= 1, start = sin(init[i]))

    for i in ids
        GRBaddgenconstrCos(backend(m), C_NULL, col(x[i]), col(c[i]), C_NULL)
        GRBaddgenconstrSin(backend(m), C_NULL, col(x[i]), col(s[i]), C_NULL)
    end

    @constraint m lift_matrix(d, r, α, M, c, s) .== 0
    @objective m Min sum(w .* lin_abs_angdiff_proxy.(c, s, θ))

    optimize!(m)
    _extract_solution(x, m)
end