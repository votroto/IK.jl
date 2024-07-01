using JuMP
using Gurobi
using SCIP


function _scip()
    optimizer_with_attributes(
        SCIP.Optimizer,
        "lp/threads" => 4
    )
end


function _default_optimizer()
    optimizer_with_attributes(
        Gurobi.Optimizer,
        #MOI.Silent() => true,
        "FuncNonlinear" => 1,
        "Nonconvex" => 2,
        "Presolve" => 2,
        "Threads" => 4
    )
end


function _quii()
    optimizer_with_attributes(
        Gurobi.Optimizer,
        MOI.Silent() => true,
        "FuncNonlinear" => 1,
        "Nonconvex" => 2,
        "Presolve" => 2,
        "Threads" => 1
    )
end

function _feasible_optimizer()
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
function zsolve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
    lift_method=lift_matrix, optimizer=_default_optimizer(), init=θ)

    ids = eachindex(d)
    m = direct_model(optimizer)

    @variable(m, W[1:3,1:4])
    @variable(m, Z >= 0 )

    @variable(m, a[ids])

    @variable(m, c[ids])
    @variable(m, s[ids])
    constrain_trig_vars.(c, s, θl, θh, init)

    @constraint m lift_method(d, r, α, M, c, s) .== W
    @constraint m c .^ 2 .+ s .^ 2 .<= 1
    @constraint m lin_angdiff_proxy.(c, s, θl) .>= 0
    @constraint m lin_angdiff_proxy.(c, s, θh) .<= 0


    @constraint m W .<= Z
    @constraint m W .>= -Z
    #GRBaddgenconstrNorm(backend(m), C_NULL, column(Z), 12, [column(W[i]) for i in eachindex(W)], GRB_INFINITY)


    @objective m Min Z
    optimize!(m)

    _extract_solution(c, s, m)
end



function bsolve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
    lift_method=lift_matrix, optimizer=_default_optimizer(), init=θ)

    ids = eachindex(d)
    m = direct_model(optimizer)

    @variable(m, W[1:3,1:4])
    @variable(m, Z >= 0)



    xstrt = init

    cstrt = cos.(init)
    sstrt = sin.(init)
    @variable(m, θl[i] <= x[i in ids] <= θh[i], start = xstrt[i])

    @variable(m, -1 <= c[i in ids] <= 1, start = cstrt[i])
    @variable(m, -1 <= s[i in ids] <= 1, start = sstrt[i])

    constrain_trig_vars.(c, s, θl, θh, init)

    @constraint m lift_method(d, r, α, M, c, s) .== W

    for i in ids
        GRBaddgenconstrCos(backend(m), C_NULL, column(x[i]), column(c[i]), C_NULL)
        GRBaddgenconstrSin(backend(m), C_NULL, column(x[i]), column(s[i]), C_NULL)
    GRBsetintattrelement(backend(m), GRB_INT_ATTR_BRANCHPRIORITY, column(x[i]), 0)
    GRBsetintattrelement(backend(m), GRB_INT_ATTR_BRANCHPRIORITY, column(c[i]), 1)
    GRBsetintattrelement(backend(m), GRB_INT_ATTR_BRANCHPRIORITY, column(s[i]), 2)
    end


    GRBaddgenconstrNorm(backend(m), C_NULL, column(Z), 12, [column(W[i]) for i in eachindex(W)], 1.0)



    @objective m Min Z
    optimize!(m)

    _extract_solution(c, s, m)
end


column(x::VariableRef) = Gurobi.c_column(backend(owner_model(x)), index(x))




#=
function ff(d, r, α, θl, θh, M, θ, w;
    lift_method=lift_matrix, optimizer=_default_optimizer(), init=θ)



    ids = eachindex(d)
    @polyvar C[ids] S[ids]

    E = build_pose_constraint_poly(d, r, α, C, S, M)

    println("minimize Objective: ", sum(w .* lin_abs_angdiff_proxy.(C, S, θ)), ";")
    for i in 1:16
        println("s.t. M$i: ", E[i], " = 0;")
    end

    for i in ids
        println("s.t. C$i: ", C[i] .^ 2 .+ S[i] .^ 2 , " = 1;")
        println("s.t. L$i: ", lin_angdiff_proxy(C[i], S[i], θl[i]) , " >= 0;")
        println("s.t. U$i: ", lin_angdiff_proxy(C[i], S[i], θh[i]) , " <= 0;")
    end

end
=#


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


function fsolve_inverse_kinematics(d, r, α, θl, θh, M, θ, w;
    lift_method=lift_matrix, optimizer=_default_optimizer(), init=θ)

    ids = eachindex(d)
    #m = Model(optimizer)
    m = direct_model(optimizer)

     @variable(m, c[ids])
     @variable(m, s[ids])
    constrain_trig_vars.(c, s, θl, θh, init)

@constraint(m, c[7:end] .== cos.(init[7:end]))
@constraint(m, s[7:end] .== sin.(init[7:end]))

    @constraint m lift_method(d, r, α, M, c, s) .== 0
    @constraint m c .^ 2 .+ s .^ 2 .== 1
    @constraint m lin_angdiff_proxy.(c, s, θl) .>= 0
    @constraint m lin_angdiff_proxy.(c, s, θh) .<= 0

    @objective m Min sum(w .* lin_abs_angdiff_proxy.(c, s, θ))

    optimize!(m)

    _extract_solution(c, s, m)
end










function ssii(d, r, α, θl, θh, M, θ, w;
    lift_method=lift_matrix, init=θ, optimizer = _default_optimizer())

    ids = eachindex(θ)
    m = direct_model(optimizer)

    xstrt = init
    cstrt = cos.(init)
    sstrt = sin.(init)
    @variable(m, θl[i] <= x[i in ids] <= θh[i], start = xstrt[i])
    @variable(m, -1 <= c[i in ids] <= 1, start = cstrt[i])
    @variable(m, -1 <= s[i in ids] <= 1, start = sstrt[i])

    constrain_trig_vars.(c, s, θl, θh, init)

    for i in ids
        GRBaddgenconstrCos(backend(m), C_NULL, column(x[i]), column(c[i]), C_NULL)
        GRBaddgenconstrSin(backend(m), C_NULL, column(x[i]), column(s[i]), C_NULL)
    end


    @constraint m lift_method(d, r, α, M, c, s) .== 0
    @objective m Min sum(w .* lin_abs_angdiff_proxy.(c, s, θ))







    function my_callback_function(cb_data, cb_where::Cint)
        # You can reference variables outside the function as normal
        # You can select where the callback is run
        if  cb_where != GRB_CB_MIPSOL
            return
        end

        println("staart")
        println(cb_where)
        # Before querying `callback_value`, you must call:
        Gurobi.load_callback_variable_primal(cb_data, cb_where)
        @show x_val = callback_value.(Ref(cb_data), x)

        nx, obj1, ret1, tim1 = fsolve_inverse_kinematics(d, r, α, θl, θh, M, θ, w; init=x_val, optimizer=_quii())

        println(sum(w .* lin_abs_angdiff_proxy.(cos.(nx), sin.(nx), θ)))
        println(MOI.submit(m, MOI.HeuristicSolution(cb_data), [x;c;s], [nx;cos.(nx);sin.(nx)]))

        println("eend")
    end
    MOI.set(m, Gurobi.CallbackFunction(), my_callback_function)





    optimize!(m)


    _extract_solution(c, s, m)
end


