using JuMP


function jump_variable_bounds(a, b)
        la = lower_bound(a)
        ua = upper_bound(a)
        lb = lower_bound(b)
        ub = upper_bound(b)
        lo = min(la * lb, la * ub, ua * lb, ua * ub)
        hi = max(la * lb, la * ub, ua * lb, ua * ub)
        lo, hi
end

function jump_quadratic_lift(a, b)
        m = a.model
        l, u = jump_variable_bounds(a, b)
        st = start_value(a) * start_value(b)

        v = @variable(m, lower_bound = l, upper_bound = u, start = st)
        @constraint(m, v == a * b)

        v
end

function jump_quadratic_product(a::QuadExpr, b::AffExpr)
        na = a.aff + sum(c * jump_quadratic_lift(x.a, x.b) for (x, c) in a.terms; init=0)
        na * b
end

function jump_quadratic_product(a::VariableRef, b::QuadExpr)
        na = b.aff + sum(c * jump_quadratic_lift(x.a, x.b) for (x, c) in b.terms; init=0)
        na * a
end

jump_quadratic_product(a::AffExpr, b::QuadExpr) = jump_quadratic_product(b, a)
function jump_quadratic_product(a::QuadExpr, b::QuadExpr)
        na = a.aff + sum(c * jump_quadratic_lift(x.a, x.b) for (x, c) in a.terms; init=0)
        nb = b.aff + sum(c * jump_quadratic_lift(x.a, x.b) for (x, c) in b.terms; init=0)
        na * nb
end


function jump_quadratic_product(a::AbstractMatrix, b::AbstractMatrix)
        [mapreduce(jump_quadratic_product, +, a[i, :], b[:, j]) for i in 1:4, j in 1:4]
end

jump_quadratic_product(a, b) = a * b