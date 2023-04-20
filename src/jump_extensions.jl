using JuMP

import Base.:*
Base.:*(a::GenericQuadExpr, b::GenericQuadExpr) = jump_quadratic_product(a, b)

function _start_value(x)
        v = start_value(x)
        isnothing(v) ? NaN : v
end

jump_quadratic_start(a, b) = _start_value(a) * _start_value(b)

function jump_quadratic_bounds(a, b)
        la = has_lower_bound(a) ? lower_bound(a) : -Inf
        ua = has_upper_bound(a) ? upper_bound(a) : Inf
        lb = has_lower_bound(b) ? lower_bound(b) : -Inf
        ub = has_upper_bound(b) ? upper_bound(b) : Inf
        lo = min(la * lb, la * ub, ua * lb, ua * ub)
        hi = max(la * lb, la * ub, ua * lb, ua * ub)
        lo, hi
end

function jump_quadmono_lift(a, b)
        model = a.model
        lower_bound, upper_bound = jump_quadratic_bounds(a, b)
        start = jump_quadratic_start(a, b)
        base_name = string(a) * string(b)

        v = @variable(model, lower_bound, upper_bound, start, base_name)
        @constraint(model, v == a * b)

        v
end

function jump_quadexpr_lift(a::QuadExpr; init=0)
        a.aff + sum(c * jump_quadmono_lift(x.a, x.b) for (x, c) in a.terms; init)
end

function jump_quadratic_product(a::QuadExpr, b::QuadExpr)
        if a == 1
                return b
        elseif isempty(a.terms) && isempty(b.terms)
                return a.aff * b.aff
        end

        jump_quadexpr_lift(a) * jump_quadexpr_lift(b)
end

function jump_quadratic_product(a::AbstractMatrix, b::AbstractMatrix)
        is, js = axes(a)
        [mapreduce(jump_quadratic_product, +, a[i, :], b[:, j]) for i in is, j in js]
end

jump_quadratic_product(a::QuadExpr, b::AffExpr) = jump_quadexpr_lift(a) * b
jump_quadratic_product(a::VariableRef, b::QuadExpr) = a * jump_quadexpr_lift(b)

jump_quadratic_product(a::AffExpr, b::QuadExpr) = jump_quadratic_product(b, a)
jump_quadratic_product(a::QuadExpr, b::VariableRef) = jump_quadratic_product(b, a)
jump_quadratic_product(a, b) = a * b