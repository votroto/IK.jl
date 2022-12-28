using JuMP

import Base.:*
Base.:*(a::GenericQuadExpr,b::GenericQuadExpr) = jump_quadratic_product(a, b)

function jump_quadratic_bounds(a, b)
        la = has_lower_bound(a) ? lower_bound(a) : -Inf
        ua = has_upper_bound(a) ? upper_bound(a) : Inf
        lb = has_lower_bound(b) ? lower_bound(b) : -Inf
        ub = has_upper_bound(b) ? upper_bound(b) : Inf
        lo = min(la * lb, la * ub, ua * lb, ua * ub)
        hi = max(la * lb, la * ub, ua * lb, ua * ub)
        lo, hi
end

function jump_quadratic_start(a, b)
	sa = start_value(a)
	sb = start_value(b)

	if isnothing(sa) || isnothing(sa)
		nothing
	else
		sa * sb
	end
end

function jump_quadratic_lift(a, b)
        m = a.model
        l, u = jump_quadratic_bounds(a, b)
        st = jump_quadratic_start(a, b)

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

function jump_quadratic_product(a::QuadExpr, b::QuadExpr)
        na = a.aff + sum(c * jump_quadratic_lift(x.a, x.b) for (x, c) in a.terms; init=0)
        nb = b.aff + sum(c * jump_quadratic_lift(x.a, x.b) for (x, c) in b.terms; init=0)
        na * nb
end

function jump_quadratic_product(a::AbstractMatrix, b::AbstractMatrix)
        [mapreduce(jump_quadratic_product, +, a[i, :], b[:, j]) for i in 1:4, j in 1:4]
end

jump_quadratic_product(a::AffExpr, b::QuadExpr) = jump_quadratic_product(b, a)

jump_quadratic_product(a, b) = a * b