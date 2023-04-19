using JuMP

import Base.:*
Base.:*(a::GenericQuadExpr, b::GenericQuadExpr) = jump_quadratic_product(a, b)

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
        model = a.model
        lb, ub = jump_quadratic_bounds(a, b)
        start = jump_quadratic_start(a, b)
        name = string(a) * string(b)

        v = @variable(model, lower_bound = lb, upper_bound = ub, start = start, base_name = name)
        @constraint(model, v == a * b)

        v
end

function jump_quadratic_envelope(x, y)
        model = x.model

 	lx = has_lower_bound(x) ? lower_bound(x) : -1
        ux = has_upper_bound(x) ? upper_bound(x) : 1
        ly = has_lower_bound(y) ? lower_bound(y) : -1
        uy = has_upper_bound(y) ? upper_bound(y) : 1
	l, h = jump_quadratic_bounds(x, y)

        z = @variable(model, lower_bound = l, upper_bound = h)

        @constraint(model, -z + lx*y + ly*x <= lx*ly)
	@constraint(model, -z + ux*y + uy*x <= ux*uy)
	@constraint(model, -z + ux*y + ly*x >= ux*ly)
	@constraint(model, -z + lx*y + uy*x >= lx*uy)

        z
end

function jump_quadratic_envelope(a::QuadExpr)
        a.aff + sum(c * jump_quadratic_envelope(x.a, x.b) for (x, c) in a.terms; init=0)
end



function jump_quadratic_product(a::QuadExpr, b::AffExpr)
        na = a.aff + sum(c * jump_quadratic_lift(x.a, x.b) for (x, c) in a.terms; init=0)
        na * b
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
        if a == 1
                return b
        elseif isempty(a.terms) && isempty(b.terms)
                return a.aff * b.aff
        end

        na = a.aff + sum(c * jump_quadratic_lift(x.a, x.b) for (x, c) in a.terms; init=0)
        nb = b.aff + sum(c * jump_quadratic_lift(x.a, x.b) for (x, c) in b.terms; init=0)
        na * nb
end

function jump_quadratic_product(a::AbstractMatrix, b::AbstractMatrix)
        [mapreduce(jump_quadratic_product, +, a[i, :], b[:, j]) for i in 1:4, j in 1:4]
end

jump_quadratic_product(a::AffExpr, b::QuadExpr) = jump_quadratic_product(b, a)
jump_quadratic_product(a::QuadExpr, b::VariableRef) = jump_quadratic_product(b, a)
jump_quadratic_product(a, b) = a * b