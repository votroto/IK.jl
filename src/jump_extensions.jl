using JuMP

import Base.:*
Base.:*(a::GenericQuadExpr, b::GenericQuadExpr) = jump_quadratic_product(a, b)
Base.:*(a::GenericAffExpr, b::GenericQuadExpr) = jump_quadratic_product(a, b)
Base.:*(a::GenericQuadExpr, b::GenericAffExpr) = jump_quadratic_product(a, b)

function _start_value(x)
    v = start_value(x)
    isnothing(v) ? NaN : v
end

_q_prod_start(a, b) = _start_value(a) * _start_value(b)

function _q_prod_bounds(a, b)
    la = has_lower_bound(a) ? lower_bound(a) : -Inf
    ua = has_upper_bound(a) ? upper_bound(a) : Inf
    lb = has_lower_bound(b) ? lower_bound(b) : -Inf
    ub = has_upper_bound(b) ? upper_bound(b) : Inf
    lo = min(la * lb, la * ub, ua * lb, ua * ub)
    hi = max(la * lb, la * ub, ua * lb, ua * ub)
    lo, hi
end

MM = Nothing
LBS = []
UBS = []
VS = []
VAAA = Dict()

function _q_mono_lift(a, b)
    global MM, LBS, UBS, VS, VAAA

    model = a.model

if MM != model
    MM = model
    LBS = []
    UBS = []
    VS = []
    VAAA = Dict()

end



    lb, ub = _q_prod_bounds(a, b)
    start = _q_prod_start(a, b)
    name = "($a$b)"

if (a,b) in keys(VAAA)
    @show "yes"
    return VAAA[(a,b)]
end


LBS = [LBS; lb]
UBS = [UBS; lb]
VS = [VS; (a,b)]


    var = @variable(model,
        lower_bound = lb,
        upper_bound = ub,
        start = start,
        base_name = name
    )

    VAAA[(a,b)] = var

    @constraint(model, var == a * b)

    var
end

function _q_expr_lift(a::QuadExpr; init=0)
    a.aff + sum(c * _q_mono_lift(x.a, x.b) for (x, c) in a.terms; init)
end

_simplify_type(a::QuadExpr) = isempty(a.terms) ? _simplify_type(a.aff) : a
_simplify_type(a::AffExpr) = isempty(a.terms) ? a.constant : a
_simplify_type(a) = a

_q_prod(a::QuadExpr, b::QuadExpr) = _q_expr_lift(a) * _q_expr_lift(b)
_q_prod(a::QuadExpr, b::AffExpr) = _q_expr_lift(a) * b
_q_prod(a::VariableRef, b::QuadExpr) = a * _q_expr_lift(b)
_q_prod(a::AffExpr, b::QuadExpr) = jump_quadratic_product(b, a)
_q_prod(a::QuadExpr, b::VariableRef) = jump_quadratic_product(b, a)
_q_prod(a, b) = a * b

function jump_quadratic_product(a, b)
    _q_prod(_simplify_type(a), _simplify_type(b))
end

function jump_quadratic_product(a::AbstractMatrix, b::AbstractMatrix)
    is, js = axes(a)
    [mapreduce(*, +, a[i, :], b[:, j]) for i in is, j in js]
end