function map_monomials(f, poly)
    sum(coefficients(poly) .* map(f, monomials(poly)), init=0)
end

function lift_poly(lift_vars, poly)
    map_monomials(e -> lift_vars[e], poly)
end

"""Linear proxy expression for minimizing abs(angdiff(x, y))"""
function lin_abs_angdiff_proxy(cosx, sinx, y)
    2 * (1 - cosx * cos(y) - sinx * sin(y))
end

"""Linear proxy for angdiff(x, y)"""
lin_angdiff_proxy(cosx, sinx, y) = -((cosx + 1) * tan(y / 2) - sinx)

function _set_vat_lb_ub_st(x, lb, ub, st)
    set_lower_bound(x, lb)
    set_upper_bound(x, ub)
    set_start_value(x, st)
end

"""Sets the box constraints and the initial guess for cos(x) and sin(x)"""
function constrain_trig_vars(c, s, θl, θh, init)
    _set_vat_lb_ub_st(c, cos_min_max(θl, θh)..., cos(init))
    _set_vat_lb_ub_st(s, sin_min_max(θl, θh)..., sin(init))
end

function _split_manipulator(ids)
    mid = div(length(ids), 2)
    f, s = take(ids, mid), drop(ids, mid)

    f, reverse(collect(s))
end

function build_pose_constraint(d, r, α, c, s)
    T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
    iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

    fwd, rev = _split_manipulator(eachindex(d))

    map(T, fwd), map(iT, rev)
	
end

function build_pose_constraint_poly(d, r, α, c, s, M)
	fwd, rev = build_pose_constraint(d, r, α, c, s)
	
    chain_poly_dirty = prod(fwd) - M * prod(rev)
    chain_poly_clean = mapcoefficients.(round_zero, chain_poly_dirty)

	chain_poly_clean
end