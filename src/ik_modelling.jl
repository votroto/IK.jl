
function lin_angdiff_approx(cosx, sinx, y)
    2 * (1 - cosx * cos(y) - sinx * sin(y))
end

function constrain_trig_vars(c, s, θl, θh; init)
    sl, sh = sin_min_max(θl, θh)
    cl, ch = cos_min_max(θl, θh)

    set_lower_bound(c, cl)
    set_lower_bound(s, sl)

    set_upper_bound(c, ch)
    set_upper_bound(s, sh)

    set_start_value(c, cos(init))
    set_start_value(s, sin(init))
end

function _split_manipulator(ids)
    mid = div(length(ids), 2)
    f, s = take(ids, mid), drop(ids, mid)
    f, reverse(collect(s))
end

function extract_solution(c, s, m)
    stat = termination_status(m)

    sol = has_values(m) ? atan.(value.(s), value.(c)) : missing
    obj = stat == OPTIMAL ? objective_value(m) : missing

    sol, obj, stat, solve_time(m)
end

function build_eqs(d, r, α, M)
    T(i) = dh_lin_t(c[i], s[i], d[i], α[i], r[i])
    iT(i) = dh_lin_inv_t(c[i], s[i], d[i], α[i], r[i])

    ids = eachindex(d)
    fwd, rev = _split_manipulator(ids)

    @polyvar c[ids] s[ids]

    prod(T, fwd) .- M * prod(iT, rev), c, s
end