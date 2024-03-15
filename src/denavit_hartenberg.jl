function _dh_RT(cα, sα, cx, sx, r, d)
    R = [
        cx -cα*sx sα*sx
        sx cα*cx -sα*cx
        0 sα cα
    ]
    T = [r * cx, r * sx, d]

    R, T
end

function _dh_RT_rat(c, s, d, α, r; tol=1e-3, T=BigInt)
    t = rationalize(T, tan(α / 2); tol)
    cα, sα = (1 - t^2) // (1 + t^2), (2t) // (1 + t^2)
    _dh_RT(cα, sα, c, s, rationalize(T, r), rationalize(T, d))
end

_dh_matrix(R, T) = [R T; zeros(eltype(T), 1, 3) 1]
_dh_matrix_inverse(R, T) = [R' -R'T; zeros(eltype(T), 1, 3) 1]

function dh_matrix(c, s, d, α, r)
    R, T = _dh_RT(cos(α), sin(α), c, s, r, d)
    _dh_matrix(R, T)
end

function dh_matrix_inverse(c, s, d, α, r)
    R, T = _dh_RT(cos(α), sin(α), c, s, r, d)
    _dh_matrix_inverse(R, T)
end

dh_matrix(x, d, α, r) = dh_matrix(cos(x), sin(x), d, α, r)
dh_matrix_inverse(x, d, α, r) = dh_matrix_inverse(cos(x), sin(x), d, α, r)

function dh_matrix_rat(c, s, d, α, r; tol=1e-3)
    R, T = _dh_RT_rat(c, s, d, α, r; tol)
    _dh_matrix(R, T)
end

function dh_matrix_rat_inverse(c, s, d, α, r; tol=1e-3)
    R, T = _dh_RT_rat(c, s, d, α, r; tol)
    _dh_matrix_inverse(R, T)
end

dh_matrix_rat(x, d, α, r) = dh_matrix_rat(cos(x), sin(x), d, α, r)
dh_matrix_rat_inverse(x, d, α, r) = dh_matrix_rat_inverse(cos(x), sin(x), d, α, r)

function dh_quaternion(c::T, s::T, d::Pd, α::Pa, r::Pr) where {T, Pd, Pa, Pr}
    Q = promote_type(T, Pd, Pa, Pr, Int)

    S = Quaternion{Q}(0, 0, 0, d)
    Z = Quaternion{Q}(c, 0, 0, s)
    A = Quaternion{Q}(0, r, 0, 0)
    X = Quaternion{Q}(cos(α / 2), sin(α / 2), 0, 0)

    ZZ = DualQuaternion(Z, S * Z / 2)
    XX = DualQuaternion(X, A * X / 2)

    ZZ * XX
end

dh_quaternion(x, d, α, r) = dh_quaternion(cos(x / 2), sin(x / 2), d, α, r)
dh_quaternion_inverse(c, s, d, α, r) = dh_quaternion(c, s, d, α, r)'
dh_quaternion_inverse(x, d, α, r) = dh_quaternion_inverse(cos(x / 2), sin(x / 2), d, α, r)