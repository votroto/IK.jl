function _dh_RT(cα, sα, cx, sx, r, d)
    R = [
        cx -cα*sx sα*sx
        sx cα*cx -sα*cx
        0 sα cα
    ]
    T = [r * cx, r * sx, d]

    R, T
end

function _cos_sin_rat(α; tol=1e-3, T=BigInt)
    t = rationalize(T, tan(α / 2); tol)
    (1 - t^2) // (1 + t^2), (2t) // (1 + t^2)
end

function _dh_RT_rat(c, s, d, α, r; tol=1e-3, T=BigInt)
    cα, sα = _cos_sin_rat(α)
    _dh_RT(cα, sα, c, s, rationalize(T, r), rationalize(T, d))
end

_dh_matrix(R, T) = [R T; 0 0 0 1]
_dh_matrix_inverse(R, T) = [R' -R'T; 0 0 0 1]

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

function dh_matrix_rat(x, d, α, r; tol=1e-3)
    dh_matrix_rat(_cos_sin_rat(x; tol)..., d, α, r; tol)
end

function dh_matrix_rat_inverse(x, d, α, r; tol=1e-3)
    dh_matrix_rat_inverse(_cos_sin_rat(x; tol)..., d, α, r; tol)
end

function _dh_quaternion(cx::T, sx::T, d::Pd, cα::Pa, sα::Pa, r::Pr) where {T,Pd,Pa,Pr}
    Q = promote_type(T, Pd, Pa, Pr, Int)

    S = Quaternion{Q}(0, 0, 0, d)
    Z = Quaternion{Q}(cx, 0, 0, sx)
    A = Quaternion{Q}(0, r, 0, 0)
    X = Quaternion{Q}(cα, sα, 0, 0)

    ZZ = DualQuaternion(Z, S * Z / 2)
    XX = DualQuaternion(X, A * X / 2)

    ZZ * XX
end

function dh_quaternion(c, s, d, α, r)
    cα, sα = cos(α / 2), sin(α / 2)
    _dh_quaternion(c, s, d, cα, sα, r)
end

dh_quaternion(x, d, α, r) = dh_quaternion(cos(x / 2), sin(x / 2), d, α, r)
dh_quaternion_inverse(c, s, d, α, r) = dh_quaternion(c, s, d, α, r)'
dh_quaternion_inverse(x, d, α, r) = dh_quaternion_inverse(cos(x / 2), sin(x / 2), d, α, r)

function dh_quaternion_rat(c, s, d, α, r; tol=1e-3, T=BigInt)
    cα, sα = _cos_sin_rat(α / 2; tol=1e-3)
    _dh_quaternion(c, s, rationalize(T, d), cα, sα, rationalize(T, r))
end

function dh_quaternion_rat(x, d, α, r; tol=1e-3)
    @show cx, sx = _cos_sin_rat(x / 2; tol)
    @show dh_quaternion_rat(cx, sx, d, α, r; tol)
end

function dh_quaternion_rat_inverse(c, s, d, α, r; tol=1e-3)
    dh_quaternion_rat(c, s, d, α, r; tol)'
end

function dh_quaternion_rat_inverse(x, d, α, r; tol=1e-3)
    cx, sx = _cos_sin_rat(x / 2; tol)
    dh_quaternion_rat_inverse(cx, sx, d, α, r; tol)
end
