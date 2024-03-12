function _dh_RTz(cα, sα, cx, sx, r, d)
    R = [
        cx -cα*sx sα*sx
        sx cα*cx -sα*cx
        0 sα cα
    ]
    T = [r * cx, r * sx, d]
    z = zeros(eltype(T), 1, 3)

    R, T, z
end

function dh_matrix(c, s, d, α, r)
    R, T, z = _dh_RTz(cos(α), sin(α), c, s, r, d)

    DH = [
        R T
        z 1
    ]
end

function dh_matrix_inverse(c, s, d, α, r)
    R, T, z = _dh_RTz(cos(α), sin(α), c, s, r, d)

    DH = [
        R' -R'T
        z 1
    ]
end

dh_matrix(x, d, α, r) = dh_matrix(cos(x), sin(x), d, α, r)
dh_matrix_inverse(x, d, α, r) = dh_matrix_inverse(cos(x), sin(x), d, α, r)

function dh_matrix_rat(c, s, d, α, r; tol=1e-3)
    t = rationalize(tan(α/2); tol)
    cα, sα = (1 - t^2) // (1 + t^2), (2t) // (1 + t^2)

    R, T, z = _dh_RTz(cα, sα, c, s, rationalize(r), rationalize(d))
    DH = [
        R T
        z 1
    ]
end

function dh_matrix_rat_inverse(c, s, d, α, r; tol=1e-3)
    t = rationalize(tan(α/2); tol)
    cα, sα = (1 - t^2) // (1 + t^2), (2t) // (1 + t^2)
    
    R, T, z = _dh_RTz(cα, sα, c, s, rationalize(r), rationalize(d))
    DH = [
        R' -R'T
        z 1
    ]
end


function dh_matrix_rat_rat(x, d, α, r; tol=1e-3)
    t = rationalize(BigInt, tan(α/2); tol)
    cα, sα = (1 - t^2) // (1 + t^2), (2t) // (1 + t^2)
    tx = rationalize(BigInt, tan(x/2); tol)
    cx, sx = (1 - tx^2) // (1 + tx^2), (2tx) // (1 + tx^2)

    R, T, z = _dh_RTz(cα, sα, cx, sx, rationalize(r), rationalize(d))
    DH = [
        R T
        z 1
    ]
end

dh_matrix_rat(x, d, α, r) = dh_matrix_rat(cos(x), sin(x), d, α, r)
dh_matrix_rat_inverse(x, d, α, r) = dh_matrix_rat_inverse(cos(x), sin(x), d, α, r)

function dh_quaternion(c, s, d, α, r)
    S = Quaternion(0, SA[0, 0, d])
    Z = Quaternion(c, SA[0, 0, s])
    A = Quaternion(0, SA[r, 0, 0])
    X = Quaternion(cos(α / 2), SA[sin(α / 2), 0, 0])

    ZZ = DualQuaternion(Z, S * Z / 2)
    XX = DualQuaternion(X, A * X / 2)

    ZZ * XX
end

dh_quaternion(x, d, α, r) = dh_quaternion(cos(x / 2), sin(x / 2), d, α, r)
dh_quaternion_inverse(c, s, d, α, r) = dh_quaternion(c, s, d, α, r)'
dh_quaternion_inverse(x, d, α, r) = dh_quaternion_inverse(cos(x / 2), sin(x / 2), d, α, r)