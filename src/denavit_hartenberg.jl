
function dh_lin_inv_t(c, s, d, α, r)
    z = [
        c s 0 0
        -s c 0 0
        0 0 1 -d
        0 0 0 1
    ]
    x = [
        1 0 0 -r
        0 cos(α) sin(α) 0
        0 -sin(α) cos(α) 0
        0 0 0 1
    ]
    x * z
end

function dh_lin_t(c, s, d, α, r)
    z = [
        c -s 0 0
        s c 0 0
        0 0 1 d
        0 0 0 1
    ]
    x = [
        1 0 0 r
        0 cos(α) -sin(α) 0
        0 sin(α) cos(α) 0
        0 0 0 1
    ]
    z * x
end

function dq_lin(c, s, d, α, r)
    S = Quaternion(0, SA[0, 0, d])
    Z = Quaternion(c, SA[0, 0, s])
    A = Quaternion(0, SA[r, 0, 0])
    X = Quaternion(cos(α / 2), SA[sin(α / 2), 0, 0])

    ZZ = DualQuaternion(Z, S * Z / 2)
    XX = DualQuaternion(X, A * X / 2)

    ZZ*XX
end

dq_lin_inv(c, s, d, α, r) = dq_lin(c, s, d, α, r)'
dq(x, d, α, r) = dq_lin(cos(x), sin(x), d, α, r)
dq_inv(x, d, α, r) = dq_lin_inv(cos(x), sin(x), d, α, r)
dh_t(x, d, α, r) = dh_lin_t(cos(x), sin(x), d, α, r)
dh_inv_t(x, d, α, r) = dh_lin_inv_t(cos(x), sin(x), d, α, r)