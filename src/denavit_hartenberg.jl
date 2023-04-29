
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

dh_t(x, d, α, r) = dh_lin_t(cos(x), sin(x), d, α, r)
dh_inv_t(x, d, α, r) = dh_lin_inv_t(cos(x), sin(x), d, α, r)