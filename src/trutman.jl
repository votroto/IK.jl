function _rotation_to_quaternion(r::Matrix)
    q = [tr(r) + 1, r[3, 2] - r[2, 3], r[1, 3] - r[3, 1], r[2, 1] - r[1, 2]] / (2 * sqrt(tr(r) + 1))
    normalize!(q)
    Quaternion(q[1], q[2:4])
end

function _quaternion_to_rotation(u::Quaternion)
    R = [
        u.q0^2 + u.q_[1]^2 - u.q_[2]^2 - u.q_[3]^2  2 * (u.q_[1] * u.q_[2] - u.q0 * u.q_[3])    2 * (u.q_[1] * u.q_[3] + u.q0 * u.q_[2])
        2 * (u.q_[1] * u.q_[2] + u.q0 * u.q_[3])    u.q0^2 - u.q_[1]^2 + u.q_[2]^2 - u.q_[3]^2  2 * (u.q_[2] * u.q_[3] - u.q0 * u.q_[1])
        2 * (u.q_[1] * u.q_[3] - u.q0 * u.q_[2])    2 * (u.q_[2] * u.q_[3] + u.q0 * u.q_[1])    u.q0^2 - u.q_[1]^2 - u.q_[2]^2 + u.q_[3]^2  
    ]
    R / (u.q0^2 + u.q_[1]^2 + u.q_[2]^2 + u.q_[3]^2)
end

function rationalize_quaternion(q::Quaternion; tol=1e-3)
    q0  = rationalize(BigInt, q.q0; tol)
    q_1 = rationalize(BigInt, q.q_[1]; tol)
    q_2 = rationalize(BigInt, q.q_[2]; tol)
    q_3 = rationalize(BigInt, q.q_[3]; tol)
    Quaternion(q0, SA[q_1, q_2, q_3])
end

function rationalize_transformation(M::Matrix; tol=1e-3)
    Rflt = M[1:3, 1:3]
    Tflt = M[1:3, 4]

    T = rationalize.(BigInt, Tflt; tol)
    q = _rotation_to_quaternion(Rflt)
    R = _quaternion_to_rotation(rationalize_quaternion(q; tol=sqrt(tol)))
    z = zeros(eltype(T), 1, 3)

    [
        R T
        z 1
    ]
end
