function _rotation_to_quaternion(r::Matrix)
    f = (2 * sqrt(tr(r) + 1))
    q = [tr(r) + 1, r[3, 2] - r[2, 3], r[1, 3] - r[3, 1], r[2, 1] - r[1, 2]] / f
    normalize!(q)
    Quaternion(q[1], q[2:4]...)
end

function _quaternion_to_rotation(u::Quaternion)
    q0, q1, q2, q3 = u.q0, u.q1, u.q2, u.q3
    R = [
        q0^2+q1^2-q2^2-q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2)
        2*(q1*q2+q0*q3) q0^2-q1^2+q2^2-q3^2 2*(q2*q3-q0*q1)
        2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) q0^2-q1^2-q2^2+q3^2
    ]
    R / (q0^2 + q1^2 + q2^2 + q3^2)
end

function rationalize_quaternion(q::Quaternion; tol=1e-3)
    q0 = rationalize(BigInt, q.q0; tol)
    q_1 = rationalize(BigInt, q.q1; tol)
    q_2 = rationalize(BigInt, q.q2; tol)
    q_3 = rationalize(BigInt, q.q3; tol)
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
