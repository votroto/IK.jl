using StaticArrays
import LinearAlgebra: dot, cross, norm

struct Quaternion{T}
    q0::T
    q1::T
    q2::T
    q3::T
end

Quaternion(O::T, v::SVector{3, T}) where T = Quaternion{T}(O, v[1], v[2], v[3])
Quaternion(O::Number, v::Number) = Quaternion(O, v, v, v)
Quaternion(O::Number) = Quaternion(O, 0)

_dot(a, b) = sum(a[i] * b[i] for i in eachindex(a))

Base.zero(::Type{Quaternion}) = Quaternion(0, 0, 0, 0)
Base.one(::Type{Quaternion}) = Quaternion(1, 0, 0, 0)
Base.:(+)(a::Quaternion, b::Quaternion) = Quaternion(a.q0 + b.q0, SA[a.q1, a.q2, a.q3] + SA[b.q1, b.q2, b.q3])
Base.:(-)(a::Quaternion, b::Quaternion) = Quaternion(a.q0 - b.q0, SA[a.q1, a.q2, a.q3] - SA[b.q1, b.q2, b.q3])
Base.:(*)(r::Number, a::Quaternion) = Quaternion(r * a.q0, r * SA[a.q1, a.q2, a.q3])
function Base.:(*)(a::Quaternion, b::Quaternion)
    q0 = a.q0 * b.q0 - _dot(SA[a.q1, a.q2, a.q3], SA[b.q1, b.q2, b.q3])
    q_ = a.q0 * SA[b.q1, b.q2, b.q3] + b.q0 * SA[a.q1, a.q2, a.q3] + cross(SA[a.q1, a.q2, a.q3], SA[b.q1, b.q2, b.q3])
    Quaternion{promote_type(typeof(q0), eltype(q_))}(q0, q_...)
end
Base.:(/)(a::Quaternion, r::Number) = 1 / r * a
Base.:(/)(a::Quaternion, r::Integer) = 1 // r * a
Base.adjoint(a::Quaternion) = Quaternion(a.q0, -a.q1, -a.q2, -a.q3)
dot(a::Quaternion, b::Quaternion) = ((a'b + b'a) / 2).q0
cross(a::Quaternion, b::Quaternion) = (a * b - b'a') / 2
norm(a::Quaternion) = sqrt(_dot(a, a))
Base.vec(a::Quaternion) = [a.q0, a.q1, a.q2, a.q3]

struct DualQuaternion
    r::Quaternion
    d::Quaternion
end

Base.zero(::Type{DualQuaternion}) = DualQuaternion(zero(Quaternion), zero(Quaternion))
Base.:(+)(a::DualQuaternion, b::DualQuaternion) = DualQuaternion(a.r + b.r, a.d + b.d)
Base.:(-)(a::DualQuaternion, b::DualQuaternion) = DualQuaternion(a.r - b.r, a.d - b.d)
Base.:(*)(n::Number, a::DualQuaternion) = DualQuaternion(n * a.r, n * a.d)
Base.:(*)(a::DualQuaternion, b::DualQuaternion) = DualQuaternion(a.r * b.r, a.d * b.r + a.r * b.d)
Base.adjoint(a::DualQuaternion) = DualQuaternion(a.r', a.d')
dot(a::DualQuaternion, b::DualQuaternion) = (a'b + b'a) / 2
cross(a::DualQuaternion, b::DualQuaternion) = (a * b - b'a') / 2
circ(a::DualQuaternion, b::DualQuaternion) = dot(a.r, b.r) + dot(a.d, b.d)
swap(a::DualQuaternion) = DualQuaternion(a.d, a.r)
norm(a::DualQuaternion) = sqrt(circ(a, a))
Base.vec(a::DualQuaternion) = [vec(a.r); vec(a.d)]