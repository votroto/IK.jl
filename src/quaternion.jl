using StaticArrays
import LinearAlgebra: dot, cross, norm

struct Quaternion
    q0::Any
    q_::SVector{3,Any}
end

Quaternion(O::Number, v::Number) = Quaternion(O, SA[v, v, v])
Quaternion(O::Number) = Quaternion(O, 0)

Base.zero(::Type{Quaternion}) = Quaternion(0, SA[0, 0, 0])
Base.one(::Type{Quaternion}) = Quaternion(1, SA[0, 0, 0])
Base.:(+)(a::Quaternion, b::Quaternion) = Quaternion(a.q0 + b.q0, a.q_ + b.q_)
Base.:(-)(a::Quaternion, b::Quaternion) = Quaternion(a.q0 - b.q0, a.q_ - b.q_)
Base.:(*)(r::Number, a::Quaternion) = Quaternion(r * a.q0, r * a.q_)
Base.:(*)(a::Quaternion, b::Quaternion) = Quaternion(a.q0 * b.q0 - dot(a.q_, b.q_), a.q0 * b.q_ + b.q0 * a.q_ + cross(a.q_, b.q_))
Base.:(/)(a::Quaternion, r::Number) = 1 / r * a
Base.adjoint(a::Quaternion) = Quaternion(a.q0, -a.q_)
dot(a::Quaternion, b::Quaternion) = (a'b + b'a) / 2
cross(a::Quaternion, b::Quaternion) = (a * b - b'a') / 2
norm(a::Quaternion) = sqrt(dot(a, a))
part_scalar(a::Quaternion) = Quaternion(a.q0, SA[0, 0, 0])
part_vector(a::Quaternion) = Quaternion(0, a.q_)
Base.vec(a::Quaternion) = [a.q0, a.q_[1], a.q_[2], a.q_[3]]

struct DualQuaternion
    r::Quaternion
    d::Quaternion
end

Base.zero(::Type{DualQuaternion}) = DualQuaternion(zero(Quaternion), zero(Quaternion))
Base.:(+)(a::DualQuaternion, b::DualQuaternion) = DualQuaternion(a.r + b.r, a.d + b.d)
Base.:(-)(a::DualQuaternion, b::DualQuaternion) = DualQuaternion(a.r - b.r, a.d - b.d)
Base.:(*)(n::Number, a::DualQuaternion) = DualQuaternion(n * a.r, n * a.d)
Base.:(*)(a::DualQuaternion, b::DualQuaternion) = DualQuaternion(a.r * b.r, a.d*b.r + a.r*b.d)
Base.adjoint(a::DualQuaternion) = DualQuaternion(a.r', a.d')
dot(a::DualQuaternion, b::DualQuaternion) = (a'b + b'a) / 2
dot(a::DualQuaternion, b::DualQuaternion) = (a * b - b'a') / 2
circ(a::DualQuaternion, b::DualQuaternion) = DualQuaternion(dot(a.r, b.r) + dot(a.d, b.d), zero(Quaternion))
swap(a::DualQuaternion) = DualQuaternion(a.d, a.r)
norm(a::DualQuaternion) = sqrt(circ(a, a))
part_scalar(a::DualQuaternion) = DualQuaternion(part_scalar(a.r), part_scalar(a.d))
part_vector(a::DualQuaternion) = DualQuaternion(part_vector(a.r), part_vector(q.d))
Base.vec(a::DualQuaternion) = [vec(a.r); vec(a.d)]
norm(a::DualQuaternion) = sqrt(a.r * a.d')