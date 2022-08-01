include("david_hartenberg.jl")

function solve_forward_kinematics(x, d, r, α)
    ids = eachindex(x)
    T(i) = dh_t(x[i], d[i], α[i], r[i])

    prod(T[i] for i in ids)
end