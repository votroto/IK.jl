using LinearAlgebra

# https://www.kinovarobotics.com/uploads/User-Guide-Gen3-R07.pdf

# The docs could not be more confusing...

function params_kinova_gen3()
    r = [0.0054, -0.2104, -0.0064, -0.2084, 0.0, -0.1059, 0]
    d = [-0.1284, -0.0064, -0.2104, -0.0064, -0.1059, 0.0, -0.0615]
    α = [π / 2, -π / 2, π / 2, -π / 2, π / 2, -π / 2, π]
    θl = [-π, -π, -π, -π, -π, -π, -π]
    θh = [π, π, π, π, π, π, π]

    w = normalize(ones(7), 1)
    θ = zeros(7)

    d, r, α, θl, θh, w, θ
end