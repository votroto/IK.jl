using LinearAlgebra

# https://www.iri.upc.edu/people/thomas/papers/PROBLEMS.pdf

function params_canadarm2()
    r = [0, 0, 6850, 6850, 0, 0, 0] ./ 1000.0
    d = [380, 635, 504, 0, 504, 635, 380] ./ 1000.0
    α = [-π / 2, π / 2, 0, 0, -π / 2, π / 2, 0]
    θl = [-π, -π, -π, -π, -π, -π, -π]
    θh = [π, π, π, π, π, π, π]

    w = normalize(ones(7), 1)
    θ = zeros(7)

    d, r, α, θl, θh, w, θ
end