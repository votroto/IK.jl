using LinearAlgebra

# https://www.kinovarobotics.com/uploads/User-Guide-Gen3-R07.pdf

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


#=
d, r, α, θl, θh, w, θ = params_kinova_gen3()

using Symbolics
@variables t[1:3] a[1:3] d[1:3] r[1:3]


z1=[cos(t[1]) -sin(t[1]) 0 0; sin(t[1]) cos(t[1]) 0 0; 0 0 1 d[1]; 0 0 0 1]
x1=[1 0 0 r[1]; 0 cos(a[1]) -sin(a[1]) 0; 0 sin(a[1]) cos(a[1]) 0; 0 0 0 1]

z2=[cos(t[2]) -sin(t[2]) 0 0; sin(t[2]) cos(t[2]) 0 0; 0 0 1 d[2]; 0 0 0 1]
x2=[1 0 0 r[2]; 0 cos(a[2]) -sin(a[2]) 0; 0 sin(a[2]) cos(a[2]) 0; 0 0 0 1]

z3=[cos(t[3]) -sin(t[3]) 0 0; sin(t[3]) cos(t[3]) 0 0; 0 0 1 d[3]; 0 0 0 1]
x3=[1 0 0 r[3]; 0 cos(a[3]) -sin(a[3]) 0; 0 sin(a[3]) cos(a[3]) 0; 0 0 0 1]

t1 = z1*x1
t2 = z2*x2
t3 = z3*x3
=#