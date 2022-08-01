using LinearAlgebra

r = [0; 0; 0; 0; 0; 0; 0]
d = [340; 0; 400; 0; 400; 0; 126]
α = [-π / 2, π / 2, -π / 2, π / 2, -π / 2, π / 2, 0]
θl = [-2.9671; -2.0944; -2.9671; -2.0944; -2.9671; -2.0944; -3.0543]
θh = [2.9671; 2.0944; 2.9671; 2.0944; 2.9671; 2.0944; 3.0543]

w = normalize(ones(7))
θ = zeros(7)