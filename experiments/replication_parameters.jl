using LinearAlgebra

function rand_between(n...; low, high)
	(high - low) .* rand(n...) .+ low
end

function random_manipulator(n_joints)
	r = rand_between(n_joints; low=0.1, high=1)
	d = rand_between(n_joints; low=0.1, high=1)
	α = rand_between(n_joints; low=-3, high=3)
	α = rand_between(n_joints; low=-3, high=3)
	α = rand([-π/2, π/2], n_joints)

	θl = vec(fill(-3.0,n_joints))
	θh = vec(fill(+3.0,n_joints))

	r, d, α, θl, θh
end

r, d, α, θl, θh = random_manipulator(7)
w = normalize(ones(7), 1)
θ = zeros(7)