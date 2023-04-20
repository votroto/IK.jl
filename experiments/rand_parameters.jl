using LinearAlgebra

function continuous_uniform(n...; low, high)
    (high - low) .* rand(n...) .+ low
end

function random_manipulator(joint_count; twist_rng, max_angle)
    r = continuous_uniform(joint_count; low=0.1, high=1)
    d = continuous_uniform(joint_count; low=0.1, high=1)
    α = twist_rng(joint_count)

    θl = vec(fill(-max_angle, joint_count))
    θh = vec(fill(+max_angle, joint_count))

	w = normalize(ones(joint_count), 1)
	θ = zeros(joint_count)

    r, d, α, θl, θh, w, θ
end

function params_random_6rad(joint_count)
	twist_rng(n) = continuous_uniform(n; low=-3, 3)
	random_manipulator(joint_count; twist_rng, 3)
end

function params_random_4rad(joint_count)
	twist_rng(n) = continuous_uniform(n; low=-3, 3) # check
	random_manipulator(joint_count; twist_rng, 2)
end

function params_random_orth(joint_count)
	twist_rng(n) = rand([-π / 2, π / 2], n)
	random_manipulator(joint_count; twist_rng, 3)
end
