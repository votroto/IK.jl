using LinearAlgebra

function continuous_uniform(n...; low, high)
    (high - low) .* rand(n...) .+ low
end

function random_manipulator(joint_count; twist_rng, max_angle)
    r = continuous_uniform(joint_count; low=0.1, high=1.0)
    d = continuous_uniform(joint_count; low=0.1, high=1.0)
    α = twist_rng(joint_count)

    θl = vec(fill(-max_angle, joint_count))
    θh = vec(fill(+max_angle, joint_count))

	w = normalize(ones(joint_count), 1)
	θ = zeros(joint_count)

    d, r, α, θl, θh, w, θ
end

function params_random_6rad(joint_count)
	twist_rng(n) = continuous_uniform(n; low=-3.0, high=3.0)
	random_manipulator(joint_count; twist_rng, max_angle=3.0)
end

function params_random_4rad(joint_count)
	twist_rng(n) = continuous_uniform(n; low=-3.0, high=3.0)
	random_manipulator(joint_count; twist_rng, max_angle=2.0)
end

function params_random_4rad_2alpha(joint_count)
	twist_rng(n) = continuous_uniform(n; low=-2.0, high=2.0)
	random_manipulator(joint_count; twist_rng, max_angle=2.0)
end

function params_random_orth(joint_count)
	twist_rng(n) = rand([-π / 2, π / 2], n)
	random_manipulator(joint_count; twist_rng, max_angle=3.0)
end
