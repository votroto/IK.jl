using Test

import IK

function _random_feasible_pose_hq(d, r, α, θl, θh)
    x = θl .+ rand(length(θl)) .* (θh .- θl)

    prod(IK.dh_matrix.(x, d, α, r)), prod(IK.dh_quaternion.(x, d, α, r))
end

function test_quat_trans_rot_part_agree(n)
    ls, hs = -3 * rand(n), 3 * rand(n)
    H, Q = _random_feasible_pose_hq(rand(n), rand(n), randn(n), ls, hs)
    rh = H[1:3, 1:3]

    @test isapprox(IK._quaternion_to_rotation(Q.r), rh; atol=1e-8)
end

@testset "quaternion" begin
    for i in 1:100
        test_quat_trans_rot_part_agree(3)
    end
end