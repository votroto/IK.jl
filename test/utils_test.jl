using Test

using IK: cos_min_max, replicate

function test_cos_min_max_symmetric()
    l = rand_between(-pi+eps(), 0)
    h = rand_between(0, pi-eps())
    @test cos_min_max(l, h) == (min(cos(l), cos(h)), 1.0)
end

function test_cos_min_max_single_sided()
    sign = rand([-1, 1])
    l = sign * rand_between(0+eps(), 0.5pi-eps())
    h = sign * rand_between(0.5pi+eps(), pi-eps())
    @test cos_min_max(l, h) == (cos(h), cos(l))
end

@testset "cos_min_max" begin
    for i in 1:100
            test_cos_min_max_symmetric()
            test_cos_min_max_single_sided()
    end
end