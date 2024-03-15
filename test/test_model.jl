using Test
using IK: lin_abs_angdiff_proxy, lin_angdiff_proxy

function gen_trig_split()
    bound = rand_between(-pi + eps(), pi - eps())
    below = rand_between(10, lo=-pi + eps(), hi=bound)
    above = rand_between(10, lo=bound, hi=pi - eps())

    bound, below, above
end

function test_angdiff_proxy_has_correct_sign()
    bnd, below, above = gen_trig_split()

    sat(x) = lin_angdiff_proxy(cos(x), sin(x), bnd)
    @test all(x -> sat(x) <= 0, below)
    @test all(x -> sat(x) >= 0, above)
end

function test_angdiff_proxy_on_same_angle_is_zero()
    xs = rand_between(10; lo=-pi + eps(), hi=pi - eps())

    diff(x) = lin_angdiff_proxy(cos(x), sin(x), x)
    @test all(x -> isapprox(diff(x), 0; atol=1e-8), xs)
end

function test_abs_angdiff_proxy_is_order_preserving()
    x = rand()
    ys = sort(rand(10); by=abs)

    obj(z) = lin_abs_angdiff_proxy(cos(z), sin(z), x, 1)
    @test issorted(x .+ ys; by=obj)
end

@testset "trig proxies" begin
    for i in 1:100
        test_angdiff_proxy_has_correct_sign()
        test_angdiff_proxy_on_same_angle_is_zero()
        test_abs_angdiff_proxy_is_order_preserving()
    end
end