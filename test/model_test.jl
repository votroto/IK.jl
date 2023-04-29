using Test
using IK: trig_lb_infeas, trig_ub_infeas

function gen_trig_split()
    bound = rand_between(-pi+eps(), pi-eps())
    below = rand_between(10, lo=-pi+eps(), hi=bound)
    above = rand_between(10, lo=bound, hi=pi-eps()) 

    bound, below, above
end

function test_lb_infeas()
    bnd, below, above = gen_trig_split()
    
    sat(x) = trig_lb_infeas(cos(x), sin(x), bnd)
    @test all(x-> sat(x) >= 0, below)
    @test all(x-> sat(x) <= 0, above)
end

function test_ub_infeas()
    bnd, below, above = gen_trig_split()
    
    sat(x) = trig_ub_infeas(cos(x), sin(x), bnd)
    @test all(x-> sat(x) >= 0, above)
    @test all(x-> sat(x) <= 0, below)
end

@testset "trig_infeas" begin
    for i in 1:100
            test_lb_infeas()
            test_ub_infeas()
    end
end