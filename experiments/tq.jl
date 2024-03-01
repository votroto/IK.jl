include("../src/quaternion.jl")
include("../src/denavit_hartenberg.jl")
include("kuka_parameters.jl")

d, r, α, θl, θh, w, θ = params_kuka_iiwa()
xa = [0.8931223845580142, 0.13110929613780642, 0.22603885976560512, 0.24755798425962494, 0.31062843917916094, 0.18434542093334116, 0.471369453581177]
M = DualQuaternion(
    Quaternion(0.8940641540612714, SA[0.7724486577875048, 0.7224934216516122, 0.642525736772827]),
    Quaternion(0.4976036764574866, SA[0.18284795540546184, 0.9046750963077218, 0.04265846807999851])
)

m1 = dq(xa[1], d[1], α[1], r[1])
m2 = dq(xa[2], d[2], α[2], r[2])
m3 = dq(xa[3], d[3], α[3], r[3])
m4 = dq(xa[4], d[4], α[4], r[4])
m5 = dq(xa[5], d[5], α[5], r[5])
m6 = dq(xa[6], d[6], α[6], r[6])
m7 = dq(xa[7], d[7], α[7], r[7])



res1 = m1 * m2 * m3 * m4 * m5 * m6 * m7 


res3 = m1 * m2 * m3 * m4 
res4 = res1 * m7' * m6' * m5'