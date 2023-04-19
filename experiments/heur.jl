using JuMP
using Gurobi

include("replication_parameters.jl")
include("../src/ahomo_kinematics.jl")
include("../src/forward_kinematics.jl")
include("../src/inverse_kinematics.jl")

r, d, α, θl, θh = random_manipulator(6)
w = normalize(ones(6), 1)
θ = zeros(6)
M = random_feasible_pose(d, r, α, θl, θh)
solve_inverse_kinematics(d, r, α, θl, θh, M, θ, w)


# initialize the variables
@var z[1:6,1:3] p[1:3]
α = randn(5)
a = randn(9)
# define the system of polynomials
f = [z[i,:] ⋅ z[i,:] for i = 2:5]
g = [z[i,:] ⋅ z[i+1,:] for i = 1:5]
h = sum(a[i] .* (z[i,:] × z[i+1,:]) for i=1:3) + sum(a[i+4] .* z[i,:] for i = 2:5)
F_ = [f .- 1; g .- cos.(α); h .- p]
# assign values to z₁ and z₆
z1 = normalize!(randn(3))
z6 = normalize!(randn(3))
F = [subs(f, z[1,:]=>z1, z[6,:]=>z6, p=>[1, 1, 0]) for f in F_]

variable_groups = [[z[2,:]; z[4,:]], [z[3,:]; z[5,:]]]
solve(System(F; variable_groups=variable_groups); start_system = :total_degree)





using DynamicPolynomials
# initialize the variables
@polyvar z[1:6, 1:3] p[1:3]
# define the system of polynomials
f = [z[i, :] ⋅ z[i, :] for i = 2:5]
g = [z[i, :] ⋅ z[i+1, :] for i = 1:5]
h = sum(a[i] .* (z[i, :] × z[i+1, :]) for i = 1:3) + sum(a[i+4] .* z[i, :] for i = 2:5)
F_ = [f .- 1; g .- cos.(α); h .- p]
E = [subs(f, z[1, :] => z1, z[6, :] => z6, p => [1, 1, 0]) for f in F_]


optimizer = _default_optimizer_poly()
m = Model(optimizer)

@variable(m, jz[1:6, 1:3])
@variable(m, jp[1:3])

E = map(e -> e([z[:]; p[:]] => [jz[:]; jp[:]]), E)

@constraint m E .== 0

optimize!(m)

#=

m = Model(Gurobi.Optimizer)

@variable m -1 <= x <= 1
@variable m -1 <= y <= 1

#@variable m xx
xx = jump_quadratic_envelope(x,x)
@constraint m xx == x * x
xx1 = jump_quadratic_envelope(x,x)
#@variable m yy
yy = jump_quadratic_envelope(x,x)
@constraint m yy == y * y

@constraint m xx + yy == 0.8
@constraint m 5*x*y - 2*xx1 - 2*x*yy - y == 0

@objective m Min -x + 0.5 + 0.5*y - 0.5

#MOI.set(m, MOI.RawOptimizerAttribute("Method"), 1)
MOI.set(m, MOI.RawOptimizerAttribute("NonConvex"), 2)
optimize!(m)
@show value(x)
@show value(y)
=#
#=


lving as a MIP

Presolve time: 0.00s
Presolved: 12 rows, 7 columns, 26 nonzeros
Presolved model has 4 bilinear constraint(s)
Variable types: 7 continuous, 0 integer (0 binary)

Root relaxation: objective -7.600000e+00, 2 iterations, 0.00 seconds (0.00 work units)

    Nodes    |    Current Node    |     Objective Bounds      |     Work
 Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time

     0     0   -7.60000    0    4          -   -7.60000      -     -    0s
H    0     0                      -3.6579366   -7.60000   108%     -    0s
     0     0   -5.80000    0    4   -3.65794   -5.80000  58.6%     -    0s
     0     0   -5.50000    0    4   -3.65794   -5.50000  50.4%     -    0s
     0     0   -5.40411    0    4   -3.65794   -5.40411  47.7%     -    0s
     0     0   -4.83691    0    3   -3.65794   -4.83691  32.2%     -    0s
     0     0   -4.73663    0    4   -3.65794   -4.73663  29.5%     -    0s
     0     0   -4.54269    0    4   -3.65794   -4.54269  24.2%     -    0s
     0     0   -4.53663    0    4   -3.65794   -4.53663  24.0%     -    0s
     0     0   -4.53663    0    4   -3.65794   -4.53663  24.0%     -    0s
H    0     0                      -3.7385825   -4.53663  21.3%     -    0s
     0     2   -4.53663    0    4   -3.73858   -4.53663  21.3%     -    0s
H   12     3                      -3.7388652   -3.74964  0.29%   3.2    0s

Cutting planes:
  RLT: 3
  PSD: 2

Explored 18 nodes (66 simplex iterations) in 0.01 seconds (0.00 work units)
Thread count was 4 (of 4 available processors)

Solution count 3: -3.73887 -3.73858 -3.65794

Optimal solution found (tolerance 1.00e-04)
Best objective -3.738865243624e+00, best bound -3.739062504539e+00, gap 0.0053%

User-callback calls 244, time in user-callback 0.00 sec
value(x) = -0.8499393771023316
value(y) = 0.5268804942782571
0.5268804942782571






function my_callback_function(cb_data)
	vars = [x, y]
	vals = [0.997, 0.003]
                c_val = callback_value(cb_data,x)
@show c_val
@show status = callback_node_status(cb_data, m)


        acc = MOI.submit(m, MOI.HeuristicSolution(cb_data), vars, vals)

	# Should be always accepted, but it never is unless queried at MIPNODE.
        #if acc == MOI.HEURISTIC_SOLUTION_ACCEPTED
        #        throw("yay")
        #end
end

MOI.set(m, MOI.RawOptimizerAttribute("NonConvex"), 2)
MOI.set(m, MOI.HeuristicCallback(), my_callback_function)
optimize!(m)
=#