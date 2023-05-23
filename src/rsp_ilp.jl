using JuMP, Gurobi
using Combinatorics
using Graphs, GraphsFlows

include("dat.jl")
include("misc.jl")
include("write_output.jl")
include("mutable_structure.jl")
include("user_cut.jl")
include("brute_force.jl")

starting_time = time()

# Initialize Sets
pars = MainPar(uc_strat = 3, transformation = true, alpha = 3)
name = "instances/small_instances/small_instance_10.dat"
n, oc, sc, rc = read_input_random(name, pars)
V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt = _declare_set(n, pars)
opening_cost, ring_cost, star_cost = oc, rc, sc
if pars.transformation
    offset, oc, rc, sc, backup = _transformation_cost(rc,sc, oc, n, V_tilt, V_certain)
end

starting_time = time()

objval_3, objsol_3, objtime_3 = three_hubs(n, V_certain, rc, sc, oc)
objval_4, objsol_4, objtime_4 = four_hubs(n, V_tilt, V_certain, rc, sc, oc)
objval_5, objsol_5, objtime_5 = five_hubs(n, V_tilt, rc, oc, sc)

main = Model(optimizer_with_attributes(Gurobi.Optimizer))

# Decision variables
@variable(main, x[i in V, j in V; i<j], Bin)
@variable(main, y[i in V,j in V], Bin)
@variable(main, x_prime[i in V,j in V; i<j], Bin)

if pars.transformation
    @variable(main, x_prime_prime[i in V,j in V; i<j], Bin)
else
    @variable(main, sigma>= 0, Int)
end

# OBJECTIVE FUNCTION
if pars.transformation
    @objective(main, Min, sum(rc[i,j]*x[i,j]+ backup[i,j]*x_prime[i,j] + ring_cost[i,j]*x_prime_prime[i,j] for (i,j) in E)+ sum(oc[i]*y[i,i] for i in V)+ sum(sc[i,j]*y[i,j] for (i,j) in A) + offset) 
else
    @objective(main, Min, sum(ring_cost[i,j]*(x[i,j] + x_prime[i,j]) for (i,j) in E)+ sum(opening_cost[i]*y[i,i] for i in V)+ sum(star_cost[i,j]*y[i,j] for (i,j) in A))
end

# CONSTRAINT
@constraint(main, degree_constr[i in V], sum(x[minmax(i,j)] for j in V if i<j) + sum(x[minmax(j,i)] for j in V if i>j)==  2*y[i,i])
@constraint(main, terminal_constr[i in V], sum(3*y[i,j] for j in V_certain if i!=j) + sum(y[i,j] for j in V_tilt if j!=i) == 3*(1-y[i,i]))
@constraint(main, [(i,j,k) in K_tilt], x[minmax(i,j)] + x[minmax(j,k)] <= 1+x_prime[minmax(i,k)])

if pars.transformation
    @constraint(main, [(i,j,k,t) in J_tilt], x[minmax(k,i)]+ x[minmax(i,j)] + x[minmax(j,t)] <= 2+x_prime_prime[minmax(k,t)])
    @constraint(main, sum(x[i,j] for (i,j) in E) >= 6)
else
    @constraint(main, sum(x[i,j] for (i,j) in E) >= 3 + sigma)
    @constraint(main, [(i,j,k,t) in J_tilt], x[minmax(k,i)]+ x[minmax(i,j)] + x[minmax(j,t)] <= 2+x_prime[minmax(k,t)])
    @constraint(main, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
end
@constraint(main, [(i,j) in A], y[i,j] <= y[j,j])
@constraint(main, y[1,1] == 1)

set_attribute(main, MOI.LazyConstraintCallback(), my_callback_subtour)
if pars.uc
    set_attribute(main, MOI.UserCutCallback(), call_back_user_cuts)
end

optimize!(main)
@show time() - starting_time
@show objval_3, objval_4, objval_5
@show objtime_3, objtime_4, objtime_5