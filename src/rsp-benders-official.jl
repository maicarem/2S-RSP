using JuMP, Gurobi
using Combinatorics

include("dat.jl")
include("transformation-cost.jl")

function master_problem()
    # offset, opening_cost, ring_cost, star_cost, backup_cost = _transformation_cost()
    @show offset, opening_cost, ring_cost, star_cost, backup_cost
    master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 1))
    
    # Decision variables
    @variable(master, x[i in V, j in V], Bin)
    @variable(master, y[i in V,j in V], Bin)
    @variable(master, x_prime[i in V,j in V], Bin)
    @variable(master, sigma>=0, Int)
    @variable(master, lambda_0>=0)
    @variable(master, lambda[i in V]>=0)

    # @objective(master, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i,i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))
    @constraint(master, degree_constr[i in V] ,sum(x[i,j] for j in V if i<j) + sum(x[j,i] for j in V if i>j)==  2*y[i,i])

    for i in 2:ceil(Int,length(V)/2)
        for S in combinations(V, i)
            @constraint(master, length(S) - 1/length(V)*sum(y[i,i] for i in S)>= sum(x[i,j] for i in S for j in S if i<j))
        end
    end
    
    # @constraint(master, degree_constr[i in V] ,sum(x[i,j] for j in V if i<j) + sum(x[j,i] for j in V if i>j)==  2*y[i,i])
    @constraint(master, [(i,j) in E], x[i,j] == x[j,i])
    @constraint(master, [(i,j) in E], x_prime[i,j] == x_prime[j,i])

    @constraint(master, terminal_constr[i in V], sum(3*y[i,j] for j in V_certain if i!=j) + sum(y[i,j] for j in V_tilt if j!=i) == 3*(1-y[i,i]))
    @constraint(master, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
    @constraint(master, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
    @constraint(master, [(i,j,k,t) in J_tilt], x[k,i]+ x[i,j] + x[j,t] <= 2+x_prime[i,k])
    @constraint(master, [(i,j,k) in K_tilt], x[i,j] + x[j,k] <= 1+x_prime[i,k])

    @constraint(master, [(i,j) in A], y[i,j] <= y[j,j])
    @constraint(master, y[1,1] == 1)
    
    # offset, opening_cost, ring_cost, star_cost, backup_cost 
    @objective(master, Min, offset + sum(opening_cost[i]* y[i,i] for i in V) + sum(ring_cost[i,j]*x[i,j] + 
                                backup_cost[i,j]*x_prime[i,j] for (i,j) in E) + sum(star_cost[i,j]*y[i,j] for (i,j) in A))
    
    return master
end

master = master_problem()
master
