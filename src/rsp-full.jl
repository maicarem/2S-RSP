using JuMP, Gurobi
using Combinatorics

include("dat.jl")

main = Model(optimizer_with_attributes(Gurobi.Optimizer))
# Decision variables
@variable(main, x[i in V, j in V; i!=j], Bin)
@variable(main, y[i in V,j in V], Bin)
@variable(main, x_prime[i in V,j in V], Bin)
@variable(main, sigma)

# # Objective function
@objective(main, Min, sum(rc[i,j]*(x[i,j]+x_prime[i,j]) for (i,j) in E)+ sum(oc[i]*y[i,i] for i in V)+ sum(sc[i,j]*y[i,j] for (i,j) in A))

for i in 2:ceil(Int,length(V)/2)
    print(i)
    for S in combinations(V, i)
        @constraint(main, length(S) - 1/length(V)*sum(y[i,i] for i in S)>= sum(x[i,j] for i in S for j in S if i<j))
    end
end
# comment
@constraint(main, degree_constr[i in V] ,sum(x[i,j] for j in V if i<j) + sum(x[j,i] for j in V if i>j)==  2*y[i,i])
@constraint(main, [(i,j) in E], x[i,j] == x[j,i])
@constraint(main, [(i,j) in E], x_prime[i,j] == x_prime[j,i])

@constraint(main, terminal_constr[i in V], sum(3*y[i,j] for j in V_certain if i!=j) + sum(y[i,j] for j in V_tilt if j!=i) == 3*(1-y[i,i]))
@constraint(main, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
@constraint(main, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
@constraint(main, [(i,j,k,t) in J_tilt], x[k,i]+ x[i,j] + x[j,t] <= 2+x_prime[i,k])
@constraint(main, [(i,j,k) in K_tilt], x[i,j] + x[j,k] <= 1+x_prime[i,k])

@constraint(main, [(i,j) in A], y[i,j] <= y[j,j])
@constraint(main, y[1,1] == 1)


optimize!(main)
@show value.(x);
@show value.(y);
@show value.(x_prime);
