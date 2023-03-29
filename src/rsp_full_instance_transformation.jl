using JuMP, Gurobi
using Combinatorics

include("dat.jl")
include("misc.jl")


main = Model(optimizer_with_attributes(Gurobi.Optimizer))
# Decision variables
@variable(main, x[i in V, j in V; i<j], Bin)
@variable(main, y[i in V,j in V], Bin)
@variable(main, x_prime[i in V,j in V; i<j], Bin)
@variable(main, sigma)

# # Objective function
@objective(main, Min, sum(rc[i,j]*x[i,j]+ backup[i,j]*x_prime[i,j] for (i,j) in E)+ sum(oc[i]*y[i,i] for i in V)+ sum(sc[i,j]*y[i,j] for (i,j) in A) + offset) 

# comment
@constraint(main, degree_constr[i in V] ,sum(x[minmax(i,j)] for j in V if i!=j)==  2*y[i,i])

@constraint(main, terminal_constr[i in V], sum(3*y[i,j] for j in V_certain if i!=j) + sum(y[i,j] for j in V_tilt if j!=i) == 3*(1-y[i,i]))
@constraint(main, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
@constraint(main, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
@constraint(main, [(i,j,k,t) in J_tilt], x[minmax(k,i)]+ x[minmax(i,j)] + x[minmax(j,t)] <= 2+x_prime[k,t])
@constraint(main, [(i,j,k) in K_tilt], x[minmax(i,j)] + x[minmax(j,k)] <= 1+x_prime[minmax(i,k)])

@constraint(main, [(i,j) in A], y[i,j] <= y[j,j])
@constraint(main, y[1,1] == 1)

function my_callback_benders_cut(cb_data)
    x_hat = Bool.(round.(callback_value.(cb_data, x)))
    y_hat = Bool.(round.(callback_value.(cb_data, y)))
    x_hat = _transform_matrix(x_hat)
    status = callback_node_status(cb_data, main)
    
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        # Check subtour in a tour
        if contain_subtour(x_hat, y_hat) == 0
            _list_hub = [i for i in 1:n if y_hat[i,i] == 1]
            for i in 2:ceil(Int,length(_list_hub)/2)
                for S in combinations(_list_hub, i)
                    con = @build_constraint(length(S) - 1/length(V)*sum(y[i,i] for i in S)>= sum(x[i,j] for i in S for j in S if i<j))
                    MOI.submit(main, MOI.LazyConstraint(cb_data), con)
                end
            end
        end
    end
end

set_attribute(main, MOI.LazyConstraintCallback(), my_callback_benders_cut)
optimize!(main)
@show value.(x);
@show value.(y);
@show value.(x_prime);