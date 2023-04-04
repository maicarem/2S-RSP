using JuMP, Gurobi
using Combinatorics

include("dat.jl")
include("misc.jl")
include("write_output.jl")

# Initialize Sets
V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt = _declare_set(n, 1)
opening_cost, ring_cost, star_cost = oc, rc, sc
offset, oc, rc, sc, backup = _transformation_cost(rc,sc, oc, n, V_tilt, V_certain)


main = Model(optimizer_with_attributes(Gurobi.Optimizer))
# Decision variables
@variable(main, x[i in V, j in V; i<j], Bin)
@variable(main, y[i in V,j in V], Bin)
@variable(main, x_prime[i in V,j in V; i<j], Bin)
@variable(main, sigma)

# # Objective function
@objective(main, Min, sum(ring_cost[i,j]*(x[i,j]+x_prime[i,j]) for (i,j) in E)+ sum(opening_cost[i]*y[i,i] for i in V)+ sum(star_cost[i,j]*y[i,j] for (i,j) in A))

# comment
@constraint(main, degree_constr[i in V], sum(x[minmax(i,j)] for j in V if i<j) + sum(x[minmax(j,i)] for j in V if i>j)==  2*y[i,i])

@constraint(main, terminal_constr[i in V], sum(3*y[i,j] for j in V_certain if i!=j) + sum(y[i,j] for j in V_tilt if j!=i) == 3*(1-y[i,i]))
@constraint(main, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
@constraint(main, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
@constraint(main, [(i,j,k,t) in J_tilt], x[minmax(k,i)]+ x[minmax(i,j)] + x[minmax(j,t)] <= 2+x_prime[minmax(k,t)])
@constraint(main, [(i,j,k) in K_tilt], x[minmax(i,j)] + x[minmax(j,k)] <= 1+x_prime[minmax(i,k)])

@constraint(main, [(i,j) in A], y[i,j] <= y[j,j])
@constraint(main, y[1,1] == 1)

function my_callback_benders_cut(cb_data)
    
    x_hat = Bool.(round.(callback_value.(cb_data, x)))
    y_hat1 = Bool.(round.(callback_value.(cb_data, y)))
    x_hat = _transform_matrix(x_hat)
    
    status = callback_node_status(cb_data, main)
    
    y_hat = zeros(Bool, n)
    for i in 1:n
        y_hat1[i,i] == 1 || continue
        y_hat[i] = 1
    end
    
    all_cycles = find_cycle(x_hat, y_hat)

    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        # Check subtour in a tour
        if length(all_cycles) > 1
            _list_hub = [i for i in 1:n if y_hat[i] == 1]
            # add subtour elimination
            for each_cycle in all_cycles
                con = @build_constraint(length(each_cycle) - 1/(length(_list_hub)- length(each_cycle))*sum(y[i,i] for i in _list_hub if i ∉ each_cycle)>= 
                sum(x[minmax(each_cycle[i], each_cycle[i+1])] for i in eachindex(each_cycle[1:end-1]))+ x[minmax(each_cycle[1], each_cycle[end])])
                MOI.submit(main, MOI.LazyConstraint(cb_data), con)
            end
        end
    end
end

set_attribute(main, MOI.LazyConstraintCallback(), my_callback_benders_cut)
optimize!(main)
_write_gurobi_log("instance", main, "Original")
@show value.(x)
@show value.(y);
@show value.(x_prime[1,2]);