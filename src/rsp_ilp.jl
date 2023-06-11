include("misc.jl")
include("user_cut.jl")


main = Model(optimizer_with_attributes(Gurobi.Optimizer))
set_time_limit_sec(main, pars.time_limit)
result_dict = Dict()
result_dict["algorithm"] = "ilp"
for param_name in ["num_constraint_ilp_including_integrality", 
                    "num_constraint_ilp_notinclude_integrality", 
                    "lower_bound", 
                    "upper_bound", 
                    "num_subtour", 
                    "num_hubs",
                    "time_sp",
                    "time_master",
                    "total_time",
                    "num_cut_sp0",
                    "num_cut_spi",
                    "obj_bf3",
                    "obj_bf4",
                    "obj_bf5",
                    "time_bf3",
                    "time_bf4",
                    "time_bf5",
                    "timestamp"]
    result_dict[param_name] = 0
end

# Decision variables
@variable(main, x[i in V, j in V; i<j], Bin)
@variable(main, y[i in V,j in V], Bin)
@variable(main, x_prime[i in V,j in V; i<j], Bin)

if pars.transformation
    @variable(main, x_prime_prime[i in V,j in V; i<j], Bin)
else
    @variable(main, sigma >= 0, Int)
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

function my_callback_subtour(cb_data)
    
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
                result_dict["num_subtour"] += 1
            end
        end
    end
end


set_attribute(main, MOI.LazyConstraintCallback(), my_callback_subtour)
if pars.uc
    set_attribute(main, MOI.UserCutCallback(), call_back_user_cuts)
end

optimize!(main)

result_dict["num_constraint_ilp_including_integrality"] = num_constraints(main; count_variable_in_set_constraints = true)
result_dict["num_constraint_ilp_notinclude_integrality"] = num_constraints(main; count_variable_in_set_constraints = false)
result_dict["lower_bound"] = objective_bound(main)
result_dict["upper_bound"] = objective_value(main)
result_dict["num_hubs"] = sum(value(y[i,i]) for i in 1:n)
result_dict["total_time"] = solve_time(main)
result_dict["timestamp"] = now()

@show name
write_ouput(pars, name, result_dict, MainPar)

# _write_gurobi_log(main, "ilp", name, pars, n, V_tilt, 0)