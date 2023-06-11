include("misc.jl")
include("dual_solution.jl")

result_dict = Dict()
result_dict["algorithm"] = "benders"
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

lb_distance = _find_lower_bound_backup(n, V_tilt, rc)
############# master_bender PROBLEM ########################################
master_bender = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    # Decision variables
@variable(master_bender, x[i in V, j in V; i<j], Bin)
@variable(master_bender, y[i in V], Bin)
@variable(master_bender, sigma>=0, Int)
if pars.one_cut
    @variable(master_bender, lambda_one>=0)
else
    @variable(master_bender, lambda_0>=0)
    @variable(master_bender, lambda[i in V]>=0)
end

function rsp_benders(pars, n, oc, sc, rc, backup, lb_distance, opening_cost, ring_cost, star_cost, V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt, result_dict)
    
    starting_time = time()
    # Objective Function
    if pars.one_cut
        if pars.transformation
            @objective(master_bender, Min, offset + sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i] for i in V)+ lambda_one)
        else
            @objective(master_bender, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i] for i in V)+ lambda_one)
        end
    else
        if pars.transformation
            @objective(master_bender, Min, offset + sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))
        else
            @objective(master_bender, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))
        end
    end

    # Constraint
    @constraint(master_bender, degree_constr[i in V] ,sum(x[minmax(i,j)] for j in V if i!=j)==  2*y[i])
    @constraint(master_bender, sum(x[i,j] for (i,j) in E) >= 6)
    @constraint(master_bender, y[1] == 1)
    
    if length(V_tilt) >= 1 && pars.one_cut == false
        @constraint(master_bender, lambda_0 >= sum(y[i]* lb_distance[i] for i in V_tilt))
    end
    
    global_upper_bound = 1e18
    iter0 = 1
    
    while time() - starting_time < pars.time_limit
        
        iter0 +=1
        println("===============Iteration ", iter0, "===============")
        
        function my_callback_subtour(cb_data)
            
            x_hat = Bool.(round.(callback_value.(cb_data, x)))
            y_hat = Bool.(round.(callback_value.(cb_data, y)))
            x_hat = _transform_matrix(x_hat)
            status = callback_node_status(cb_data, master_bender)
            all_cycles = find_cycle(x_hat, y_hat)

            if status == MOI.CALLBACK_NODE_STATUS_INTEGER
                # Check subtour in a tour
                if length(all_cycles) > 1
                    _list_hub = [i for i in 1:n if y_hat[i] == 1]
                    # add subtour elimination
                    for each_cycle in all_cycles
                        con = @build_constraint(length(each_cycle) - 1/(length(_list_hub)- length(each_cycle))*sum(y[i] for i in _list_hub if i ∉ each_cycle)>= 
                        sum(x[minmax(each_cycle[i], each_cycle[i+1])] for i in eachindex(each_cycle[1:end-1]))+ x[minmax(each_cycle[1], each_cycle[end])])
                        MOI.submit(master_bender, MOI.LazyConstraint(cb_data), con)
                        result_dict["num_subtour"] +=1
                    end
                end
            end
        end

        set_attribute(master_bender, MOI.LazyConstraintCallback(), my_callback_subtour)
        # set_attribute(master_bender, MOI.UserCutCallback(), call_back_user_cuts_benders)
        master_time = time()
        optimize!(master_bender)
        result_dict["time_master"] += time() - master_time

        lower_bound = objective_value(master_bender)
        println("Objective value at iteration $(iter0) is $(lower_bound)")
        x_hat_1, y_hat = Bool.(round.(value.(x))), Bool.(round.(value.(y)))
        x_hat = _transform_matrix(x_hat_1)
        
        if pars.one_cut
            lambda_one_hat = value(lambda_one)
        else
            lambda_0_hat, lambda_hat = value(lambda_0), round.(value.(lambda))
        end
        
        sp_time = time()
        if pars.transformation
            (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat, V_tilt, n, backup, ring_cost, sc)
        else
            (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat, V_tilt, n, ring_cost, ring_cost, star_cost)
        end
        
        # Objective value, and add cut
        obj_sp0 = cal_obj_sp0(alpha, beta, x_hat)
        obj_spi = cal_obj_spi(φ, γ, y_hat)
        
        if pars.transformation
            upper_bound = offset + sum(rc[i,j]*x_hat[i,j] for (i,j) in E)+ sum(oc[i]*y_hat[i] for i in V) + obj_sp0 + obj_spi
        else
            upper_bound =  sum(ring_cost[i,j]*x_hat[i,j] for (i,j) in E)+ sum(opening_cost[i]*y_hat[i] for i in V) + obj_sp0 + obj_spi
        end
        
        if upper_bound < global_upper_bound
            global_upper_bound = upper_bound
            @info "New global upper bound = $(global_upper_bound)"
        end
        @show "Current upper bound $(upper_bound), Global Upper bound = $(global_upper_bound)"

        gap = (upper_bound - lower_bound)/upper_bound
        
        println("Upper bound is $(upper_bound)")
        
        
        if  gap < 1e-10
            println("This is optimal with objective $(lower_bound)")
            break
        end
        
        if pars.one_cut
            _add_one_cut(master_bender, obj_sp0, obj_spi, alpha, beta, lambda_one_hat, n, φ, γ, result_dict)
        else
            _add_cut_SP0(master_bender, alpha, beta, lambda_0_hat, obj_sp0, result_dict)
            for i in 1:n
                y_hat[i] == 0 && φ[i]!= 0|| continue
                _add_cut_SPi(master_bender, φ, γ, i, lambda_hat, y_hat, result_dict)
            end
        end

        result_dict["time_sp"] += time() - sp_time
        result_dict["total_time"] = time() - starting_time
        result_dict["num_hubs"] = floor(Int, sum(y_hat[i] for i in 1:n))
        result_dict["timestamp"] = now()
        result_dict["lower_bound"] = lower_bound
        result_dict["upper_bound"] = global_upper_bound
        write_ouput(pars, name, result_dict, MainPar)
        # _write_log_benders(name, n, V_tilt, "benders", lower_bound, global_upper_bound, time() - starting_time, pars)
    end

end

rsp_benders(pars, n, oc, sc, rc, backup, lb_distance, opening_cost, ring_cost, star_cost, V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt, result_dict)

