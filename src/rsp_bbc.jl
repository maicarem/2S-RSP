import GLPK
import MathOptInterface as MOI

include("dual_solution.jl")
include("misc.jl")

result_dict = initializeResult_dict("bbc")


lb_distance = _find_lower_bound_backup(n, V_tilt, rc)
master_bbc = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 1))
set_time_limit_sec(master_bbc, pars.time_limit)

@variable(master_bbc, x[i in V, j in V; i<j], Bin)
@variable(master_bbc, y[i in V], Bin)
if pars.one_cut
    @variable(master_bbc, lambda_one >=0)
else
    @variable(master_bbc, lambda_0>=0)
    @variable(master_bbc, lambda[i in V]>=0)
end

if pars.one_cut
    if pars.transformation
        @objective(master_bbc, Min, offset + sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i] for i in V)+ lambda_one)
    else
        @objective(master_bbc, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i] for i in V)+ lambda_one)
    end
else
    if pars.transformation
        @objective(master_bbc, Min, offset + sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))
    else
        @objective(master_bbc, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))
    end
end

@constraint(master_bbc, degree_constr[i in V] ,sum(x[minmax(i,j)] for j in V if i!=j)==  2*y[i])
@constraint(master_bbc, sum(x[i,j] for (i,j) in E) >= 6)
@constraint(master_bbc, y[1] == 1)

if pars.one_cut == false
    if length(V_tilt) >= 1
        @constraint(master_bbc, lambda_0 >= sum(y[i]* lb_distance[i] for i in V_tilt))
    end
end

function main_rsp_bbc()
    global_upper_bound = 1e18
    function my_callback_benders_cut(cb_data)
        x_hat = Bool.(round.(callback_value.(cb_data, x)))
        y_hat = Bool.(round.(callback_value.(cb_data, y)))
        x_hat = _transform_matrix(x_hat)
        
        status = callback_node_status(cb_data, master_bbc)
        all_cycles = find_cycle(x_hat, y_hat)

        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            if length(all_cycles) > 1
                _list_hub = [i for i in 1:n if y_hat[i] == 1]
                # add subtour elimination
                for each_cycle in all_cycles
                    con = @build_constraint(length(each_cycle) - 1/(length(_list_hub)- length(each_cycle))*sum(y[i] for i in _list_hub if i ∉ each_cycle)>= 
                    sum(x[minmax(each_cycle[i], each_cycle[i+1])] for i in eachindex(each_cycle[1:end-1]))+ x[minmax(each_cycle[1], each_cycle[end])])
                    MOI.submit(master_bbc, MOI.LazyConstraint(cb_data), con)
                    result_dict["num_subtour"] += 1
                end
            elseif length(all_cycles) == 1
                
                if pars.one_cut
                    lambda_one_hat = round(callback_value(cb_data, lambda_one))
                else
                    lambda_0_hat = round(callback_value(cb_data, lambda_0))
                    lambda_hat = round.(callback_value.(cb_data, lambda))
                end
                
                time_sp = time()
                if pars.transformation
                    (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat, V_tilt, n, backup, ring_cost, sc)
                else
                    (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat, V_tilt, n, ring_cost, ring_cost, star_cost)
                end

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
                
                if pars.one_cut
                    if pars.transformation
                        lower_bound = offset + sum(rc[i,j]*x_hat[i,j] for (i,j) in E)+ sum(oc[i]*y_hat[i] for i in V) + lambda_one_hat
                    else
                        lower_bound = sum(ring_cost[i,j]*x_hat[i,j] for (i,j) in E)+ sum(opening_cost[i]*y_hat[i] for i in V)+ lambda_one_hat
                    end
                else
                    if pars.transformation
                        lower_bound = offset + sum(rc[i,j]*x_hat[i,j] for (i,j) in E)+ sum(oc[i]*y_hat[i] for i in V) + lambda_0_hat + sum(lambda_hat[i] for i in V)
                    else
                        lower_bound = sum(ring_cost[i,j]*x_hat[i,j] for (i,j) in E)+ sum(opening_cost[i]*y_hat[i] for i in V)+ lambda_0_hat + sum(lambda_hat[i] for i in V)
                    end
                end
                
                gap = (upper_bound - lower_bound)/upper_bound
                # ADD SP0
                if pars.one_cut
                    if !(lambda_one_hat >= obj_sp0 + obj_spi)
                        cut = @build_constraint(lambda_one >= sum((x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ x[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* beta[(val[1],val[2],val[3])] for val in keys(beta) if length(keys(beta))>0) + sum((x[minmax(val[3],val[1])[1],minmax(val[3],val[1])[2]] + x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+x[minmax(val[2],val[4])[1],minmax(val[2],val[4])[2]] - 2) * alpha[(val[1],val[2],val[3],val[4])] for val in keys(alpha) if length(keys(alpha))>0)
                                                + sum(3(1-y[indice])*φ[indice] for indice in 1:n) - sum(y[j]*γ[indice,j] for indice in 1:n for j in 1:n if j!=indice))
                        MOI.submit(master_bbc, MOI.LazyConstraint(cb_data), cut)
                        result_dict["num_cut_sp0"] += 1
                    end
                else
                    if !(lambda_0_hat >= obj_sp0)
                        cut = @build_constraint(lambda_0 >= sum((x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ x[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* beta[(val[1],val[2],val[3])] for val in keys(beta)) + sum((x[minmax(val[3],val[1])[1],minmax(val[3],val[1])[2]] + x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+x[minmax(val[2],val[4])[1],minmax(val[2],val[4])[2]] - 2) * alpha[(val[1],val[2],val[3],val[4])] for val in keys(alpha)))
                        MOI.submit(master_bbc, MOI.LazyConstraint(cb_data), cut)
                        result_dict["num_cut_sp0"] += 1
                    end
                    # Add SPi
                    for i in 1:n
                        y_hat[i] == 0 && φ[i]!= 0|| continue
                        if !(lambda_hat[i] >= 3(1-y_hat[i])*φ[i] - sum(y_hat[j]*γ[i,j] for j in V if j!=i))
                            con = @build_constraint(lambda[i] >= 3(1-y[i])*φ[i] - sum(y[j]*γ[i,j] for j in V if j!=i))
                            MOI.submit(master_bbc, MOI.LazyConstraint(cb_data), con)
                            result_dict["num_cut_spi"] += 1     
                        end
                    end
                end
                result_dict["time_sp"] = time() - time_sp
            end
        end
    end
    
    # # # # # # # # # # # # 
    # MAIN PROGRAM
    # # # # # # # # # # # # 
    set_attribute(master_bbc, MOI.LazyConstraintCallback(), my_callback_benders_cut)
    set_optimizer_attribute(master_bbc, "Cutoff", bf_global_obj)
    optimize!(master_bbc)
    if has_values(master_bbc)
        result_dict["num_constraint_ilp_including_integrality"] = num_constraints(master_bbc; count_variable_in_set_constraints = true)
        result_dict["num_constraint_ilp_notinclude_integrality"] = num_constraints(master_bbc; count_variable_in_set_constraints = false)
        result_dict["lower_bound"] = objective_bound(master_bbc)
        result_dict["upper_bound"] = objective_value(master_bbc)
        result_dict["num_hubs"] = sum(value(y[i,i]) for i in 1:n)
        result_dict["total_time"] = solve_time(master_bbc)
        result_dict["timestamp"] = now()
    else
        result_dict["lower_bound"] = objective_bound(master_bbc)
        result_dict["upper_bound"] = bf_global_obj
        result_dict["num_hubs"] = bf_global_num_hubs
        result_dict["total_time"] = solve_time(master_bbc)
        result_dict["timestamp"] = now()
    end
    result_dict["node_count"] = MOI.get(master_bbc, MOI.NodeCount())
    @info "Node count = $(result_dict["node_count"])"
    result_dict["gap"] = (result_dict["upper_bound"] - result_dict["lower_bound"])/result_dict["upper_bound"]
    @info "Completed ILP ...Gap = $(result_dict["gap"]), Lower_bound = $(result_dict["lower_bound"])"
    write_ouput(pars, name, result_dict, MainPar)
end

main_rsp_bbc()