import GLPK
import MathOptInterface as MOI

result_dict = Dict()
result_dict["algorithm"] = "bbc"
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

include("dual_solution.jl")
include("misc.jl")

lb_distance = _find_lower_bound_backup(n, V_tilt, rc)
master_bbc = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 1))
set_time_limit_sec(master_bbc, pars.time_limit)

@variable(master_bbc, x[i in V, j in V; i<j], Bin)
@variable(master_bbc, y[i in V], Bin)
@variable(master_bbc, lambda_0>=0)
@variable(master_bbc, lambda[i in V]>=0)

if pars.transformation
    @objective(master_bbc, Min, offset + sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))
else
    @objective(master_bbc, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))
end

@constraint(master_bbc, degree_constr[i in V] ,sum(x[minmax(i,j)] for j in V if i!=j)==  2*y[i])
@constraint(master_bbc, sum(x[i,j] for (i,j) in E) >= 6)
@constraint(master_bbc, y[1] == 1)

if length(V_tilt) >= 1
    @constraint(master_bbc, lambda_0 >= sum(y[i]* lb_distance[i] for i in V_tilt))
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
                
                lambda_0_hat = round(callback_value(cb_data, lambda_0))
                lambda_hat = round.(callback_value.(cb_data, lambda))
                
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
                
                if pars.transformation
                    lower_bound = offset + sum(rc[i,j]*x_hat[i,j] for (i,j) in E)+ sum(oc[i]*y_hat[i] for i in V) + lambda_0_hat + sum(lambda_hat[i] for i in V)
                else
                    lower_bound = sum(ring_cost[i,j]*x_hat[i,j] for (i,j) in E)+ sum(opening_cost[i]*y_hat[i] for i in V)+ lambda_0_hat + sum(lambda_hat[i] for i in V)
                end
                
                gap = (upper_bound - lower_bound)/upper_bound
                # ADD SP0
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
                result_dict["time_sp"] = time() - time_sp
                result_dict["num_hubs"] = floor(Int, sum(y_hat[i] for i in 1:n))
            end

                # add cut SP0

        end
    end
    
    # # # # # # # # # # # # 
    # MAIN PROGRAM
    # # # # # # # # # # # # 
    set_attribute(master_bbc, MOI.LazyConstraintCallback(), my_callback_benders_cut)
    optimize!(master_bbc)

    result_dict["num_constraint_ilp_including_integrality"] = num_constraints(master_bbc; count_variable_in_set_constraints = true)
    result_dict["num_constraint_ilp_notinclude_integrality"] = num_constraints(master_bbc; count_variable_in_set_constraints = false)
    result_dict["lower_bound"] = objective_bound(master_bbc)
    result_dict["upper_bound"] = objective_value(master_bbc)
    result_dict["total_time"] = solve_time(master_bbc)
    result_dict["timestamp"] = now()

    @show result_dict
    _write_log_bbc(name, master_bbc, n, V_tilt, "bbc", pars)
end

main_rsp_bbc()