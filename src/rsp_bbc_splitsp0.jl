using JuMP, Gurobi
using Combinatorics
using Graphs, GraphsFlows

import GLPK
import MathOptInterface as MOI
import Dates

include("dat.jl")
include("dual_solution.jl")
include("misc.jl")
include("mutable_structure.jl")
include("user_cut.jl")


pars = MainPar(uc_strat = 3, transformation = true , alpha = 3, benders = true)
name = "instances/small_instances/small_instance_10.dat"
n, oc, sc, rc = read_input_random(name, pars)
V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt = _declare_set(n, pars)
opening_cost, ring_cost, star_cost = oc, rc, sc
if pars.transformation
    offset, oc, rc, sc, backup = _transformation_cost(rc,sc, oc, n, V_tilt, V_certain)
end
lb_distance = _find_lower_bound_backup(n, V_tilt, rc)


master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 1))

@variable(master, x[i in V, j in V; i<j], Bin)
@variable(master, y[i in V], Bin)
@variable(master, lambda_01>=0)
@variable(master, lambda_02>=0)
@variable(master, lambda[i in V]>=0)

if pars.transformation
    @objective(master, Min, offset + sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i] for i in V)+ lambda_01 + lambda_02 + sum(lambda[i] for i in V))
else
    @objective(master, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i] for i in V)+ lambda_01 + lambda_02 + sum(lambda[i] for i in V))
end

@constraint(master, degree_constr[i in V] ,sum(x[minmax(i,j)] for j in V if i!=j)==  2*y[i])
@constraint(master, sum(x[i,j] for (i,j) in E) >= 6)
@constraint(master, y[1] == 1)
if length(V_tilt) >= 1
    if pars.transformation
        @constraint(master, lambda_01 >= sum(y[i]* lb_distance[i] for i in V_tilt))
    else
        @constraint(master, lambda_01 + lambda_02 >= sum(y[i]* lb_distance[i] for i in V_tilt))
    end
end


function main()
    global_upper_bound = 1e18
    function my_callback_benders_cut(cb_data)
        x_hat = Bool.(round.(callback_value.(cb_data, x)))
        y_hat = Bool.(round.(callback_value.(cb_data, y)))
        x_hat = _transform_matrix(x_hat)
        
        status = callback_node_status(cb_data, master)
        all_cycles = find_cycle(x_hat, y_hat)

        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            if length(all_cycles) > 1
                _list_hub = [i for i in 1:n if y_hat[i] == 1]
                # add subtour elimination
                for each_cycle in all_cycles
                    con = @build_constraint(length(each_cycle) - 1/(length(_list_hub)- length(each_cycle))*sum(y[i] for i in _list_hub if i ∉ each_cycle)>= 
                    sum(x[minmax(each_cycle[i], each_cycle[i+1])] for i in eachindex(each_cycle[1:end-1]))+ x[minmax(each_cycle[1], each_cycle[end])])
                    MOI.submit(master, MOI.LazyConstraint(cb_data), con)
                end
            elseif length(all_cycles) == 1
                lambda_01_hat = round(callback_value(cb_data, lambda_01))
                lambda_02_hat = round(callback_value(cb_data, lambda_02))
                lambda_hat = round.(callback_value.(cb_data, lambda))
                
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
                    lower_bound = offset + sum(rc[i,j]*x_hat[i,j] for (i,j) in E)+ sum(oc[i]*y_hat[i] for i in V) + lambda_01_hat + lambda_02_hat + sum(lambda_hat[i] for i in V)
                else
                    lower_bound = sum(ring_cost[i,j]*x_hat[i,j] for (i,j) in E)+ sum(opening_cost[i]*y_hat[i] for i in V)+ lambda_01_hat + lambda_02_hat + sum(lambda_hat[i] for i in V)
                end
                
                gap = (upper_bound - lower_bound)/upper_bound
                # ADD SP0
                if !(lambda_01_hat >= sum((x_hat[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ x_hat[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* beta[(val[1],val[2],val[3])] for val in keys(beta)))
                    cut = @build_constraint(lambda_01 >= sum((x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ x[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* beta[(val[1],val[2],val[3])] for val in keys(beta)))
                    MOI.submit(master, MOI.LazyConstraint(cb_data), cut)
                end

                if !(lambda_02_hat >= sum((x_hat[minmax(val[3],val[1])[1],minmax(val[3],val[1])[2]] + x_hat[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+x_hat[minmax(val[2],val[4])[1],minmax(val[2],val[4])[2]] - 2) * alpha[(val[1],val[2],val[3],val[4])] for val in keys(alpha)))
                    cut = @build_constraint(lambda_02 >= sum((x[minmax(val[3],val[1])[1],minmax(val[3],val[1])[2]] + x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+x[minmax(val[2],val[4])[1],minmax(val[2],val[4])[2]] - 2) * alpha[(val[1],val[2],val[3],val[4])] for val in keys(alpha)))
                    MOI.submit(master, MOI.LazyConstraint(cb_data), cut)
                end
                # Add SPi
                for i in 1:n
                    y_hat[i] == 0 && φ[i]!= 0|| continue
                    if !(lambda_hat[i] >= 3(1-y_hat[i])*φ[i] - sum(y_hat[j]*γ[i,j] for j in V if j!=i))
                        con = @build_constraint(lambda[i] >= 3(1-y[i])*φ[i] - sum(y[j]*γ[i,j] for j in V if j!=i))
                        MOI.submit(master, MOI.LazyConstraint(cb_data), con)     
                    end
                end
            end

                # add cut SP0

        end
    end
    
    # # # # # # # # # # # # 
    # MAIN PROGRAM
    # # # # # # # # # # # # 
    set_attribute(master, MOI.LazyConstraintCallback(), my_callback_benders_cut)
    optimize!(master)
end

main()