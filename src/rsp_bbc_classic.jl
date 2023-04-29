using JuMP, Gurobi
using Combinatorics
using Graphs, GraphsFlows

import GLPK
import MathOptInterface as MOI
import Dates

include("dat.jl")
include("dual_solution.jl")
include("write_output.jl")
include("misc.jl")
include("mutable_structure.jl")
include("user_cut.jl")

# Initialize Sets
V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt = _declare_set(n, 0)
opening_cost, ring_cost, star_cost = oc, rc, sc
pars = MainPar(uc_strat = 2)
benders_enabled = true

add_SP0 = false 
name1 = "cut_bbc_classic"
name2 = "log_file_BBC_classic"
if add_SP0
    name1 = "cut_bbc_classic_SP0y"
    name2 = "log_file_BBC_classic_SP0y"
end
############# MASTER PROBLEM ########################################

master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 1))
# Decision variables
@variable(master, x[i in V, j in V; i<j], Bin)
@variable(master, y[i in V], Bin)
@variable(master, sigma>=0, Int)
@variable(master, lambda_0>=0)
@variable(master, lambda[i in V]>=0)

# Objective Function
@objective(master, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))

# Constraint
@constraint(master, degree_constr[i in V] ,sum(x[minmax(i,j)] for j in V if i!=j)==  2*y[i])
@constraint(master, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
@constraint(master, [(i,j) in T_tilt], sigma >= y[i] + y[j])
@constraint(master, y[1] == 1)

function main_program()

    function my_callback_benders_cut(cb_data)
        
        x_hat_1 = Bool.(round.(callback_value.(cb_data, x)))
        y_hat = Bool.(round.(callback_value.(cb_data, y)))
        x_hat = _transform_matrix(x_hat_1)
        
        status = callback_node_status(cb_data, master)
        _list_hub = [i for i in 1:n if y_hat[i] == 1]
        
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            # Check subtour in a tour, if yes -> add benders cut, if no add subtour elimination
            all_cycles = find_cycle(x_hat, y_hat)
            if length(all_cycles) == 1
                
                lambda_0_hat = round(callback_value(cb_data, lambda_0))
                lambda_hat = round.(callback_value.(cb_data, lambda))
                
                lower_bound = sum(ring_cost[i,j]*x_hat[i,j] for (i,j) in E)+ sum(opening_cost[i]*y_hat[i] for i in V)+ lambda_0_hat + sum(lambda_hat[i] for i in V)
                
                (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat, ring_cost, ring_cost, star_cost)
                
                obj_sp0 = cal_obj_sp0(alpha, beta, x_hat)
                obj_spi = cal_obj_spi(φ, γ, y_hat)
                
                upper_bound = (sum(ring_cost[i,j] * x_hat[i,j] for (i,j) in E)+ sum(opening_cost[i] * y_hat[i] for i in V)) + obj_sp0 + obj_spi
                
                gap = (upper_bound - lower_bound)/upper_bound
                if gap > 1e-10
                    
                    open("$(name1).txt", "a") do io
                        println(io, "=======BENDERS CUTS==============")
                        println(io, "Route: $(transform_route(x_hat))")
                    end

                    # add cut SP0
                    
                    if !(lambda_0_hat >= sum((x_hat_1[minmax(k,i)] + x_hat_1[minmax(i,j)]+ x_hat_1[minmax(j,t)]-2)*alpha[i,j,k,t] for (i,j,k,t) in find_index(alpha))+ sum((x_hat_1[minmax(i,j)]+x_hat_1[minmax(j,k)]-1)*beta[i,j,k] for (i,j,k) in find_index(beta)))
                        con = @build_constraint(lambda_0 >= sum((x[minmax(k,i)]+x[minmax(i,j)]+x[minmax(j,t)]-2)*alpha[i,j,k,t] for (i,j,k,t) in find_index(alpha))+ sum((x[minmax(i,j)]+x[minmax(j,k)]-1)*beta[i,j,k] for (i,j,k) in find_index(beta)))
                        MOI.submit(master, MOI.LazyConstraint(cb_data), con)
                        write_to_log("$(name1).txt", "Cut SP0(x)", con)
                        
                        if add_SP0
                            distance = _find_lower_bound_backup(_list_hub, ring_cost)
                            if !(lambda_0_hat >= sum(y_hat[i]*distance[i] for i in _list_hub))
                                con2 = @build_constraint(lambda_0 >= sum(y[i]*distance[i] for i in _list_hub)) 
                                MOI.submit(master, MOI.LazyConstraint(cb_data), con2)
                                write_to_log("$(name1).txt", "Cut SP0(y[i])", con2)
                            end
                        end
                        
                        
                    end
                
                    # add cut SP_i
                    for i in 1:n
                        y_hat[i] == 0 && φ[i] != 0 || continue
                        
                        if !(lambda_hat[i] >= 3(1-y_hat[i])*φ[i] - sum(y_hat[j]*γ[i,j] for j in V if j!=i))
                            con = @build_constraint(lambda[i] >= 3(1-y[i])*φ[i] - sum(y[j]*γ[i,j] for j in V if j!=i))
                            MOI.submit(master, MOI.LazyConstraint(cb_data), con)                            
                            
                            write_to_log("$(name1).txt", "Cut SPi", con)

                        end
                    end
                end
            else
                
                open("$(name1).txt", "a") do io
                    println(io, "=======SUBTOUR==========")
                    for each_cycle in all_cycles
                        println(io, "Subtour: $(each_cycle)")
                    end
                end
                # add subtour elimination
                for each_cycle in all_cycles
                    con = @build_constraint(length(each_cycle) - 1/(length(_list_hub)- length(each_cycle))*sum(y[i] for i in _list_hub if i ∉ each_cycle)>= 
                    sum(x[minmax(each_cycle[i], each_cycle[i+1])] for i in eachindex(each_cycle[1:end-1]))+ x[minmax(each_cycle[1], each_cycle[end])])
                    MOI.submit(master, MOI.LazyConstraint(cb_data), con)
    
                    write_to_log("$(name1).txt", "Subtour = ", con)
                end
                
            end
        end
    end

    open("$(name1).txt","w") do io
        println(io, "BBC_Classic")
    end

    set_attribute(master, MOI.LazyConstraintCallback(), my_callback_benders_cut)
    set_attribute(master, MOI.UserCutCallback(), call_back_user_cuts_benders)
    optimize!(master)
    
    x_hat = value.(x)
    x_hat = _transform_matrix(x_hat)
    lambda_0_hat = value(lambda_0)
    lambda_hat = value.(lambda)
    
    open("$(name2).txt","w") do io
        println(io, "Obj = $(objective_value(master))")
        println(io, "Running time = $(solve_time(master))")
        println(io, "Gap = $(relative_gap(master))")
        println(io, "Route =  $(transform_route(x_hat))")
        println(io, "lambda SP0 =  $(lambda_0_hat)")
        println(io, "lambda SPi =  $(sum(lambda_hat))")
        println(io, master)
        println(io, "===================")
        println(io, solution_summary(master))
    end
    
end

main_program()