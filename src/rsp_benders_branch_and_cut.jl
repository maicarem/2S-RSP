using JuMP, Gurobi
using Combinatorics
import GLPK
import MathOptInterface as MOI

include("dat.jl")
include("dual_solution.jl")
include("write_output.jl")

############# MASTER PROBLEM ########################################

master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 1))
# Decision variables
@variable(master, x[i in V, j in V; i<j], Bin)
@variable(master, y[i in V,j in V], Bin)
@variable(master, sigma>=0, Int)
@variable(master, lambda_0>=0)
@variable(master, lambda[i in V]>=0)

# Objective Function
@objective(master, Min, offset + sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i,i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))

# Constraint
@constraint(master, degree_constr[i in V] ,sum(x[minmax(i,j)] for j in V if i!=j)==  2*y[i,i])
@constraint(master, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
@constraint(master, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
@constraint(master, y[1,1] == 1)


function main_program()

    function my_callback_benders_cut(cb_data)
        x_hat = Bool.(round.(callback_value.(cb_data, x)))
        y_hat = Bool.(round.(callback_value.(cb_data, y)))
        x_hat = _transform_matrix(x_hat)
        status = callback_node_status(cb_data, master)
        
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            # Check subtour in a tour, if yes -> add benders cut, if no add subtour elimination
            if contain_subtour(x_hat, y_hat)

                lambda_0_hat = callback_value(cb_data, lambda_0)
                lambda_hat = callback_value.(cb_data, lambda)
                
                lower_bound = offset + sum(rc[i,j]*x_hat[i,j] for (i,j) in E)+ sum(oc[i]*y_hat[i,i] for i in V)+ lambda_0_hat + sum(lambda_hat[i] for i in V)
                (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat)

                obj_sp0 = cal_obj_sp0(alpha, beta, x_hat)
                obj_spi = cal_obj_spi(φ, γ, y_hat)
                
                # write results
                open("debug.txt","w") do io
                    for i in 1:n
                        if y_hat[i,i] == 0 && φ[i]!= 0
                            println(io, "varphi[$(i)] = $(φ[i])")
                        end
                        for j in 1:n
                            if x_hat[i,j] ==1
                                println(io, "x[$(i), $(j)] = 1, cost = $(rc[i,j])")
                            end
                            for k in 1:n
                                if beta[i,j,k] != 0
                                    println(io, "beta[$(i),$(j),$(k)] = $(beta[i,j,k]), cost = $(backup[i,k])")
                                end
                                for t in 1:n
                                    if alpha[i,j,k,t] != 0
                                        println(io, "alpha[$(i),$(j),$(k), $(t)] = $(alpha[i,j,k,t]), cost = $(backup[k,t])")
                                    end
                                end
                            end
                        end
                    end  
                    println(io,"Result for SP0: $(obj_sp0)")   
                    println(io,"Result for SPi: $(obj_spi)")   
                end

                upper_bound = (offset + sum(rc[i,j] * x_hat[i,j] for (i,j) in E)+ sum(oc[i] * y_hat[i,i] for i in V)) + obj_sp0 + obj_spi
                
                gap = (upper_bound - lower_bound)/upper_bound
                
                if gap > 1e-10
                    
                    # add cut SP0
                    con = @build_constraint(lambda_0 >= sum((x[minmax(k,i)]+2*x[minmax(i,j)]+x[minmax(j,t)]-3)*alpha[i,j,k,t] for (i,j,k,t) in J_tilt)+ sum((x[minmax(i,j)]+x[minmax(j,k)]-1)*beta[i,j,k] for (i,j,k) in K_tilt))
                    MOI.submit(master, MOI.LazyConstraint(cb_data), con)
                    
                    # add cut SP_i
                    for i in 1:n
                        y_hat[i,i] == 0 || continue
                        con = @build_constraint(lambda[i] >= 3(1-y[i,i])*φ[i] - sum(y[j,j]*γ[i,j] for j in V if j!=i))
                        MOI.submit(master, MOI.LazyConstraint(cb_data), con)
                    end
                end
            else
                # add subtour elimination
                _list_hub = [i for i in 1:n if y_hat[i,i] == 1]
                for i in 2:ceil(Int,length(_list_hub)/2)
                    for S in combinations(_list_hub, i)
                        con = @build_constraint(length(S) - 1/length(V)*sum(y[i,i] for i in S)>= sum(x[i,j] for i in S for j in S if i<j))
                        MOI.submit(master, MOI.LazyConstraint(cb_data), con)
                    end
                end
            end
        end
    end
    
    set_attribute(master, MOI.LazyConstraintCallback(), my_callback_benders_cut)
    optimize!(master)
    _write_gurobi_log("instance", master, "Branch and cut")
    # @show solution_summary(master)
end

main_program()