using JuMP, Gurobi
using Combinatorics
import GLPK
import MathOptInterface as MOI

include("dat.jl")
include("dual_solution.jl")
include("ga_warm_start.jl")

# Add cut constraint
function _add_cut_SP0(_alpha,_beta)
    cut = @constraint(master, lambda_0 >= sum((x[k,i]+2*x[i,j]+x[j,t]-3)*_alpha[i,j,k,t]
                        for (i,j,k,t) in J_tilt)+ sum((x[i,j]+x[j,k]-1)*_beta[i,j,k] for (i,j,k) in K_tilt))
    @info "Adding the cut $(cut)"
end

function _add_cut_SPi(_varphi, _gamma, indice)
    cut = @constraint(master, lambda[indice] >= 3(1-y[indice,indice])*_varphi[indice] - sum(y[j,j]*_gamma[indice,j] for j in V if j!=indice))
    @info "Adding the cut $(cut)"
end

function cal_obj_sp0(alpha, beta, x_hat)
    return sum((x_hat[k,i] + x_hat[i,j]+x_hat[j,t] - 2) * alpha[i,j,k,t] for (i,j,k,t) in J_tilt) + sum((x_hat[i,j]+ x_hat[j,k] - 1)* beta[i,j,k] for (i,j,k) in K_tilt)
end

function cal_obj_spi(varphi, gamma, y_hat)
    obj = zeros(Float64, n)
    for i in 1:n
        y_hat[i,i] == 0 || continue
        obj[i] = 3(1-y_hat[i,i])*varphi[i] + sum([y_hat[j,j]*gamma[i,j] for j in V if j!=i])
    end
    return sum(obj)
end

############# MASTER PROBLEM ########################################

master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 1, "SubMIPNodes" => 20))
# Decision variables
@variable(master, x[i in V, j in V], Bin)
@variable(master, y[i in V,j in V], Bin)
@variable(master, sigma>=0, Int)
@variable(master, lambda_0>=0)
@variable(master, lambda[i in V]>=0)

# Objective Function
@objective(master, Min, offset + sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i,i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))

# Constraint
@constraint(master, degree_constr[i in V] ,sum(x[i,j] for j in V if i<j) + sum(x[j,i] for j in V if i>j)==  2*y[i,i])
@constraint(master, [(i,j) in E], x[i,j] == x[j,i])
@constraint(master, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
@constraint(master, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
@constraint(master, y[1,1] == 1)


function main_program()

    function my_callback_benders_cut(cb_data)
        x_hat = Bool.(round.(callback_value.(cb_data, x)))
        y_hat = Bool.(round.(callback_value.(cb_data, y)))
        status = callback_node_status(cb_data, master)
        
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            # Check subtour in a tour
            if contain_subtour(x_hat, y_hat)
                pass
            else
                _list_hub = [i for i in 1:n if y_hat[i,i] == 1]
                for i in 2:ceil(Int,length(_list_hub)/2)
                    for S in combinations(_list_hub, i)
                        con = @build_constraint(length(S) - 1/length(V)*sum(y[i,i] for i in S)>= sum(x[i,j] for i in S for j in S if i<j))
                        MOI.submit(master, MOI.LazyConstraint(cb_data), con)
                    end
                end
            end

            # Add cut from subproblems
            
            lambda_0_hat = callback_value(cb_data, lambda_0)
            lambda_hat = callback_value.(cb_data, lambda)
            lower_bound = offset + sum(rc[i,j]*x_hat[i,j] for (i,j) in E)+ sum(oc[i]*y_hat[i,i] for i in V)+ lambda_0_hat + sum(lambda_hat[i] for i in V)
            
            (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat)
        
            obj_sp0 = cal_obj_sp0(alpha, beta, x_hat)
            obj_spi = cal_obj_spi(φ, γ, y_hat)

            upper_bound = (offset + sum(rc[i,j] * x_hat[i,j] for (i,j) in E)+ sum(oc[i] * y_hat[i,i] for i in V)) + obj_sp0 + obj_spi
            gap = (upper_bound - lower_bound)/upper_bound
            
            if gap > 1e-10
                
                # add cut SP0
                con = @build_constraint(lambda_0 >= sum((x[k,i]+2*x[i,j]+x[j,t]-3)*alpha[i,j,k,t] for (i,j,k,t) in J_tilt)+ sum((x[i,j]+x[j,k]-1)*beta[i,j,k] for (i,j,k) in K_tilt))
                MOI.submit(master, MOI.LazyConstraint(cb_data), con)
                
                # add cut SP_i
                for i in 1:n
                    y_hat[i,i] == 0 || continue
                    con = @build_constraint(lambda[i] >= 3(1-y[i,i])*φ[i] - sum(y[j,j]*γ[i,j] for j in V if j!=i))
                    MOI.submit(master, MOI.LazyConstraint(cb_data), con)
                end
            end
        end
    end

    
    # ring_list = genetic_algorithm(100, 3)
    # length_ring_list = length(ring_list)
    
    # x_ga = zeros(n,n)
    # y_ga = zeros(n,n)
    

    # for i in 1:length_ring_list-1
    #     x_ga[ring_list[i],ring_list[i+1]] = 1
    #     x_ga[ring_list[i+1],ring_list[i]] = 1
    #     y_ga[ring_list[i], ring_list[i]] = 1
    # end

    # @show ring_list
    # @show y_ga
    # @show x_ga

    # (beta, alpha), (φ, γ) = dual_solution(x_ga, y_ga)
    # @constraint(master, lambda_0 >= sum((x[k,i]+2*x[i,j]+x[j,t]-3)*alpha[i,j,k,t] for (i,j,k,t) in J_tilt)+ sum((x[i,j]+x[j,k]-1)*beta[i,j,k] for (i,j,k) in K_tilt))
    # for i in 1:n
    #     y_ga[i,i] == 0 || continue
    #     @constraint(master, lambda[i] >= 3(1-y[i,i])*φ[i] - sum(y[j,j]*γ[i,j] for j in V if j!=i))
    # end

    # for i in 1:n
    #     set_start_value(y[i,i], y_ga[i,i])
    #     for j in 1:n
    #         set_start_value(x[i,j], x_ga[i,j])

    #     end
    # end
    
    set_attribute(master, MOI.LazyConstraintCallback(), my_callback_benders_cut)
    optimize!(master)
    @show value.(x)
end

main_program()