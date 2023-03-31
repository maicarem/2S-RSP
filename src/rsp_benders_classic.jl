using JuMP, Gurobi
using Combinatorics

include("dat.jl")
include("dual_solution.jl")
include("misc.jl")

# Add cut constraint
function _add_cut_SP0(_alpha,_beta,_lambda_0, _x_hat)
    if !(_lambda_0 >= sum((_x_hat[minmax(k,i)]+2*_x_hat[minmax(i,j)]+_x_hat[minmax(j,t)]-3)*_alpha[i,j,k,t] for (i,j,k,t) in J_tilt)+ sum((_x_hat[minmax(i,j)]+_x_hat[minmax(j,k)]-1)*_beta[i,j,k] for (i,j,k) in K_tilt))
        cut = @constraint(master, lambda_0 >= sum((x[minmax(k,i)]+2*x[minmax(i,j)]+x[minmax(j,t)]-3)*_alpha[i,j,k,t]
                            for (i,j,k,t) in J_tilt)+ sum((x[minmax(i,j)]+x[minmax(j,k)]-1)*_beta[i,j,k] for (i,j,k) in K_tilt))
        @info "Adding the cut $(cut)"
    end
end

function _add_cut_SPi(_varphi, _gamma, indice, _lambda, _y_hat)
    if !(_lambda[indice] >= 3(1-_y_hat[indice,indice])*_varphi[indice] - sum(_y_hat[j,j]*_gamma[indice,j] for j in V if j!=indice))
        cut = @constraint(master, lambda[indice] >= 3(1-y[indice,indice])*_varphi[indice] - sum(y[j,j]*_gamma[indice,j] for j in V if j!=indice))
        @info "Adding the cut $(cut)"
    end
end

############# MASTER PROBLEM ########################################

master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
# Decision variables
@variable(master, x[i in V, j in V; i<j], Bin)
@variable(master, y[i in V,j in V], Bin)
@variable(master, sigma>=0, Int)
@variable(master, lambda_0>=0)
@variable(master, lambda[i in V]>=0)

# Objective Function
@objective(master, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i,i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))

# Constraint
@constraint(master, degree_constr[i in V] ,sum(x[minmax(i,j)] for j in V if i!=j)==  2*y[i,i])
@constraint(master, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
@constraint(master, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
@constraint(master, y[1,1] == 1)

# @constraint(master, x[1,2] == 1)
# @constraint(master, x[1,9] == 1)
# @constraint(master, x[7,9] == 1)
# @constraint(master, x[2,3] == 1)
# @constraint(master, x[3,7] == 1)

function main_program()
    for iter0 in 1:1000
        println("===============Iteration ", iter0, "===============")
        
        function my_callback_benders_cut(cb_data)
            
            x_hat = Bool.(round.(callback_value.(cb_data, x)))
            y_hat = Bool.(round.(callback_value.(cb_data, y)))
            x_hat = _transform_matrix(x_hat)
            status = callback_node_status(cb_data, master)
            
            if status == MOI.CALLBACK_NODE_STATUS_INTEGER
                # Check subtour in a tour
                if contain_subtour(x_hat, y_hat) == 0
                    _list_hub = [i for i in 1:n if y_hat[i,i] == 1]
                    for i in 2:ceil(Int,length(_list_hub)/2)
                        for S in combinations(_list_hub, i)
                            con = @build_constraint(length(S) - 1/length(_list_hub)*sum(y[i,i] for i in S)>= sum(x[minmax(i,j)] for i in S for j in S if i<j))
                            MOI.submit(master, MOI.LazyConstraint(cb_data), con)
                        end
                    end
                end
            end
        end

        set_attribute(master, MOI.LazyConstraintCallback(), my_callback_benders_cut)
        optimize!(master)
        lower_bound = objective_value(master)
        println("Objective value at iteration $(iter0) is $(lower_bound)")
        x_hat_1, y_hat = Bool.(round.(value.(x))), Bool.(round.(value.(y)))
        x_hat = _transform_matrix(x_hat_1)
        lambda_0_hat, lambda_hat = value(lambda_0), round.(value.(lambda))

        (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat, ring_cost, star_cost)
        # Objective value, and add cut
        
        obj_sp0 = cal_obj_sp0(alpha, beta, x_hat)
        obj_spi = cal_obj_spi(φ, γ, y_hat)
        
        

        upper_bound =  sum(ring_cost[i,j]*x_hat[i,j] for (i,j) in E)+ sum(opening_cost[i]*y_hat[i,i] for i in V) + obj_sp0 + obj_spi

        open("result/bender/debug_$(iter0).txt","w") do io
            println(io, "Lower bound: $(lower_bound)")
            println(io, "Route: $(transform_route(x_hat))")
            
            println(io, "For master prob backup edges: $(lambda_0_hat)")
            println(io, "For master prob star edges: $(sum(lambda_hat))")
            
            println(io, "For sp backup edges: $(obj_sp0)")
            println(io, "For sp star edges: $(obj_spi)")
            
            println(io, master)
            
        end

        gap = (upper_bound - lower_bound)/upper_bound
        
        println("Upper bound is $(upper_bound)")
        
        
        if  gap < 1e-10
            println("This is optimal with objective $(lower_bound)")
            break
        end

        _add_cut_SP0(alpha,beta, lambda_0_hat, x_hat_1)
        for i in 1:n
            y_hat[i,i] == 0 && φ[i]!= 0|| continue
            _add_cut_SPi(φ, γ, i, lambda_hat, y_hat)
        end
    end
end

main_program()