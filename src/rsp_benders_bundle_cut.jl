using JuMP, Gurobi
using Combinatorics

include("dat.jl")
include("dual_solution.jl")
include("misc.jl")

# Add cut constraint
V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt = _declare_set(n, 1)
opening_cost, ring_cost, star_cost = oc, rc, sc
add_SP0 = false

############# MASTER PROBLEM ########################################

master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
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
    for iter0 in 1:2
        println("===============Iteration ", iter0, "===============")
        
        function my_callback_benders_cut(cb_data)
            
            x_hat = Bool.(round.(callback_value.(cb_data, x)))
            y_hat = Bool.(round.(callback_value.(cb_data, y)))
            x_hat = _transform_matrix(x_hat)
            status = callback_node_status(cb_data, master)
            all_cycles = find_cycle(x_hat, y_hat)

            if status == MOI.CALLBACK_NODE_STATUS_INTEGER
                # Check subtour in a tour
                if length(all_cycles) > 1
                    _list_hub = [i for i in 1:n if y_hat[i] == 1]
                    # add subtour elimination
                    for each_cycle in all_cycles
                        con = @build_constraint(length(each_cycle) - 1/(length(_list_hub)- length(each_cycle))*sum(y[i] for i in _list_hub if i ∉ each_cycle)>= 
                        sum(x[minmax(each_cycle[i], each_cycle[i+1])] for i in eachindex(each_cycle[1:end-1]))+ x[minmax(each_cycle[1], each_cycle[end])])
                        MOI.submit(master, MOI.LazyConstraint(cb_data), con)
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
        @show x_hat_1
        (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat, ring_cost, ring_cost, star_cost)
        # Objective value, and add cut
        
        obj_sp0 = cal_obj_sp0(alpha, beta, x_hat)
        obj_spi = cal_obj_spi(φ, γ, y_hat)
        
        upper_bound =  sum(ring_cost[i,j]*x_hat[i,j] for (i,j) in E)+ sum(opening_cost[i]*y_hat[i] for i in V) + obj_sp0 + obj_spi

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
        
        _add_cut_SP0(master, alpha, beta, lambda_0_hat, x_hat_1)

        for i in 1:n
            y_hat[i] == 0 && φ[i]!= 0|| continue
            _add_cut_SPi(master, φ, γ, i, lambda_hat, y_hat)
        end

        subproblem0 = Model(optimizer_with_attributes(Gurobi.Optimizer))        
        @variable(subproblem0, α[(i,j,k,t) in J_tilt] >= 0, Int)
        @variable(subproblem0, β[(i,j,k) in K_tilt] >= 0, Int)
        
        @show find_index(alpha)

        # @objective(subproblem0, Max, sum(α[(i,j,k,t)]*(x_hat[minmax(k,i)]+ x_hat[minmax(i,j)] + x_hat[minmax(j,t) - 2]) for (i,j,k,t) in find_index(alpha)) + sum(β[(i,j,k)]*(x_hat[minmax(i,j)]+x_hat[minmax(j,k)]-1) for (i,j,k) in find_index(beta)))
        @objective(subproblem0, Max, sum(α[(i,j,k,t)] for (i,j,k,t) in find_index(alpha)) + sum(β[(i,j,k)] for (i,j,k) in find_index(beta)))
        
        @constraint(subproblem0, [(i,j) in E], sum(α[(m,n,i,j)] for m in V_tilt for n in V_tilt if m!=n && m!= i && m!= j && n!= i && n!= j) + sum(β[(i,k,j)] for k in V_tilt if k!= i && k!= j) <= ring_cost[i,j])
        
        @constraint(subproblem0, -α[(3,4,1,2)] <= -calculate_eta(alpha))
        @constraint(subproblem0, α[(3,4,1,2)] <= calculate_eta(alpha))
        @constraint(subproblem0, α[(5,6,1,2)] <= calculate_eta(alpha)/0.1)
        
        optimize!(subproblem0)
        α_hat = value.(α)
        for (index, value) in enumerate(α_hat)
            value != 0 || continue
            println("Index = $(index), Value = $(value)")
        end

        @show obj_sp0


    end
end

main_program()