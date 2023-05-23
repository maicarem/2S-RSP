using JuMP, Gurobi
using Combinatorics
using Graphs, GraphsFlows


include("dat.jl")
include("dual_solution.jl")
include("misc.jl")
include("mutable_structure.jl")
include("user_cut.jl")


pars = MainPar(uc_strat = 3, transformation = true, alpha = 3, benders = false)
name = "instances/small_instances/small_instance_10.dat"
n, oc, sc, rc = read_input_random(name, pars)
V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt = _declare_set(n, pars)
opening_cost, ring_cost, star_cost = oc, rc, sc
if pars.transformation
    offset, oc, rc, sc, backup = _transformation_cost(rc,sc, oc, n, V_tilt, V_certain)
end
lb_distance = _find_lower_bound_backup(n, V_tilt, rc)
############# MASTER PROBLEM ########################################

master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
# Decision variables
@variable(master, x[i in V, j in V; i<j], Bin)
@variable(master, y[i in V], Bin)
@variable(master, sigma>=0, Int)
@variable(master, lambda_0>=0)
@variable(master, lambda[i in V]>=0)

# Objective Function
if pars.transformation
    @objective(master, Min, offset + sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))
else
    @objective(master, Min, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))
end

# Constraint
@constraint(master, degree_constr[i in V] ,sum(x[minmax(i,j)] for j in V if i!=j)==  2*y[i])
@constraint(master, sum(x[i,j] for (i,j) in E) >= 6)
@constraint(master, y[1] == 1)
if length(V_tilt) >= 1
    @constraint(master, lambda_0 >= sum(y[i]* lb_distance[i] for i in V_tilt))
end

function main_program()
    global_upper_bound = 1e18
    for iter0 in 1:200
        println("===============Iteration ", iter0, "===============")
        
        function my_callback_subtour(cb_data)
            
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

        set_attribute(master, MOI.LazyConstraintCallback(), my_callback_subtour)
        # set_attribute(master, MOI.UserCutCallback(), call_back_user_cuts_benders)
        
        optimize!(master)
        lower_bound = objective_value(master)
        println("Objective value at iteration $(iter0) is $(lower_bound)")
        x_hat_1, y_hat = Bool.(round.(value.(x))), Bool.(round.(value.(y)))
        x_hat = _transform_matrix(x_hat_1)
        lambda_0_hat, lambda_hat = value(lambda_0), round.(value.(lambda))
        
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

        _add_cut_SP0(master, alpha,beta, lambda_0_hat, x_hat)
        for i in 1:n
            y_hat[i] == 0 && φ[i]!= 0|| continue
            _add_cut_SPi(master, φ, γ, i, lambda_hat, y_hat)
        end
    end
end

main_program()

