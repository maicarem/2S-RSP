using JuMP, Gurobi
using Combinatorics

include("dat.jl")

############# MASTER PROBLEM ########################################
master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
# Decision variables
@variable(master, x[i in V, j in V], Bin)
@variable(master, y[i in V,j in V], Bin)
@variable(master, sigma>=0, Int)
@variable(master, lambda_0>=0)
@variable(master, lambda[i in V]>=0)

# Objective Function
@objective(master, Min, sum(rc[i,j]*x[i,j] for (i,j) in E)+ sum(oc[i]*y[i,i] for i in V)+ lambda_0 + sum(lambda[i] for i in V))

# Constraint
@constraint(master, degree_constr[i in V] ,sum(x[i,j] for j in V if i<j) + sum(x[j,i] for j in V if i>j)==  2*y[i,i])

for i in 2:ceil(Int,length(V)/2)
    for S in combinations(V, i)
        @constraint(master, length(S) - 1/length(V)*sum(y[i,i] for i in S)>= sum(x[i,j] for i in S for j in S if i<j))
    end
end

@constraint(master, [(i,j) in E], x[i,j] == x[j,i])
@constraint(master, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
@constraint(master, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
@constraint(master, y[1,1] == 1)


############# SUBPROBLEM SP_0 ########################################
sub_0 = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

@variable(sub_0, x_prime[i in V, j in V; i!=j] >=0)
@objective(sub_0, Min, sum(rc[i,j]*x_prime[i,j] for (i,j) in E))

@constraint(sub_0, th_constr[(i,j,k,t) in J_tilt], x_prime[k,t]>= 0)
@constraint(sub_0, oh_constr[(i,j,k) in K_tilt], x_prime[i,k]>= 0)

############# SUBPROBLEM SP_i ########################################

sub_probs= Dict()

function solve_SP0(x_hat)
    for (i,j,k,t) in J_tilt
        set_normalized_rhs(th_constr[(i,j,k,t)], x_hat[k,i]+ x_hat[i,j]-2)
    end

    for (i,j,k) in K_tilt
        set_normalized_rhs(oh_constr[(i,j,k)], x_hat[i,j] + x_hat[j,k] -1)
    end

    optimize!(sub_0)
    return (obj = objective_value(sub_0), y = value.(x_prime), α = dual.(th_constr), β = dual.(oh_constr))
end

# Find Omega i (indices where values <= second smallest value)
function find_Omega(y_hat, sc, H_tilt,  rank0)
    smallest = Dict()
    idx_smallest = Dict()
    for i in V
        if y_hat[i,i] == 0
            if length(H_tilt) == 1
                smallest[i] = sort(sc[i,H_tilt])[1]
            else
                smallest[i] = sort(sc[i,H_tilt])[rank0]
            end
            idx_smallest[i] = findall(x -> x <= smallest[i], sc[i,H_tilt]) #note
        end
    end
    return smallest, idx_smallest # (value of s_imi in uncertain hubs, indices of these uncertain hubs)
end

function solve_SPi(y_hat)
    H_tilt = []
    for i in V
        if y_hat[i,i] == 1 && i in V_tilt
            append!(H_tilt, i)
        end
    end

    Omega = find_Omega(y_hat, sc, H_tilt,  2)[2]
    
    φ = Dict()
    γ = Dict()
    obj = Dict()
    y_val = Dict()
    
    for i in V
        if y_hat[i,i] == 0
            sub_probs[i] = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
            
            # Declare variables
            @variable(sub_probs[i], y_prime[i,j in V; j!=i]>=0)
            # Objective function for SP i
            @objective(sub_probs[i], Min, sum(sc[i,j]*y_prime[i,j] for j in V if i!=j))
            
            # Constraint
            sub_constr_i0 = @constraint(sub_probs[i], 3*sum(y_prime[i,j] for j in V_certain if i!=j) + sum(y_prime[i,j] for j in V_tilt if i!=j) == 3*(1-y_hat[i,i]))
            @constraint(sub_probs[i], sub_constr_i1[j in V; j!=i], -y_prime[i,j] >= -y_hat[j,j])
            @constraint(sub_probs[i], sub_constr_i2[j in Omega[i]; j!=i], -2*y_prime[i,j] + sum(y_prime[i,k] for k in H_tilt if k!=j) >=0)
            @info "Adding constraint $(sub_constr_i2)"
            optimize!(sub_probs[i])
            obj[i] = objective_value(sub_probs[i])
            y_val[i] = value.(y_prime)

            # Obtain dual values
            φ[i] = dual.(sub_constr_i0)
            γ[i] = dual.(sub_constr_i1)
        end
    end
    return obj, y_val, Omega, φ, γ
end

function add_cut_SP0(ret)
    cut = @constraint(master, lambda_0 >= sum((x[k,i]+2*x[i,j]+x[j,t]-3)*ret[3][(i,j,k,t)] 
                        for (i,j,k,t) in J_tilt)+ sum((x[i,j]+x[j,k]-1)*ret[4][(i,j,k)] for (i,j,k) in K_tilt))
    @info "Adding the cut $(cut)"
end

function add_cut_SPi(ret, indice)
    cut = @constraint(master, lambda[indice] >= 3*(1-y[indice,indice])*ret[4][indice] - sum(y[j,j]*ret[5][indice][j] for j in V if j!=indice))
    @info "Adding the cut $(cut)"
end

################ MAIN PROGRAM ######################
function main_program()
    iter0 = 0 
    cut_found = 1
    while cut_found == 1
        println("===============Iteration ", iter0, "===============")
        iter0 +=1
        cut_found = 0 
        optimize!(master)
        @show value.(x);
        @show value.(y);

        lower_bound = objective_value(master)
        
        println("Current Master Problem = ", objective_value(master))
        @assert termination_status(master) == OPTIMAL
        x_hat, y_hat = value.(x), value.(y)
        sp0_res = solve_SP0(x_hat)
        spi_res = solve_SPi(y_hat)
        
        if spi_res[1] == Dict{Any, Any}()
            upper_bound = sp0_res[1]
        else
            # upper_bound = sp0_res[1]
            upper_bound = sp0_res[1] + sum(i[2] for i in spi_res[1])
        end
        
        gap = (upper_bound - lower_bound) / upper_bound
        
        if gap < 1e-18
            println("Terminating with the optimal solution")
            println("Objective value ", upper_bound + objective_value(master))
            break
        end

        add_cut_SP0(sp0_res)
        
        for i in V
            if y_hat[i,i] == 0
                add_cut_SPi(spi_res, i)
            end
        end
        cut_found = 1
        
    end
end

# main_program()