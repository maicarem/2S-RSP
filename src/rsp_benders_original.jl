using JuMP, Gurobi
using Combinatorics

include("dat.jl")
include("dual_solution.jl")

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

master = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
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

for i in 2:ceil(Int,length(V)/2)
    for S in combinations(V, i)
        @constraint(master, length(S) - 1/length(V)*sum(y[i,i] for i in S)>= sum(x[i,j] for i in S for j in S if i<j))
    end
end

@constraint(master, [(i,j) in E], x[i,j] == x[j,i])
@constraint(master, sum(x[i,j] for (i,j) in E) >= 3+ sigma)
@constraint(master, [(i,j) in T_tilt], sigma >= y[i,i] + y[j,j])
@constraint(master, y[1,1] == 1)

function main_program()
    for iter0 in 1:1000
        println("===============Iteration ", iter0, "===============")
        optimize!(master)
        lower_bound = objective_value(master)
        println("Objective value at iteration $(iter0) is $(lower_bound)")
        x_hat, y_hat = value.(x), value.(y)
        (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat)
        

        # Objective value, and add cut
        
        obj_sp0 = cal_obj_sp0(alpha, beta, x_hat)
        obj_spi = cal_obj_spi(φ, γ, y_hat)
        upper_bound = (offset + sum(rc[i,j]*x_hat[i,j] for (i,j) in E)+ sum(oc[i]*y_hat[i,i] for i in V)) + obj_sp0 + obj_spi
        gap = (upper_bound - lower_bound)/upper_bound
        
        println("Upper bound is $(upper_bound)")
        
        
        if  gap < 1e-10
            println("This is optimal with objective $(lower_bound)")
            break
        end

        _add_cut_SP0(alpha,beta)
        for i in 1:n
            y_hat[i,i] == 0 || continue
            _add_cut_SPi(φ, γ, i)
        end
    end
end

main_program()