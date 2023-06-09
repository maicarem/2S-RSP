using LightGraphs

# Find list of indices j having s_ij<= s_imi2
# unhub0: [(ind[i], cost[i]) for i in 1:3]
function find_omega_i(unhub0, hub_list, indice, star_cost0)
    _cost = unhub0[2][2]
    if _cost == 1e18
        _cost = unhub0[1][2]
        if _cost == 1e18
            _cost = 0
        end
    end
    return [i for i in 1:length(star_cost0[indice,:]) if star_cost0[indice,i] <= _cost && i ∈ V_tilt && i ∈ hub_list]
end

# Find first, second and third smallest distance of s_ij, j is uncertain and a hub:: return [(hub1, distance1),...]
# Input: hub_list, V_tilt
function uncertain_smallest(hub_list, indice, star_cost0)
    sorted0 = sort([star_cost0[indice, i] for i in 1:n if i in V_tilt && i in hub_list && i != indice])
    ind, cost = zeros(Int, 3), ones(3)*1e18
    if length(sorted0) != 0
        for i in 1: min(3, length(sorted0))
            ind[i] = [j for j in 1:n if star_cost0[indice, j] == sorted0[i] && j in hub_list && j in V_tilt && j ∉ ind && j != indice][1]
            cost[i] = star_cost0[indice, ind[i]]
        end
    end
    return [(ind[i], cost[i]) for i in 1:3]
end

# Find indice mi^star and smallest cost 
function certain_smallest(hub_list, indice, star_cost0)
    sorted0 = sort([star_cost0[indice, i] for i in 1:n if i ∉ V_tilt && i in hub_list && i != indice])
    if length(sorted0)!=0
        ind = Int([j for j in 1:n if star_cost0[indice, j] == sorted0[1] && j in hub_list && j ∉ V_tilt && j != indice][1])
        cost = star_cost0[indice, ind]
    else
        ind, cost = 0, 1e18
    end
    return (ind, cost)
end

function _uncertain_dict(hub_list, star_cost0)
    _smallest_uncertain_dict = Dict()
    _Omega_dict = Dict()
    for i in 1:n
        if i ∉ hub_list
            _smallest_uncertain_dict[i] = uncertain_smallest(hub_list, i, star_cost0)
            _Omega_dict[i] = find_omega_i(_smallest_uncertain_dict[i], hub_list, i, star_cost0)
        end
    end
    return _smallest_uncertain_dict, _Omega_dict
end

function _list_hub(_y_hat)
    hub = []
    hub_tilt = []
    for i in 1:n
        _y_hat[i] == 1 || continue
        push!(hub, i)
        i ∈ V_tilt || continue
        push!(hub_tilt, i)
    end   
   return hub, hub_tilt 
end



function cal_obj_spi(varphi, gamma, y_hat)
    obj = 0
    for i in 1:n
        y_hat[i] == 0 || continue
        obj+= 3*varphi[i]
    end
    return obj
    # return 3*sum(varphi[i] for i in 1:n if y_hat[i,i] == 0)
end

# function minmax(i,j)
#     return min(i,j), max(i,j)
# end

function _transform_matrix(x)
    x_hat = zeros(Bool, n,n)
    for (i,j) in E
        x_hat[i,j] = Bool(round(x[i,j]))
        x_hat[j,i] = Bool(round(x[i,j]))
    end
    return x_hat
end

function transform_route(x_hat)
    path = [1]
    num0 = size(x_hat)[1]
    while length(path) < 2 || path[1]!=path[end]
        if length(path) == sum(x_hat[i,j] for i in 1:num0 for j in 1:num0 if i<j)
            break
        else
            for vertex in 2:num0
                if x_hat[path[end], vertex] == 1 && vertex ∉ path[2:end]
                    push!(path, vertex)
                    break
                end
            end
        end
    end
    return push!(path,1)
end

function dfs(node, adjacent_list, visited, cycle)
    visited[node] = 1
    for next in adjacent_list[node]
        !(visited[next]) || continue
        append!(cycle, next)
        dfs(next, adjacent_list, visited, cycle)
    end
end

function find_cycle(x_hat, y_hat)
    
    adjacent_list = Dict()
    num0 = size(x_hat)[1]
    visited = zeros(Bool, num0)
    hub_list = []
    all_cycles = []

    for i in 1: num0
        adjacent_list[i] = []
        
        if y_hat[i]
            append!(hub_list, i)
        else
            visited[i] = 1
        end

        for j in 1: num0
            x_hat[i,j] || continue
            append!(adjacent_list[i], j)
        end
    
    end
    
    for i in hub_list
        !(visited[i]) || continue
        cycle = [i]
        dfs(i, adjacent_list, visited, cycle)
        push!(all_cycles, cycle)
    end
    return all_cycles

end

function _find_lower_bound_backup(n, V_tilt, edge_matrix)
    distance = zeros(Float64, n)
    for i in 1:n
        i ∈ V_tilt || continue
        distance[i] = minimum(edge_matrix[[j for j in 1:n if j!=i], [j for j in 1:n if j!=i]])
    end
    return distance
end

function find_index(alpha)
    return Tuple.(findall(!iszero ,alpha))
end

function write_to_log(name, text, con)
    open(name, "a") do io
        println(io, "$(text): $(con)")
    end
end

# function _add_cut_SP0(master, _alpha,_beta,_lambda_0, _x_hat)
#     if !(_lambda_0 >= sum((_x_hat[minmax(k,i)]+2*_x_hat[minmax(i,j)]+_x_hat[minmax(j,t)]-3)*_alpha[i,j,k,t] for (i,j,k,t) in find_index(_alpha))+ sum((_x_hat[minmax(i,j)]+_x_hat[minmax(j,k)]-1)*_beta[i,j,k] for (i,j,k) in find_index(_beta)))
#         cut = @constraint(master, lambda_0 >= sum((x[minmax(k,i)]+x[minmax(i,j)]+x[minmax(j,t)]-2)*_alpha[i,j,k,t] for (i,j,k,t) in find_index(_alpha))+ sum((x[minmax(i,j)]+x[minmax(j,k)]-1)*_beta[i,j,k] for (i,j,k) in find_index(_beta)))
#         @info "Adding the cut $(cut)"
#     end
# end

function _add_cut_SPi(master, _varphi, _gamma, indice, _lambda, _y_hat, result_dict)
    if !(_lambda[indice] >= 3(1-_y_hat[indice])*_varphi[indice] - sum(_y_hat[j]*_gamma[indice,j] for j in V if j!=indice))
        cut = @constraint(master, lambda[indice] >= 3(1-y[indice])*_varphi[indice] - sum(y[j]*_gamma[indice,j] for j in V if j!=indice))
        # @info "Adding the cut $(cut)"
        result_dict["num_cut_spi"] += 1
    end
end

function _add_cut_SP0_test(master, _alpha,_beta,_lambda_01, _lambda_02, _x_hat)
    if length(keys(_beta)) > 0 && !(_lambda_01 >= sum((_x_hat[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ _x_hat[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* _beta[(val[1],val[2],val[3])] for val in keys(_beta)))
        cut = @constraint(master, lambda_01 >= sum((x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ x[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* _beta[(val[1],val[2],val[3])] for val in keys(_beta)))
        # @info "Adding the cut lambda_01 $(cut)"
    end

    if length(keys(_alpha)) > 0 !(_lambda_02 >= sum((_x_hat[minmax(val[3],val[1])[1],minmax(val[3],val[1])[2]] + _x_hat[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+_x_hat[minmax(val[2],val[4])[1],minmax(val[2],val[4])[2]] - 2) * _alpha[(val[1],val[2],val[3],val[4])] for val in keys(_alpha)))
        cut = @constraint(master, lambda_02 >= sum((x[minmax(val[3],val[1])[1],minmax(val[3],val[1])[2]] + x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+x[minmax(val[2],val[4])[1],minmax(val[2],val[4])[2]] - 2) * _alpha[(val[1],val[2],val[3],val[4])] for val in keys(_alpha)))
        # @info "Adding the cut lambda_02 $(cut)"
    end
end



# function _add_knapsack_inequalities(master, _alpha, _beta, _gamma, _varphi, x_hat, y_hat, global_upper_bound)
#     con = @constraint(master, sum(ring_cost[i,j]*x[i,j] for (i,j) in E)+ sum(opening_cost[i]*y[i] for i in V) + sum(lambda[i] for i in V) + sum((x[minmax(k,i)]+x[minmax(i,j)]+x[minmax(j,t)] - 2)*_alpha[i,j,k,t]
#                         for (i,j,k,t) in find_index(_alpha))+ sum((x[minmax(i,j)]+x[minmax(j,k)] - 1)*_beta[i,j,k] for (i,j,k) in find_index(_beta)) + sum(3(1-y[indice])*_varphi[indice] - sum(y[j]*_gamma[indice,j] for j in V if j!=indice) for indice in V[2:end] if _varphi[indice]!= 0)
#                         <= global_upper_bound)
#     @info "Use global_upper_bound = $(global_upper_bound) \n Constraint: $(con)"
# end

function cal_obj_sp0(alpha, beta, x_hat)
    
    K_tilt_1 = keys(beta)
    J_tilt_1 = keys(alpha)

    if length(K_tilt_1) == 0
        return 0
    elseif length(J_tilt_1) == 0
        return  sum((x_hat[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ x_hat[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* beta[(val[1],val[2],val[3])] for val in K_tilt_1)
    else
        return sum((x_hat[minmax(val[3],val[1])[1],minmax(val[3],val[1])[2]] + x_hat[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+x_hat[minmax(val[2],val[4])[1],minmax(val[2],val[4])[2]] - 2) * alpha[(val[1],val[2],val[3],val[4])] for val in J_tilt_1) + sum((x_hat[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ x_hat[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* beta[(val[1],val[2],val[3])] for val in K_tilt_1)
    end
end
# val[1],val[2],val[3],val[4] = i,j,k,t

# function cal_obj_sp0(alpha, beta, x_hat)
#     K_tilt_1 = find_index(beta)
#     J_tilt_1 = find_index(alpha)
#     if length(V_tilt) == 1
#         return sum((x_hat[i,j]+ x_hat[j,k] - 1)* beta[i,j,k] for (i,j,k) in K_tilt_1)
#     elseif length(V_tilt) == 0
#         return 0
#     else
#         return sum((x_hat[k,i] + x_hat[i,j]+x_hat[j,t] - 2) * alpha[i,j,k,t] for (i,j,k,t) in J_tilt_1) + sum((x_hat[i,j]+ x_hat[j,k] - 1)* beta[i,j,k] for (i,j,k) in K_tilt_1)
#     end
# end


# Add cut SP0 (no splitting)
function _add_cut_SP0(master, _alpha,_beta,_lambda, sp0, result_dict)
    if !(_lambda >= sp0)
        cut = @constraint(master, lambda_0 >= sum((x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ x[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* _beta[(val[1],val[2],val[3])] for val in keys(_beta) if length(keys(_alpha))>0) + sum((x[minmax(val[3],val[1])[1],minmax(val[3],val[1])[2]] + x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+x[minmax(val[2],val[4])[1],minmax(val[2],val[4])[2]] - 2) * _alpha[(val[1],val[2],val[3],val[4])] for val in keys(_alpha) if length(keys(_alpha))>0))
        # @info "Adding the cut lambda_0 $(cut)"
        result_dict["num_cut_sp0"] += 1
    end
end

function _add_one_cut(master, sp0, sp1, _alpha, _beta, _lambda_one_hat, n, _varphi, _gamma, result_dict)
    if !(_lambda_one_hat >= sp0 + sp1)
        con = @constraint(master, lambda_one >= sum((x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+ x[minmax(val[2],val[3])[1],minmax(val[2],val[3])[2]] - 1)* _beta[(val[1],val[2],val[3])] for val in keys(_beta) if length(keys(_beta))>0) + sum((x[minmax(val[3],val[1])[1],minmax(val[3],val[1])[2]] + x[minmax(val[1],val[2])[1],minmax(val[1],val[2])[2]]+x[minmax(val[2],val[4])[1],minmax(val[2],val[4])[2]] - 2) * _alpha[(val[1],val[2],val[3],val[4])] for val in keys(_alpha) if length(keys(_alpha))>0)
                                                + sum(3(1-y[indice])*_varphi[indice] for indice in 1:n) - sum(y[j]*_gamma[indice,j] for indice in 1:n for j in 1:n if j!=indice))
        # @info "Adding only one cut $(con)"
        result_dict["num_cut_sp0"] += 1
    end
end

function initializeResult_dict(name)
    result_dict = Dict()
    result_dict["algorithm"] = name
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
                        "timestamp",
                        "gap",
                        "node_count"] # get the number of node 
        result_dict[param_name] = 0
    end
    return result_dict
end