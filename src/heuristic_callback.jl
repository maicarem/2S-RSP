"""
Input: structure of the ring 
Output: submit heuristic callback 
"""

function minmax(i,j)
    return min(i,j), max(i,j)
end

function length_delta(route, v1, v2, ring_cost)
    return -ring_cost[route[v1], route[v1+1]] - ring_cost[route[v2], route[v2+1]] + ring_cost[route[v1+1], route[v2+1]] + ring_cost[route[v1], route[v2]]
end

function two_opt_procedure(route, offset, ring_cost, )
    best_obj = cal_route(route, offset, ring_cost, backup_cost, V_tilt, transformation::Bool, cal_init::Bool)
    for i in 1:length(route)-1
        for j in i+1:length(route)-1
            length_delta(route, i, j, ring_cost) <0 || continue
            best_obj += length_delta(route, i, j, ring_cost)
            route = two_optSwap(route,i,j)
        end
    end
    return route
end

# Input: route1 = [1,2,3,4,1], cal_init: calculated as first hand
function cal_route(route1, offset, ring_cost, backup_cost, V_tilt, transformation::Bool, cal_init::Bool)
    
    _ring = 0
    _lambda_0 = cal_lambda_0(route1, backup_cost, V_tilt)
    
    if cal_init
        _ring = sum([ring_cost[minmax(route1[i],route1[i+1])] for i in 1:length(route1)-1])
    end

    return ifelse(transformation, _ring + _lambda_0+ offset,  _ring + _lambda_0)
end

function cal_lambda_0(route1, backup_cost, V_tilt)
    lambda_0 = 0
    saved0 = Dict{Tuple{Int, Int}, Float64}()

    for idx in 2:length(route1)-1
        
        (route1[idx] in V_tilt) && (minmax(route1[idx-1], route1[idx+1]) ∉ keys(saved0))| continue
        lambda_0 += backup_cost[minmax(route1[idx-1], route1[idx+1])]
        saved0[minmax(route1[idx-1], route1[idx+1])] = 1
        
        (route1[idx+1] in V_tilt) && (minmax(route1[idx-1], route1[idx+2]) ∉ keys(saved0)) || continue
        lambda_0 += backup_cost[minmax(route1[idx-1], route1[idx+2])]
        saved0[minmax(route1[idx-1], route1[idx+1])] = 1
    
    end

    return lambda_0
end
