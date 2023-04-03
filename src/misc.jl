using LightGraphs
include("dat.jl")

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

function contain_subtour(x_hat, y_hat)

    visited = zeros(Int, n)
    parent = zeros(Int, n)
    vistime = zeros(Int, n)
    adjacent_list = Dict()
    
    timer = 0
    cycle_count = 0

    for i in 1:n
        if y_hat[i] == 1
            adjacent_list[i] = []
            for j in 1:n
                if x_hat[i,j] == 1
                    append!(adjacent_list[i], j)
                end
            end
        else
            visited[i] = 1
        end
    end
    
    visited[1] = 1

    function dfs(node)
        
        for next in adjacent_list[node]
            if visited[next] == 0
                visited[next] = 1
                timer +=1
                vistime[next] = timer
                parent[next] = node
                dfs(next)
            elseif visited[next] == 1 && parent[node] == next
                continue
            elseif visited[next] == 1 && parent[node] != next && vistime[next]< timer
                cycle_count +=1
                continue
            end
        end
    end
    dfs(1)

    if sum(visited) != length(visited)
        return false
    else
        return cycle_count == 1
    end
    
end

function cal_obj_sp0(alpha, beta, x_hat)
    if length(V_tilt) == 1
        return sum((x_hat[i,j]+ x_hat[j,k] - 1)* beta[i,j,k] for (i,j,k) in K_tilt)
    elseif length(V_tilt) == 0
        return 0
    else
        return sum((x_hat[k,i] + x_hat[i,j]+x_hat[j,t] - 2) * alpha[i,j,k,t] for (i,j,k,t) in J_tilt) + sum((x_hat[i,j]+ x_hat[j,k] - 1)* beta[i,j,k] for (i,j,k) in K_tilt)
    end
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

function minmax(i,j)
    return min(i,j), max(i,j)
end

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
