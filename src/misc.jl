
# Find list of indices j having s_ij<= s_imi2
# unhub0: [(ind[i], cost[i]) for i in 1:3]
function find_omega_i(unhub0, hub_list, indice)
    _cost = unhub0[2][2]
    if _cost == 1e18
        _cost = unhub0[1][2]
        if _cost == 1e18
            _cost = 0
        end
    end
    return [i for i in 1:length(sc[indice,:]) if sc[indice,i] <= _cost && i ∈ V_tilt && i ∈ hub_list]
end

# Find first, second and third smallest distance of s_ij, j is uncertain and a hub:: return [(hub1, distance1),...]
# Input: hub_list, V_tilt
function uncertain_smallest(hub_list, indice)
    n = length(sc[indice,:])
    sorted0 = sort([sc[indice, i] for i in 1:n if i in V_tilt && i in hub_list && i != indice])
    ind, cost = zeros(Int, 3), ones(3)*1e18
    if length(sorted0) != 0
        for i in 1: min(3, length(sorted0))
            ind[i] = [j for j in 1:n if sc[indice, j] == sorted0[i] && j in hub_list && j in V_tilt && j ∉ ind && j != indice][1]
            cost[i] = sc[indice, ind[i]]
        end
    end
    return [(ind[i], cost[i]) for i in 1:3]
end

# Find indice mi^star and smallest cost 
function certain_smallest(hub_list, indice)
    n = length(sc[indice,:])
    sorted0 = sort([sc[indice, i] for i in 1:n if i ∉ V_tilt && i in hub_list && i != indice])
    if length(sorted0)!=0
        ind = Int([j for j in 1:n if sc[indice, j] == sorted0[1] && j in hub_list && j ∉ V_tilt && j != indice][1])
        cost = sc[indice, ind]
    else
        ind, cost = 0, 1e18
    end
    return (ind, cost)
end

function _uncertain_dict(hub_list)
    _smallest_uncertain_dict = Dict()
    _Omega_dict = Dict()
    for i in 1:n
        if i ∉ hub_list
            _smallest_uncertain_dict[i] = uncertain_smallest(hub_list, i)
            _Omega_dict[i] = find_omega_i(_smallest_uncertain_dict[i], hub_list, i)
        end
    end
    return _smallest_uncertain_dict, _Omega_dict
end

function _list_hub(_y_hat)
    hub = []
    hub_tilt = []
    for i in 1:n
        _y_hat[i,i] == 1 || continue
        push!(hub, i)
        i ∈ V_tilt || continue
        push!(hub_tilt, i)
    end   
   return hub, hub_tilt 
end

# include("dat.jl")
# hub = [1,2,6,3,4,5]
# @show _uncertain_dict(hub)