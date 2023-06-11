using LinearAlgebra
using DelimitedFiles

function read_input_random(filepath, pars0)
    
    f = open(filepath, "r")
    _n,_,_ = [parse(Float64, i) for i in split(readline(f)," ")]
    _n = Int(_n)

    coor = zeros(_n,2)
    dist = ones(_n,_n) * 1e18
    for i in 1:_n
        _, coor[i,1], coor[i,2] = [parse(Float64, i) for i in split(readline(f)," ")]
    end

    # calculate star cost
    for i in 1:_n
        for j in 1:_n
            i<j || continue
            dist[i,j] = round(norm([coor[i,1] - coor[j,1], coor[i,2] - coor[j,2]], 2))
            dist[j,i] = dist[i,j]
        end
    end
    star_cost = ceil.(Int, dist * pars.alpha)
    ring_cost = ceil.(Int, (10- pars.alpha) * dist)
    # star_cost = dist
    # ring_cost = dist * 2
    readline(f)
    opening_cost = [parse(Float64, i) for i in split(readline(f)," ")]
    close(f)
    return _n, opening_cost, star_cost, ring_cost # number of nodes, opening cost, star cost, ring cost
end

# Declare sets
function _declare_set(n, pars)
    V = 1:n
    V_tilt = 2:n
    V_certain = [i for i in V if i ∉ V_tilt]
    A = [(i,j) for i in V for j in V if i !=j]
    A_prime = [(i,j) for i in V for j in V]
    E = [(i,j) for i in V for j in V if i<j]
    T_tilt = [(i,j) for i in V_tilt for j in V_tilt if i<j]
    if pars.benders
        J_tilt = []
        K_tilt = []
    else
        J_tilt = [(i,j,k,t) for k in V for t in V for i in V_tilt for j in V_tilt if k<t && i!=j && i!=k && i!=t && j!=k && j!=t]
        K_tilt = [(i,j,k) for k in V for i in V for j in V_tilt if i<k && i!=j && j!=k]
    end
    return V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt
end

# Instance transformation
function _find_offset_star(n, V_tilt, V_certain, sc)
    offset = zeros(Float64, n)
    for i in 1:n
        _V_tilt = [j for j in V_tilt if j!=i]
        _V_certain = [j for j in V_certain if j!=i]
        lb1, lb2 = 1e18, 1e18
        if length(_V_certain) == 0
            lb2 = 3 * minimum([sc[i,k] for k in _V_tilt])
        elseif length(_V_tilt) == 0
            lb1 = minimum([sc[i,k] for k in _V_certain])
        else
            lb1 = minimum([sc[i,k] for k in _V_certain])
            lb2 = 3 * minimum([sc[i,k] for k in _V_tilt])
        end
        offset[i] = min(lb1, lb2)
    end
    return offset
end

function _find_offset_backup(n, rc)
    offset_backup = zeros(Float64, n)
    for i in 1:n
        offset_backup[i] = minimum(rc[1:n .!= i, i])
    end
    return offset_backup
end

function _transformation_cost(rc, sc, oc, n, V_tilt, V_certain)
    
    star_offset_list = _find_offset_star(n, V_tilt, V_certain, sc)
    backup_offset_list = _find_offset_backup(n, rc)
    offset = sum([star_offset_list[i] for i in 1:n])
    
    opening_cost = zeros(Float64, n)
    ring_cost = ones(Float64, n,n) * 1e18
    star_cost = ones(Float64, n,n) * 1e18
    backup_cost = ones(Float64, n,n) * 1e18

    for i in 1:n
        
        opening_cost[i] = oc[i] - star_offset_list[i]
        
        for j in 1:n
            
            i<j || continue
            if i in V_tilt && j in V_certain
                ring_cost[i,j] = rc[i,j] + 1/2 * backup_offset_list[j]
            elseif i in V_certain && j in V_tilt
                ring_cost[i,j] = rc[i,j] + 1/2 * backup_offset_list[i]
            elseif i in V_tilt && j in V_tilt
                ring_cost[i,j] = rc[i,j] + 1/2 * (backup_offset_list[i] + backup_offset_list[j])
            else
                ring_cost[i,j] = rc[i,j]
            end
            
            backup_cost[i,j] = rc[i,j] - 1/2 * backup_offset_list[i] - 1/2 * backup_offset_list[j]
            ring_cost[j,i] = ring_cost[i,j]
            backup_cost[j,i] = backup_cost[i,j]
        end
        
        for j in V_certain
            j!=i || continue
            star_cost[i,j] = sc[i,j] - star_offset_list[i]
        end
        
        for j in V_tilt
            j!=i || continue
            star_cost[i,j] = sc[i,j] - 1/3 * star_offset_list[i]
        end
    end
    # Check again the instance transformation for the backup cost, now revert to unchanged
    return offset, opening_cost, ring_cost, star_cost, backup_cost
    # return offset, opening_cost, rc, star_cost, rc # <- this one is unchanged
end
# V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt = _declare_set(n, 1)
# opening_cost, ring_cost, star_cost = oc, rc, sc
# offset, oc, rc, sc, backup = _transformation_cost(rc,sc, oc, n, V_tilt, V_certain)

function read_input_journal(filepath, pars)
    
    if readdlm(filepath)[7] == "EDGE_WEIGHT_SECTION"
        n = length(readdlm(filepath)[8, 1:end])
        dist = readdlm(filepath)[8:7+n, 1:end]
        
        for i in 1:n
            for j in i+1:n
                dist[i,j] = dist[j,i]
            end
            dist[i,i] = 1e18
        end

        star_cost = ceil.(Int, dist * pars.alpha)
        ring_cost = ceil.(Int, (10- pars.alpha)* dist)
        opening_cost = zeros(Int, n)
    else
        coor = readdlm(filepath)[7:end-1, 2:3]
        n = size(coor)[1]
        dist = ones(n,n) * 1e18
        
        for i in 1:n
            for j in 1:n
                i<j || continue
                dist[i,j] = round(norm([coor[i,1] - coor[j,1], coor[i,2] - coor[j,2]], 2))
                dist[j,i] = dist[i,j]
            end
        end

        star_cost = ceil.(Int, dist * pars.alpha)
        ring_cost = ceil.(Int, (10- pars.alpha)* dist)
        opening_cost = zeros(Int, n)
    end
    return n, opening_cost, star_cost, ring_cost # number of nodes, opening cost, star cost, ring cost
end

# filepath = "instances/small_instances/Instances_journal_article/EUC_2D/berlin52.tsp"
# filepath = "instances/small_instances/Instances_journal_article/EUC_2D/bier127.tsp"

# @show read_input_journal(filepath, 7)
# filepath = "instances/small_instances/Instances_journal_article/non-EUC_2D/ulysses16.tsp"
# include("mutable_structure.jl")
# pars = MainPar(transformation = false, one_cut = true)
# read_input_journal(filepath, pars)
