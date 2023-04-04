include("read_input.jl")
n, oc, sc, rc = read_input("instances/small_instances/small_instance_10.dat")


ABSOLUTE_OPTIMALITY_GAP = 1e-6
bench = 0.5

# Declare sets
V = 1:n
V_tilt = 2:n
V_certain = [i for i in V if i ∉ V_tilt]
A = [(i,j) for i in V for j in V if i !=j]
A_prime = [(i,j) for i in V for j in V]
E = [(i,j) for i in V for j in V if i<j]
J_tilt = [(i,j,k,t) for k in V for t in V for i in V_tilt for j in V_tilt if k<t && i!=j && i!=k && i!=t && j!=k && j!=t]
K_tilt = [(i,j,k) for k in V for i in V for j in V_tilt if i<k && i!=j && j!=k]
T_tilt = [(i,j) for i in V_tilt for j in V_tilt if i<j]

# Instance transformation
function _find_offset_star()
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

function _find_offset_backup()
    offset_backup = zeros(Float64, n)
    for i in 1:n
        offset_backup[i] = minimum(rc[1:n .!= i, i])
    end
    return offset_backup
end

function _transformation_cost()
    
    star_offset_list = _find_offset_star()
    backup_offset_list = _find_offset_backup()
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

opening_cost, ring_cost, star_cost = oc, rc, sc
offset, oc, rc, sc, backup = _transformation_cost()