# Data Raw
n = 6
sc = ones(n,n)
for i = 1:n
    sc[i,4] = 0.2
    sc[4,i] = 0.2
    sc[i,i] = 1e18
end
rc = sc*2
for i = 1:n
    rc[i,4] = 0.1
    rc[4,i] =0.1
    rc[i,i] = 1e18
end
oc = zeros(n)
ABSOLUTE_OPTIMALITY_GAP = 1e-6

# Declare sets
V = 1:n
V_tilt = [1,2,3,4]
V_certain = [i for i in V if i ∉ V_tilt]
A = [(i,j) for i in V for j in V if i !=j]
A_prime = [(i,j) for i in V for j in V]
E = [(i,j) for i in V for j in V if i<j]
J_tilt = [(i,j,k,t) for k in V for t in V for i in V_tilt for j in V_tilt if k<t && i!=j && i!=k && i!=t && j!=k && j!=t]
K_tilt = [(i,j,k) for k in V for i in V for j in V_tilt if i<k && i!=j && j!=k]
T_tilt = [(i,j) for i in V_tilt for j in V_tilt if i<j]

# Instance transformation
function _find_offset_star()
    offset = Dict{Int, Float64}()
    for i in 1:n
        if length(V_tilt) <= 1
            offset[i] = minimum(3 * sc[i,[j for j in 1:n if i!=j && j in V_certain]])
        elseif length(V_certain) <= 1
            offset[i] = minimum(sc[i,[j for j in 1:n if i!=j && j in V_tilt]])
        else
            offset[i] = min(minimum(3 * sc[i,[j for j in 1:n if j!=i && j in V_certain]]), minimum(sc[i,[j for j in 1:n if j!=i && j in V_tilt]]))
        end
    end
    return offset
end

function _find_offset_backup()
    offset_backup = Dict{Int, Float64}()
    for i in 1:n
        offset_backup[i] = minimum(rc[1:n .!= i, i])
    end
    return offset_backup
end

function _transformation_cost()
    
    star_offset_dict = _find_offset_star()
    backup_offset_dict = _find_offset_backup()
    offset = sum([star_offset_dict[i] for i in 1:n])
    
    opening_cost = Dict{Int, Float64}()
    ring_cost = ones(Float64, n,n) * 1e18
    star_cost = ones(Float64, n,n) * 1e18
    backup_cost = ones(Float64, n,n) * 1e18

    for i in 1:n
        opening_cost[i] = oc[i] - star_offset_dict[i]
        for j in 1:n
            i<j || continue
            if i in V_tilt && j in V_certain
                ring_cost[i,j] = rc[i,j] + 1/2 * backup_offset_dict[j]
            elseif i in V_certain && j in V_tilt
                ring_cost[i,j] = rc[i,j] + 1/2 * backup_offset_dict[i]
            elseif i in V_tilt && j in V_tilt
                ring_cost[i,j] = rc[i,j] + 1/2 * backup_offset_dict[i] + 1/2 * backup_offset_dict[j]
            else
                ring_cost[i,j] = rc[i,j]
            end
            
            backup_cost[i,j] = rc[i,j] - 1/2 * backup_offset_dict[i] - 1/2 * backup_offset_dict[j]
            backup_cost[j,i] = backup_cost[i,j]
            if j in V_certain
                j!=i || continue
                star_cost[i,j] = sc[i,j] - star_offset_dict[i]
            else
                j!=i || continue
                star_cost[i,j] = sc[i,j] - 1/2 * star_offset_dict[i]
            end
            ring_cost[j,i] = ring_cost[i,j]
            star_cost[j,i] = star_cost[i,j]
        end
    end
    return offset, opening_cost, ring_cost, star_cost, backup_cost
end

offset, oc, rc, sc, backup = _transformation_cost()
@show backup