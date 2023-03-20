# Data Raw
n = 10
sc = [0.286709  0.621484   0.634555  0.549575   0.702369   0.272015   0.4792      0.616595    0.114736  0.817952;
0.913911  0.0825983  0.263735  0.562878   0.423548   0.80518    0.792842    0.0409161   0.436656  0.158712;
0.701553  0.579373   0.882798  0.165368   0.255035   0.941467   0.232989    0.0288358   0.337281  0.724993;
0.297348  0.187299   0.111605  0.0534248  0.774593   0.245285   0.131914    0.571874    0.784561  0.269658;
0.425547  0.799424   0.432703  0.743771   0.591218   0.30639    0.0192837   0.568855    0.248702  0.371714;
0.107852  0.649824   0.151641  0.701445   0.0981913  0.365841   0.834072    0.16675     0.27794   0.831993;
0.242974  0.197402   0.395388  0.842056   0.615189   0.778668   0.00807703  0.465702    0.406278  0.477218;
0.89664   0.406671   0.869833  0.812098   0.295635   0.634037   0.355542    0.93294     0.869041  0.0110682;
0.736589  0.91335    0.422185  0.379827   0.893815   0.781795   0.0732774   0.0864879   0.512766  0.988308;
0.453407  0.370305   0.953327  0.0295641  0.493021   0.0193349  0.0137156   0.00421785  0.995167  0.0307335]
sc = sc*10
# ones(n,n)
for i = 1:n
    sc[i,4] = 0.2
    sc[4,i] = 0.2
    sc[i,i] = 1e18
end
rc = sc*2
for i = 1:n
    rc[i,4] = 0.1
    rc[4,i] = 0.1
    rc[i,i] = 1e18
end
oc = zeros(n)
ABSOLUTE_OPTIMALITY_GAP = 1e-6

# Declare sets
V = 1:n
V_tilt = [1,2,3,4,5,6,7]
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
    return offset, opening_cost, ring_cost, star_cost, backup_cost
end

opening_cost, ring_cost, star_cost = oc, rc, sc
offset, oc, rc, sc, backup = _transformation_cost()
lll=1