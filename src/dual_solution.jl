"""
Find dual solution
Input: y_hat, x_hat: structure of the ring
"""

function _dual_backup_edges(y_hat, x_hat, backup_edge_1, backup_edge_2)
    
    alpha = zeros(Float64, n,n,n,n)
    beta = zeros(Float64, n,n,n)
    
    beta_prime = Dict()
    saved0 = Dict{Tuple{Int, Int}, Float64}()
    
    for j in V_tilt
        beta_prime[j] = []
    end
    # Find adjacent list in uncertain hubs
    for i in V
        for j in V
            i!=j || continue
            if round(x_hat[i,j]) == 1 && i ∈ V_tilt && j ∉ beta_prime[i]
                append!(beta_prime[i], j)
            elseif round(x_hat[i,j]) == 1 && j ∈ V_tilt && i ∉ beta_prime[j]
                append!(beta_prime[j], i)
            end
        end
    end

    # Find dual beta[i,j,k], i<k
    for i in V_tilt
        y_hat[i] == 1 && (beta_prime[i][1], beta_prime[i][2]) ∉ keys(saved0) || continue
        beta[beta_prime[i][1],i, beta_prime[i][2]] = backup_edge_1[beta_prime[i][1], beta_prime[i][2]]
        saved0[(beta_prime[i][1], beta_prime[i][2])] = 1
    end
    
    # Find alpha[i,j,k,t], k<t
    for i in V_tilt
        y_hat[i] == 1 || continue
        for j in beta_prime[i]
            j ∈ V_tilt || continue
            _m0 = [k for k in beta_prime[i] if k!=j][1]
            _n0 = [k for k in beta_prime[j] if k!=i][1]
            m0 = min(_m0, _n0)
            n0 = max(_m0, _n0)
            (m0,n0) ∉ keys(saved0) || continue
            if x_hat[m0,i] == 1 && x_hat[j,n0] == 1
                alpha[i,j,m0,n0] = backup_edge_2[m0,n0]
            elseif x_hat[m0,j] == 1 && x_hat[i,n0] == 1
                alpha[j,i,m0,n0] = backup_edge_2[m0,n0]
            end
            saved0[(m0,n0)] = 1
        end
    end
    return beta, alpha
end
# Find φ_i in SP_i
function _cal_varphi(i, star_cost0)
    if length(V_tilt) == n
        return minimum(star_cost0[i,V_tilt])
    elseif length(V_certain) == n
        return 1/3*minimum(star_cost0[i, V_certain])
    else
        return min(1/3*minimum(star_cost0[i, V_certain]), minimum(star_cost0[i,V_tilt])) 
    end
end

function _dual_star_structure(_y_hat, star_cost0)
    H, H_tilt = _list_hub(_y_hat)
    uncertain_dist, _ = _uncertain_dict(H, star_cost0)
    # @show uncertain_dist
    φ, γ = zeros(Float64, n), zeros(Float64, n,n)

    for i in 1:n
        if i in H
            φ[i] = _cal_varphi(i, star_cost0) # generate cut if i is a hub
        else
            x = min(certain_smallest(H, i, star_cost0)[2], uncertain_dist[i][1][2]+uncertain_dist[i][2][2]+uncertain_dist[i][3][2])
            φ[i] = x/3
            for j in V_certain # find γ and φ if i is a terminal
                j!= i || continue
                γ[i,j] = max(0, x-star_cost0[i,j])
            end
            for j in V_tilt
                j ∉ H_tilt && j!= i || continue
                γ[i,j] = max(0, 1/3 * x - star_cost0[i,j])
            end
        end
    end
    return φ, γ
end

function dual_solution(y_hat, x_hat, backup_edge_1, backup_edge_2, star_cost0) 
    # Output: beta, alpha, φ, γ
    return _dual_backup_edges(y_hat, x_hat, backup_edge_1, backup_edge_2), _dual_star_structure(y_hat, star_cost0)
end

# y_hat = [1.0  0.0  0.0  0.0  0.0  0.0;
# 0.0  0.0  0.0  0.0  0.0  0.0;
# 0.0  0.0  1.0  0.0  0.0  0.0;
# 0.0  0.0  0.0  1.0  0.0  0.0;
# 0.0  0.0  0.0  0.0  1.0  0.0;
# 0.0  0.0  0.0  0.0  0.0  1.0]

# x_hat = [0.0  0.0  1.0  1.0  0.0  0.0;
# 0.0  0.0  0.0  0.0  0.0  0.0;
# 1.0  0.0  0.0  0.0  0.0  1.0;
# 1.0  0.0  0.0  0.0  1.0  0.0;
# 0.0  0.0  0.0  1.0  0.0  1.0;
# 0.0  0.0  1.0  0.0  1.0  0.0]

# Test 2


y_hat = [
    1.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0;
    0.0   1.0  -0.0  -0.0  -0.0  -0.0   0.0  -0.0   0.0  -0.0;
   -0.0  -0.0   1.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0;
    1.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0;
    1.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0;
    1.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0;
   -0.0  -0.0  -0.0  -0.0  -0.0  -0.0   1.0  -0.0  -0.0  -0.0;
    1.0  -0.0   0.0  -0.0  -0.0  -0.0  -0.0   0.0  -0.0  -0.0;
   -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0   1.0  -0.0;
    1.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0  -0.0]

# @show _list_hub(y_hat)
# @show _dual_star_structure(y_hat, star_cost)
# @show minimum(sc[2, V_certain])
# @show dual_solution(y_hat, x_hat)
# @show _dual_backup_edges(y_hat, x_hat, backup)