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

function dual_solution(y_hat, x_hat, V_tilt, n, backup_edge_1, backup_edge_2, star_cost0) 
    # Output: beta, alpha, φ, γ
    # return _dual_backup_edges(y_hat, x_hat, backup_edge_1, backup_edge_2), _dual_star_structure(y_hat, star_cost0)
    return dual_solution_sp0(x_hat, n, V_tilt, backup_edge_1, backup_edge_2), _dual_star_structure(y_hat, star_cost0)
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



# input: x_hat
# output: alpha, beta
function dual_solution_sp0(x_hat, n, V_tilt, backup1, backup2)
    # Initialize
    tour = transform_route(x_hat)
    saved0 = zeros(n,n)
    alpha = Dict()
    beta = Dict()

    # Calculate beta 
    for idx_node in range(2,length(tour)-1)
        if tour[idx_node] in V_tilt
            node1, node2 = minmax(tour[idx_node-1],tour[idx_node+1])
            if saved0[node1, node2] == 0
                beta[(node1, tour[idx_node], node2)] = backup1[node1, node2]
                saved0[node1, node2] = 1
            end
        end
    end

    for idx_node in range(2,length(tour)-2)
        if tour[idx_node] in V_tilt && tour[idx_node+1] in V_tilt
            node1, node2 = minmax(tour[idx_node - 1], tour[idx_node + 2])
            if saved0[node1, node2] == 0
                alpha[(tour[idx_node], tour[idx_node + 1], tour[idx_node - 1], tour[idx_node + 2])] = backup2[node1, node2]
                saved0[node1, node2] = 1
            end
        end
    end

    # for idx_node in range(2,length(tour)-1)
    #     tour[idx_node] in V_tilt || continue
    #     for node_prime in range(1, n)
    #         if node_prime ∉ [tour[idx_node+1], tour[idx_node-1], tour[idx_node]]
    #             node1_prime, node2_prime = minmax(tour[idx_node - 1], node_prime)
    #             if saved0[node1_prime, node2_prime] == 0
    #                 beta[(node1_prime, tour[idx_node], node2_prime)] = backup1[node1_prime, node2_prime]
    #                 saved0[node1_prime, node2_prime] = 1
    #             end
    #         end
    #     end
    # end

    # for idx_node in range(2,length(tour)-2)
    #     if tour[idx_node] in V_tilt && tour[idx_node+1] in V_tilt
    #         for node_prime in range(1,n)
    #             if node_prime ∉ [tour[idx_node], tour[idx_node-1], tour[idx_node+1], tour[idx_node+2]]
    #                 node1, node2 = minmax(node_prime, tour[idx_node-1])
    #                 if saved0[node1, node2] == 0
    #                     alpha[(tour[idx_node], tour[idx_node+1], tour[idx_node-1], node_prime)] = backup[node1, node2]
    #                     saved0[node1, node2] = 1
    #                 end
    #             end
    #         end
    #     end
    # end
    return beta, alpha
end

# x_hat = zeros(Int, 8,8)
# for (i,j) in [(1,3), (3,4), (4,6), (6,5), (5,2), (2,1)]
#     x_hat[i,j] = 1
#     x_hat[j,i] = 1
# end
# backup = ones(Int, 8, 8)
# beta, alpha = dual_solution_sp0(x_hat, 8, [2,3,4,5,6,7,8], backup, ones(Int,8,8))
# # test = sum((x_hat[minmax(i,j)[1], minmax(i,j)[2]]+ x_hat[minmax(j,k)[1], minmax(j,k)[2]]-1)*beta[(i,j,k)] for (i,j,k) in keys(beta))
# test = sum((x_hat[minmax(i,k)[1], minmax(i,k)[2]] + x_hat[minmax(k,t)[1], minmax(k,t)[2]] + x_hat[minmax(t,j)[1], minmax(t,j)[2]] -2) * alpha[(k,t,i,j)] for (k,t,i,j) in keys(alpha))