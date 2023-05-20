function three_hub(n, V_certain, rc, sc, oc)
    global_obj = 1e18
    global_sol = []
    if length(V_certain) < 3
        return global_obj, global_sol
    else
        obj = oc[1]
        for i1 in V_certain[V_certain .!= 1] 
            for i2 in V_certain
                i2!= 1 && i2 > i1 || continue
                obj += oc[i1] + oc[i2] + rc[1, i1]+ rc[1, i2] + rc[i1, i2]
                for v in 2:n
                    v ∉ [i1, i2] || continue
                    obj += min(sc[v,1], sc[v,i1], sc[v, i2])
                end
                if obj < global_obj
                    global_obj = obj
                    global_sol = [1, i1, i2]
                end
            end
        end 
    end
    return global_obj, global_sol
end

function four_hubs(n, V_tilt, V_certain, rc, sc, oc)
    global_obj = 1e18
    global_sol = []
    if length(V_certain) < 3 
        return global_obj, global_sol
    else
        for i1 in V_certain[V_certain .!= 1] 
            for i2 in V_certain
                i2 != 1 && i2!= i1 || continue
                for i3 in 2:n
                    i3 ∉ [i1, i2] || continue
                    obj = oc[i1]+ oc[i2] + oc[i3] + oc[1]
                    # STAR cost
                    for v in 2:n
                        v ∉ [i1, i2, i3] || continue
                        if i3 ∈ V_tilt
                            obj += min(sc[v,1], sc[v,i1], sc[v,i2])
                        else 
                            obj += min(sc[v,1], sc[v,i1], sc[v,i2], sc[v,i3])
                        end
                    end

                    # RING cost
                    function cal_cost(v1,v2,v3,v4, rc, V_tilt)
                        obj1 = rc[v1, v2] + rc[v2, v3] + rc[v3, v4] + rc[v4, v1]
                        if v1 ∈ V_tilt
                            obj1 += rc[v2, v4]
                        end
                        return obj1
                    end
                    m = 1e18
                    local_solution = []
                    for (i,j, k,t) in [(i3, 1, i2, i1), (i3, i2, i1, 1), (i3, i1, 1, i2)]
                        if cal_cost(i,j,k,t, rc, V_tilt) < m
                            m = cal_cost(i,j,k,t, rc, V_tilt)
                            local_solution = [i,j,k,t]
                        end
                    end
                    obj += m
                    if obj < global_obj
                        global_obj = obj
                        global_sol = local_solution
                    end
                end
            end
        end

    end
    return global_obj, global_sol
end

function cal_star_cost_3uncertains(sc, v, H_tilt)
    return sum(sort(sc[v, H_tilt])[1:3])
end

function cal_tour_uncertain(rc, tour)
    return sum([rc[tour[i],tour[j]] for i in 1:4 for j in i+1:5])
end

function cal_tour_certain(rc, tour)
    return sum([rc[tour[i],tour[i+1]] for i in 1:4])+ rc[tour[5], tour[1]]
end

function cal_tour_mixed(rc, tour, H_tilt)
    cost0 = cal_tour_certain(rc, tour)

    for k in 3:4
        if tour[k] ∈ H_tilt
            cost0 += rc[k-1, k+1]
        end
    end
    if (tour[2] ∈ H_tilt && tour[3] ∈ H_tilt) || tour[5] ∈ H_tilt
        cost0 += rc[tour[1], tour[4]]
    end
    if tour[3] ∈ H_tilt && tour[4] ∈ H_tilt
        cost0 += rc[tour[2], tour[5]]
    end
    if tour[4] ∈ H_tilt && tour[5] ∈ H_tilt || tour[2] ∈ H_tilt
        cost0 += rc[tour[1], tour[3]]
    end
    return cost0
end



function five_hubs(n, V_tilt, rc, oc, sc)
    
    global_cost = 1e18
    global_sol = []

    for i2 in 1:n-3
        for i3 in i2+1:n-2
            for i4 in i2+1:n-1
                for i5 in i3+1:n
                    obj = oc[1] + oc[i2] + oc[i3] + oc[i4] + oc[i5]
                    H_tilt = []
                    H_certain = [1]
                    for m in [i2, i3, i4, i5]
                        if m ∈ V_tilt
                            push!(H_tilt, m)
                        else
                            push!(H_certain, m)
                        end
                    end
                    if length(H_tilt) < 3
                        for k in 2:n
                            if k ∉ [i2, i3, i4, i5]
                                obj += minimum(sc[k, v] for v in H_certain)   
                            end
                        end
                    elseif length(H_certain) == 0
                        for k in 2:n
                            if k ∉ [i2, i3, i4, i5]
                                obj += cal_star_cost_3uncertains(sc, k, H_tilt)
                            end
                        end
                    else
                        for k in 2:n
                            if k ∉ [i2, i3, i4, i5]
                                obj += min(cal_star_cost_3uncertains(sc, k, H_tilt), minimum(sc[k, v] for v in H_certain))
                            end
                        end
                    end
                    
                    local_cost = 1e18
                    local_solution = []
                    
                    if length(H_tilt) == 4
                        for tour in [[1, i2, i3, i4, i5],
                                    [1, i2, i4, i3, i5],
                                    [1, i2, i5, i3, i4],
                                    [1, i3, i2, i4, i5],
                                    [1, i3, i4, i2, i5],
                                    [1, i3, i5, i2, i4],
                                    [1, i4, i2, i3, i5],
                                    [1, i4, i3, i2, i5],
                                    [1, i4, i5, i2, i3],
                                    [1, i5, i2, i3, i4],
                                    [1, i5, i3, i2, i4],
                                    [1, i5, i4, i2, i3]]
                            cost0 = cal_tour_uncertain(rc, tour)
                            if cost0 < local_cost
                                local_cost = cost0
                            end
                        end
                    elseif length(H_tilt) == 0
                        for tour in [[1, i2, i3, i4, i5],
                            [1, i2, i4, i3, i5],
                            [1, i2, i5, i3, i4],
                            [1, i3, i2, i4, i5],
                            [1, i3, i4, i2, i5],
                            [1, i3, i5, i2, i4],
                            [1, i4, i2, i3, i5],
                            [1, i4, i3, i2, i5],
                            [1, i4, i5, i2, i3],
                            [1, i5, i2, i3, i4],
                            [1, i5, i3, i2, i4],
                            [1, i5, i4, i2, i3]]
                            cost0 = cal_tour_certain(rc, tour)
                            if cost0 < local_cost
                                local_cost = cost0
                            end
                        end
                    else
                        for tour in [[1, i2, i3, i4, i5],
                            [1, i2, i4, i3, i5],
                            [1, i2, i5, i3, i4],
                            [1, i3, i2, i4, i5],
                            [1, i3, i4, i2, i5],
                            [1, i3, i5, i2, i4],
                            [1, i4, i2, i3, i5],
                            [1, i4, i3, i2, i5],
                            [1, i4, i5, i2, i3],
                            [1, i5, i2, i3, i4],
                            [1, i5, i3, i2, i4],
                            [1, i5, i4, i2, i3]]
                            cost0 = cal_tour_mixed(rc, tour, H_tilt)
                            if cost0 < local_cost
                                local_cost = cost0
                            end
                        end
                    end
                    local_cost += obj
                    if local_cost < global_cost
                        global_cost = local_cost
                        global_sol = [1, i2, i3, i4, i5]
                    end
                end
            end
        end
    end
    return global_cost, global_sol
end


using JuMP, Gurobi
using Combinatorics
using Graphs, GraphsFlows

include("dat.jl")
include("misc.jl")
include("write_output.jl")
include("mutable_structure.jl")
include("user_cut.jl")

# Initialize Sets
pars = MainPar(uc_strat = 3, transformation = false, alpha = 3)
name = "instances/small_instances/small_instance_10.dat"
n, oc, sc, rc = read_input_random(name, pars)
V, V_tilt, V_certain, A, A_prime, E, T_tilt, J_tilt, K_tilt = _declare_set(n, pars)
@show five_hubs(n, V_tilt, rc, oc, sc)
@show four_hubs(n, V_tilt, V_certain, rc, sc, oc)
@show three_hub(n, V_certain, rc, sc, oc)
