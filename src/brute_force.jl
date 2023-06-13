function three_hubs(n, V_certain, ring_cost, star_cost, opening_cost)
    t = time()
    bf_global_obj = 1e18
    global_sol = []
    if length(V_certain) < 3
        return bf_global_obj, global_sol, time() - t
    else
        obj = opening_cost[1]
        for i1 in V_certain[V_certain .!= 1] 
            for i2 in V_certain
                i2!= 1 && i2 > i1 || continue
                obj += opening_cost[i1] + opening_cost[i2] + ring_cost[1, i1]+ ring_cost[1, i2] + ring_cost[i1, i2]
                for v in 2:n
                    v ∉ [i1, i2] || continue
                    obj += min(star_cost[v,1], star_cost[v,i1], star_cost[v, i2])
                end
                if obj < bf_global_obj
                    bf_global_obj = obj
                    global_sol = [1, i1, i2]
                end
            end
        end 
    end
    return bf_global_obj, global_sol, time() - t
end

function four_hubs(n, V_tilt, V_certain, ring_cost, star_cost, opening_cost)
    t = time()
    bf_global_obj = 1e18
    global_sol = []
    if length(V_certain) < 3 
        return bf_global_obj, global_sol, time() - t
    else
        for i1 in V_certain[V_certain .!= 1] 
            for i2 in V_certain
                i2 != 1 && i2!= i1 || continue
                for i3 in 2:n
                    i3 ∉ [i1, i2] || continue
                    obj = opening_cost[i1]+ opening_cost[i2] + opening_cost[i3] + opening_cost[1]
                    # STAR cost
                    for v in 2:n
                        v ∉ [i1, i2, i3] || continue
                        if i3 ∈ V_tilt
                            obj += min(star_cost[v,1], star_cost[v,i1], star_cost[v,i2])
                        else 
                            obj += min(star_cost[v,1], star_cost[v,i1], star_cost[v,i2], star_cost[v,i3])
                        end
                    end

                    # RING cost
                    function cal_cost(v1,v2,v3,v4, ring_cost, V_tilt)
                        obj1 = ring_cost[v1, v2] + ring_cost[v2, v3] + ring_cost[v3, v4] + ring_cost[v4, v1]
                        if v1 ∈ V_tilt
                            obj1 += ring_cost[v2, v4]
                        end
                        return obj1
                    end
                    m = 1e18
                    lopening_costal_solution = []
                    for (i,j, k,t) in [(i3, 1, i2, i1), (i3, i2, i1, 1), (i3, i1, 1, i2)]
                        if cal_cost(i,j,k,t, ring_cost, V_tilt) < m
                            m = cal_cost(i,j,k,t, ring_cost, V_tilt)
                            lopening_costal_solution = [i,j,k,t]
                        end
                    end
                    obj += m
                    if obj < bf_global_obj
                        bf_global_obj = obj
                        global_sol = lopening_costal_solution
                    end
                end
            end
        end

    end
    return bf_global_obj, global_sol, time() - t
end

function cal_star_cost_3uncertains(star_cost, v, H_tilt)
    if length(H_tilt)<3
        return 1e18
    else
        return sum(sort(star_cost[v, H_tilt])[1:3])
    end
end

function cal_tour_uncertain(ring_cost, tour)
    return sum([ring_cost[tour[i],tour[j]] for i in 1:4 for j in i+1:5])
end

function cal_tour_certain(ring_cost, tour)
    return sum([ring_cost[tour[i],tour[i+1]] for i in 1:4])+ ring_cost[tour[5], tour[1]]
end

function cal_tour_mixed(ring_cost, tour, H_tilt)
    cost0 = cal_tour_certain(ring_cost, tour)

    for k in 3:4
        if tour[k] ∈ H_tilt
            cost0 += ring_cost[k-1, k+1]
        end
    end
    if (tour[2] ∈ H_tilt && tour[3] ∈ H_tilt) || tour[5] ∈ H_tilt
        cost0 += ring_cost[tour[1], tour[4]]
    end
    if tour[3] ∈ H_tilt && tour[4] ∈ H_tilt
        cost0 += ring_cost[tour[2], tour[5]]
    end
    if tour[4] ∈ H_tilt && tour[5] ∈ H_tilt || tour[2] ∈ H_tilt
        cost0 += ring_cost[tour[1], tour[3]]
    end
    return cost0
end



function five_hubs(n, V_tilt, ring_cost, opening_cost, star_cost)
    t = time()
    global_cost = 1e18
    global_sol = []

    for i2 in 1:n-3
        for i3 in i2+1:n-2
            for i4 in i3+1:n-1
                for i5 in i4+1:n
                    obj = opening_cost[1] + opening_cost[i2] + opening_cost[i3] + opening_cost[i4] + opening_cost[i5]
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
                                obj += minimum(star_cost[k, v] for v in H_certain)   
                            end
                        end
                    elseif length(H_certain) == 0
                        for k in 2:n
                            if k ∉ [i2, i3, i4, i5]
                                obj += cal_star_cost_3uncertains(star_cost, k, H_tilt)
                            end
                        end
                    else
                        for k in 2:n
                            if k ∉ [i2, i3, i4, i5]
                                obj += min(cal_star_cost_3uncertains(star_cost, k, H_tilt), minimum(star_cost[k, v] for v in H_certain))
                            end
                        end
                    end
                    
                    lopening_costal_cost = 1e18
                    lopening_costal_solution = []
                    
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
                            cost0 = cal_tour_uncertain(ring_cost, tour)
                            if cost0 < lopening_costal_cost
                                lopening_costal_cost = cost0
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
                            cost0 = cal_tour_certain(ring_cost, tour)
                            if cost0 < lopening_costal_cost
                                lopening_costal_cost = cost0
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
                            cost0 = cal_tour_mixed(ring_cost, tour, H_tilt)
                            if cost0 < lopening_costal_cost
                                lopening_costal_cost = cost0
                            end
                        end
                    end
                    lopening_costal_cost += obj
                    if lopening_costal_cost < global_cost
                        global_cost = lopening_costal_cost
                        global_sol = [1, i2, i3, i4, i5]
                    end
                end
            end
        end
    end
    return global_cost, global_sol, time() - t
end

objval_3, objsol_3, objtime_3 = three_hubs(n, V_certain, ring_cost, star_cost, opening_cost)
objval_4, objsol_4, objtime_4 = four_hubs(n, V_tilt, V_certain, ring_cost, star_cost, opening_cost)
objval_5, objsol_5, objtime_5 = five_hubs(n, V_tilt, ring_cost, opening_cost, star_cost)

global bf_global_obj = objval_5
global bf_global_num_hubs = 5
if objval_3 < bf_global_obj
    bf_global_obj = objval_3
    bf_global_num_hubs = 3
end

if  objval_4 < bf_global_obj
    bf_global_obj = objval_4
    bf_global_num_hubs = 4
end


result_dict = Dict()
result_dict["algorithm"] = "brute_force"
for param_name in ["num_constraint_ilp_including_integrality", 
                    "num_constraint_ilp_notinclude_integrality", 
                    "lower_bound", 
                    "upper_bound", 
                    "num_subtour", 
                    "num_hubs",
                    "time_sp",
                    "time_master",
                    "total_time",
                    "num_cut_sp0",
                    "num_cut_spi",
                    "obj_bf3",
                    "obj_bf4",
                    "obj_bf5",
                    "time_bf3",
                    "time_bf4",
                    "time_bf5",
                    "timestamp"]
    result_dict[param_name] = 0
end

result_dict["obj_bf3"] = objval_3
result_dict["obj_bf4"] = objval_4
result_dict["obj_bf5"] = objval_5

result_dict["time_bf3"] = objtime_3
result_dict["time_bf4"] = objtime_4
result_dict["time_bf5"] = objtime_5

result_dict["timestamp"] = now()
result_dict["total_time"] = objtime_3 + objtime_4 + objtime_5
result_dict["lower_bound"] = bf_global_obj
result_dict["upper_bound"] = bf_global_obj

result_dict["num_hubs"] = bf_global_num_hubs
# println(result_dict) 
write_ouput(pars, name, result_dict, MainPar)


# _write_log_bruteforce(name, n, V_tilt, "brute_force", objval_3, objsol_3, objtime_3, objval_4, objsol_4, objtime_4, objval_5, objsol_5, objtime_5)
