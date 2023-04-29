function create_connectivity_cut_strategy_1(cb_data, x, y, V, n, pars, benders_enabled)

    ε = 1e-6
    x_m = zeros(Float64, n, n)
    if benders_enabled
        y_m = Float64[callback_value(cb_data, y[i]) for i in 1:n]
    else
        y_m = Float64[callback_value(cb_data, y[i,i]) for i in 1:n]
    end
    flow_graph = Graphs.DiGraph(n)
    capacity_matrix = zeros(Float64, n, n)
    for i in 1:n
        for j in i+1:n
            x_m[i,j] = callback_value(cb_data, x[i, j])
            if x_m[i,j] > ε
                Graphs.add_edge!(flow_graph, i, j) ; Graphs.add_edge!(flow_graph, j, i)
                capacity_matrix[i,j] = x_m[i,j]
                capacity_matrix[j,i] = x_m[i,j]
            end
        end
    end
    max_violated_node = -1
    max_current_violation = .0
    bestf, bestpart1, bestpart2 = .0, Int[], Int[]
    for i in 2:n
        if y_m[i] > ε
            part1, part2, flow = GraphsFlows.mincut(flow_graph, 1, i, capacity_matrix, EdmondsKarpAlgorithm())
            if max_current_violation < 2y_m[i] - flow
                max_violated_node = i
                max_current_violation = 2y_m[i] - flow
                bestf = flow
                bestpart1 = part1
                bestpart2 = part2
            end
        end
    end

    if max_current_violation > pars.uc_tolerance
        if benders_enabled
            con = @build_constraint(sum(sum(x[minmax(j,k)] for j in bestpart1) for k in bestpart2) >= 2y[max_violated_node])
        else
            con = @build_constraint(sum(sum(x[minmax(j,k)] for j in bestpart1) for k in bestpart2) >= 2y[max_violated_node,max_violated_node])
        end
        @info "User cut: $(con)"
        return max_current_violation, con
    end
    return -1, []
end
function create_connectivity_cut_strategy_2(cb_data, x, y, V, n, pars, benders_enabled)

    ε = 10e-16
    x_m = zeros(Float64, n, n)
    if benders_enabled
        y_m = Float64[callback_value(cb_data, y[i]) for i in 1:n]
    else
        y_m = Float64[callback_value(cb_data, y[i,i]) for i in 1:n]
    end
    flow_graph = Graphs.DiGraph(n)
    capacity_matrix = zeros(Float64, n, n)
    for i in 1:n
        for j in i+1:n
            x_m[i,j] = callback_value(cb_data, x[i, j])
            if x_m[i,j] > ε
                Graphs.add_edge!(flow_graph, i, j) ; Graphs.add_edge!(flow_graph, j, i)
                capacity_matrix[i,j] = x_m[i,j]
                capacity_matrix[j,i] = x_m[i,j]
            end
        end
    end
    max_violated_node = -1
    max_current_violation = .0

    bestf, bestpart1, bestpart2 = .0, Int[], Int[]
    for i in 2:n
        if y_m[i] > ε
            part1, part2, flow = GraphsFlows.mincut(flow_graph, 1, i, capacity_matrix, EdmondsKarpAlgorithm())
            if max_current_violation < 2y_m[i] - flow
                max_violated_node = i
                max_current_violation = 2y_m[i] - flow
                bestf = flow
                bestpart1 = part1
                bestpart2 = part2
                break
            end
        end
    end

    if max_current_violation > pars.uc_tolerance
        if benders_enabled
            con = @build_constraint(sum(sum(x[minmax(j,k)] for j in bestpart1) for k in bestpart2) >= 2y[max_violated_node])
        else
            con = @build_constraint(sum(sum(x[minmax(j,k)] for j in bestpart1) for k in bestpart2) >= 2y[max_violated_node,max_violated_node])
        end
        return max_current_violation, con
    end
    return -1, []
end
function create_connectivity_cut_strategy_3(cb_data, x, y, V, n, pars, benders_enabled)
    ε = 10e-16
    x_m = zeros(Float64, n, n)
    if benders_enabled
        y_m = Float64[callback_value(cb_data, y[i]) for i in 1:n]
    else
        y_m = Float64[callback_value(cb_data, y[i,i]) for i in 1:n]
    end
    flow_graph = Graphs.DiGraph(n)
    capacity_matrix = zeros(Float64, n, n)
    for i in 1:n
        for j in i+1:n
            x_m[i,j] = callback_value(cb_data, x[i, j])
            if x_m[i,j] > ε
                Graphs.add_edge!(flow_graph, i, j) ; Graphs.add_edge!(flow_graph, j, i)
                capacity_matrix[i,j] = x_m[i,j]
                capacity_matrix[j,i] = x_m[i,j]
            end
        end
    end
    max_violated_node = -1
    max_current_violation = .0

    bestf, bestF, bestlabels = .0, zeros(Float64, n, n), ones(Int64, n)
    for i in 2:n
        if y_m[i] > ε
            flow, F, labels = maximum_flow(flow_graph, 1, i, capacity_matrix, algorithm=BoykovKolmogorovAlgorithm())
            if max_current_violation < 2y_m[i] - flow
                if sum(labels .== 1) <= sum(bestlabels .== 1)
                    max_violated_node = i
                    max_current_violation = 2y_m[i] - flow
                    bestf = flow
                    bestF = F
                    bestlabels = labels
                end
            end
        end
    end
    if max_current_violation > pars.uc_tolerance
        if benders_enabled
            con = @build_constraint(sum(sum(x[j,k] for j in V
            if j < k && ((bestlabels[j] == 1 && bestlabels[k] != 1) ||
            (bestlabels[k] == 1 && bestlabels[j] != 1))) for k in V) >= 2y[max_violated_node])
        else
            con = @build_constraint(sum(sum(x[j,k] for j in V
            if j < k && ((bestlabels[j] == 1 && bestlabels[k] != 1) ||
            (bestlabels[k] == 1 && bestlabels[j] != 1))) for k in V) >= 2y[max_violated_node, max_violated_node])
        end
        
        return max_current_violation, con
    end
    return -1, []
end
function create_connectivity_cut_strategy_4(cb_data, x, y, V, n, nconnectivity_cuts, pars, benders_enabled)
    if nconnectivity_cuts < pars.uc_strat_4_limit
        return create_connectivity_cut_strategy_1(cb_data, x, y, V, n, pars, benders_enabled)
    end
    return -1, []
end

function call_back_user_cuts(cb_data)
    if pars.uc_strat == 1
        max_current_value, con = create_connectivity_cut_strategy_1(cb_data, x, y, V, n, pars, benders_enabled)
        if max_current_value > pars.uc_tolerance
            MOI.submit(main, MOI.UserCut(cb_data), con)
            # nconnectivity_cuts += 1
        end
    elseif pars.uc_strat == 2
        max_current_value, con = create_connectivity_cut_strategy_2(cb_data, x, y, V, n, pars, benders_enabled)
        if max_current_value > pars.uc_tolerance
            MOI.submit(main, MOI.UserCut(cb_data), con)
            # nconnectivity_cuts += 1
        end
    elseif pars.uc_strat == 3
        max_current_value, con = create_connectivity_cut_strategy_3(cb_data, x, y, V, n, pars, benders_enabled)
        if max_current_value > pars.uc_tolerance
            MOI.submit(main, MOI.UserCut(cb_data), con)
            # nconnectivity_cuts += 1
        end
    elseif pars.uc_strat == 4
        max_current_value, con = create_connectivity_cut_strategy_4(cb_data, x, y, V, n, nconnectivity_cuts, pars, benders_enabled)
        if max_current_value > pars.uc_tolerance
            MOI.submit(main, MOI.UserCut(cb_data), con)
            # nconnectivity_cuts += 1
        end
    end
end