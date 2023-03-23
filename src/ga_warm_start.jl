include("dat.jl")
include("dual_solution.jl")

function initialize(size)
    return rand(Float64, size)
end

function find_hub(encoding)
    
    chromosome = copy(encoding)
    encoding .= ifelse.(encoding .> bench, 1, 0)
    chromosome .= ifelse.(chromosome .> bench, chromosome, 0)
    
    while sum(encoding) < 5 || encoding[1] == 0
        encoding = rand(Float64, n)
        chromosome = copy(encoding)

        encoding .= ifelse.(encoding .> bench, 1, 0)
        chromosome .= ifelse.(chromosome .> bench, chromosome, 0)

    end
    
    return chromosome
end

# output the ring [1,2,3,5,3,1]
function init_ring(chromosome)
    encoding = copy(chromosome)
    encoding .= ifelse.(encoding .> bench, 1, 0)
    hub = [i for i in 1: length(encoding) if encoding[i]==1]
    encoding1 = [chromosome[i] for i in 1: length(encoding) if encoding[i] == 1]
    return append!(hub[sortperm(encoding1)], hub[sortperm(encoding1)][1])
end

# output: x_ij, i<j
function decoding_ring(_hub)
    _x = zeros(Int, n,n)
    idx = 1
    while idx + 1 <= length(_ring)
        _x[_hub[idx], _hub[idx+1]] = 1
        _x[_hub[idx+1], _hub[idx]] = 1
        idx +=1
    end
    return _x
end

function two_optSwap(route, v1::Int, v2::Int)
    new_route = zeros(Int, length(route))
    # take route[0] to route[v1] and add them in order to new_route
    new_route[1:v1] = route[1:v1]
    # take route[v1+1] to route[v2] and add them in reverse order to new_route
    new_route[v1+1:v2] = reverse(route[v1+1:v2])
    new_route[v2+1:length(route)] = route[v2+1:length(route)]
    return new_route
end

function cal_obj_sp0(alpha, beta, x_hat)
    return sum((x_hat[k,i] + x_hat[i,j]+x_hat[j,t] - 2) * alpha[i,j,k,t] for (i,j,k,t) in J_tilt) + sum((x_hat[i,j]+ x_hat[j,k] - 1)* beta[i,j,k] for (i,j,k) in K_tilt)
end

function cal_obj_spi(varphi, gamma, y_hat)
    obj = zeros(Float64, n)
    for i in 1:n
        y_hat[i,i] == 0 || continue
        obj[i] = 3(1-y_hat[i,i])*varphi[i] + sum([y_hat[j,j]*gamma[i,j] for j in V if j!=i])
    end
    return sum(obj)
end

function cal_backup_and_star(route)
    
    y_hat = zeros(n,n)
    x_hat = zeros(n,n)
    length_route = length(route)
    
    for i in 1:length_route-1
        x_hat[route[i],route[i+1]] = 1
        x_hat[route[i+1],route[i]] = 1
        y_hat[route[i],route[i]] = 1
    end

    obj_0 = sum(rc[i,j]*x_hat[i,j] for (i,j) in E)+ sum(oc[i]*y_hat[i,i] for i in V)
    (beta, alpha), (φ, γ) = dual_solution(y_hat, x_hat)
        
    obj_sp0 = cal_obj_sp0(alpha, beta, x_hat)
    obj_spi = cal_obj_spi(φ, γ, y_hat)
    return obj_0 + obj_sp0 + obj_spi
end

function cal_route(route)
    return offset + cal_backup_and_star(route)
end

function cal_objective_val(chromosome)
    route = init_ring(chromosome)
    return cal_route(route)
end

# lengthDelta = - dist(route[v1], route[v1+1]) - dist(route[v2], route[v2+1]) + dist(route[v1+1], route[v2+1]) + dist(route[v1], route[v2])
function length_delta(route, v1, v2)
    return -rc[route[v1], route[v1+1]] - rc[route[v2], route[v2+1]] + rc[route[v1+1], route[v2+1]] + rc[route[v1], route[v2]]
end

function two_opt_procedure(route)
    best_obj = cal_route(route)
    for i in 1:length(route)-1
        for j in i+1:length(route)-1
            length_delta(route, i, j) <0 || continue
            best_obj += length_delta(route, i, j)
            route = two_optSwap(route,i,j)
        end
    end
    return route
end

function initialize_population(size_pop)
    popu = []
    pop_obj = ones(Float64, size_pop) * 1e18
    for i in 1: size_pop
        _init0 = initialize(n)
        chromosome = find_hub(_init0)
        push!(popu, chromosome)
        pop_obj[i] = cal_route(init_ring(chromosome))
    end
    return popu, pop_obj
end

function mutation(child)
    hub_list = init_ring(child)
    v1, v2 = rand(hub_list[2:end], 2)
    
    child_orin = copy(child)
    a = child[v2]
    child[v2] = child[v1]
    child[v1] = a
    
    enc1 = copy(child)
    enc1 .= ifelse.(enc1 .> bench, 1, 0)
    
    # while sum(enc1) < 3
    #     v1, v2 = rand(2:n, 2)

    #     a = child_orin[v2]
    #     child[v2] = child_orin[v1]
    #     child[v1] = a
        
    #     enc1 = copy(child)
    #     enc1 .= ifelse.(enc1 .> bench, 1, 0)
    # end

    return child
end

function crossover(child1, child2)
    position = rand(2:n-1)
    c1 = append!(child1[1:position],child2[position+1:n])
    c2 = append!(child2[1:position],child1[position+1:n])
    return c1, c2
end

function tournament_selection(population, obj_popu, size_pop)
    selection_idx = rand(1:size_pop)
    for i in rand(1:size_pop,2)
        obj_popu[i] < obj_popu[selection_idx] || continue
        selection_idx = i
    end
    return population[selection_idx], obj_popu[selection_idx]
end

function genetic_algorithm(size_pop, size_iter)
    
    popu, obj_pop = initialize_population(size_pop)
    # idx_smallest = argmin(obj_pop)
    # best_chromosome, best_obj_chromosome = popu[idx_smallest], obj_pop[idx_smallest]

    for _ in 1:size_iter
        # tournament selection
        for i in 1:size_pop
            popu[i], obj_pop[i] = tournament_selection(popu, obj_pop, size_pop)
        end

        for i in 1:size_pop
            # mutation
            if rand() > 0.3
                popu[i] = mutation(popu[i])
                obj_pop[i] = cal_objective_val(popu[i])
            end
            
            # crossover
            # if rand() > 0.6
            #     j = rand(1:size_pop)
            #     while j == i
            #         j = rand(1:size_pop)
            #     end
            
            #     popu[i], popu[j] = crossover(popu[i], popu[j])
            #     obj_pop[i], obj_pop[j] = cal_objective_val(popu[i]), cal_objective_val(popu[j])
            
            # end

            # 2-opt
            if rand() > 0.1
                # @info "Use two-opt"
                _ring = init_ring(popu[i])
                _length_ring = length(_ring)
                
                _new_ring= two_opt_procedure(_ring)
                obj_pop[i] = cal_route(_new_ring)

                for idx in 1:_length_ring
                    popu[i][_ring[idx]] = popu[i][_new_ring[idx]]
                end
            end
        end
    end
    min0 = argmin(obj_pop)
    return init_ring(popu[min0]), obj_pop[min0]
end

@show genetic_algorithm(20,20)
# route = [5,1,2,3,4,5]
# @show cal_route(route)