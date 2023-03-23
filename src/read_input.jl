using LinearAlgebra

function read_input(filepath)
    
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
    ring_cost = dist * 2
    readline(f)
    opening_cost= [parse(Float64, i) for i in split(readline(f)," ")]
    
    close(f)
    return _n, opening_cost, dist, ring_cost # number of nodes, opening cost, star cost, ring cost
end

# @show read_input("instances/small_instances/small_instance_20.dat")

