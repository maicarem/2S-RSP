# using GZip
# using LinearAlgebra

# function read_dat(file_name)
#     file = GZip.open(file_name)
#     n = 51
#     x = zeros(Int, n)
#     y = zeros(Int, n)

#     for _ in 1:6
#         readline(file)
#     end
    
#     for i in 1:n
#         _, x[i], y[i] = [parse(Int, i) for i in split(readline(file)," ")]
#     end
#     return x, y
# end

# function create_dat(x, y)
#     n = length(x)
    
#     star_cost = zeros(Float32, n, n)
    
#     for i in 1:n
#         for j in 1:n
#             if i!=j
#                 star_cost[i,j] = norm([x[i] - x[j], y[i] - y[j]], 2)
#             end
#         end
#     end
    
#     ring_cost = star_cost*(10-2)
#     opening_cost = zeros(Float32,n)
#     return opening_cost, star_cost, ring_cost

# end


# x, y = read_dat("ALL_tsp/eil51.tsp.gz")
# @show create_dat(x, y)
 
