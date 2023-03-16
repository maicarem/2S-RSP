n = 6
sc = ones(n,n)
# sc = [0 7 4 2 1 9 1 3 8 2;
#     5 0 5 9 6 7 6 7 6 4; 
#     8 5 0 3 10 3 5 9 1 1; 
#     8 9 2 0 6 7 1 6 10 5; 
#     10 6 8 4 0 3 1 10 2 9; 
#     3 9 3 4 9 0 5 7 6 2; 
#     1 7 6 7 1 4 0 9 10 4; 
#     2 2 10 9 6 9 10 0 7 5; 
#     2 2 3 8 6 5 7 10 0 5; 
#     2 1 2 1 7 6 3 9 9 0]
# sc = ones(n,n)
# sc = rand(1:10, n, n)
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