# Test 1: Check that it runs
A_braess = [0 1 1 0 
            0 0 1 1
            0 0 0 1 
            0 0 0 0]
A_braess = Float64.(A_braess)
g = SimpleWeightedDiGraph(A_braess)

p_ins = StatsBase.Weights(ones(nv(g)))

#Geodesic distance matrix for braess network
d_mat = geod_dist_mat(g)

# Inverse tempr for metropolis hastings
μ = 0.5
o = 1
d = 4
p_splice = 0.3

mh = MHInstance(g, o, d, μ, p_splice)
mh.history[1] = state = MicroState([1,2,3,4], 1, 2, 4)
#println("Init history: ", mh.history)
# Sample

num_steps = 10

mh_evolve!(mh, num_steps, weight_func)
#println("Final history: ", mh.history)


##########################################################################
##########################################################################



# Test 2: Check that it does things properly
points = [3 0;
          2 2;
          2 4;
          5 4;
          5 2;
          8 3]
points = Float64.(points)          
          
edges_array = [(1, 2),
         (1, 5),
         (2, 3),
         (2, 5),
         (3, 4),
         (4, 5),
         (4, 6),
         (5, 6)]
         
edge_lens = [norm(points[edges_array[i][2],:] - points[edges_array[i][1],:]) for i in 1:length(edges_array)]

g = SimpleWeightedDiGraph(size(points)[1])
let g = g       
    for i in 1:length(edges_array)
        add_edge!(g, edges_array[i]..., edge_lens[i])
        add_edge!(g, reverse(edges_array[i])..., edge_lens[i])
    end
end

p_ins = StatsBase.Weights(ones(nv(g)))

d_mat = geod_dist_mat(g)

# Inverse tempr for metropolis hastings
μ = 0.5
o = 1
d = 4
p_splice = 0.3

mh = MHInstance(g, o, d, μ, p_splice)

num_steps = 10

mh_evolve!(mh, num_steps, weight_func)
