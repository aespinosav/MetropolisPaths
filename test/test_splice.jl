using LinearAlgebra, LightGrqaphs, SimpleWeightedGraph

# Test 1: Braess directed graph (uniform p_insert)

A_braess = [0 1 1 0 
            0 0 1 1
            0 0 0 1 
            0 0 0 0]
A_braess = Float64.(A_braess)
g = SimpleWeightedDiGraph(A_braess)
p_in = StatsBase.Weights(ones(nv(g)))
state = MicroState([1,2,3,4], 1, 2, 4)
splice(g, p_in, state)

########################################################


# Test 2: Braess directed graph (p_insert from Flötteröd & Bierlaire paper)

A_braess = [0 1 1 0 
            0 0 1 1
            0 0 0 1 
            0 0 0 0]
A_braess = Float64.(A_braess)
g = SimpleWeightedDiGraph(A_braess)
d_mat = geod_dist_mat(g)

o = 1
d = 4

μ = 0.5

p_insert = make_p_insert_with_denom(o, d, g, μ, d_mat)
p_dist = StatsBase.Weights([p_insert(v) for v in vertices(g)])

state = MicroState([1,2,3,4], 1, 2, 4)

splice(g, p_dist, state)

###############################################################


# Test 3: Larger network
points = [3 0;
          2 2;
          2 4;
          5 4;
          5 2;
          8 3]
points = Float64.(points)          
          
edges = [(1, 2),
         (1, 5),
         (2, 3),
         (2, 5),
         (3, 4),
         (4, 5),
         (4, 6),
         (5, 6)]
         
edge_lens = [norm(points[edges[i][2],:] - points[edges[i][1],:]) for i in 1:length(edges)]

g = SimpleWeightedDiGraph(size(points)[1])
let g = g       
    for i in 1:length(edges)
        add_edge!(g, edges[i]..., edge_lens[i])
        add_edge!(g, reverse(edges[i])..., edge_lens[i])
    end
end
# Seems like the edge ordering is done sorting by origin first then destination
