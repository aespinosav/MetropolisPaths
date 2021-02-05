# Test 1: Make sure it gives an answer
A_braess = [0 1 1 0 
            0 0 1 1
            0 0 0 1 
            0 0 0 0]
A_braess = Float64.(A_braess)
g = SimpleWeightedDiGraph(A_braess)
d_mat = geod_dist_mat(g)

p_ins = StatsBase.Weights(ones(nv(g)))
p_spl = 0.3

state1 = MicroState([1,2,3,4], 1, 2, 4)
state2 = splice(g, p_ins, state1)

state3 = MicroState([1,3,4], 1, 2, 3)
state4 = MicroState([1,2,4], 1, 2, 3)

# state 1 is not slpiceable
proposal_probability(state1, state2, p_spl, p_ins, g, d_mat)

# state 2 is spliceable
proposal_probability(state2, state1, p_spl, p_ins, g, d_mat)

# state 3 and 4 are both spliceable
proposal_probability(state3, state4, p_spl, p_ins, g, d_mat)

# Same state
proposal_probability(state2, state3, p_spl, p_ins, g, d_mat)

# Should actually check this properly!!


###################################################################
###################################################################


# Test 2: Larger Network

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

s = mh.history[end]
s2 = splice(g, p_ins, s)
is_spliceable(s, g, mh.sp_dist_mat)
is_spliceable(s2, g, mh.sp_dist_mat)


