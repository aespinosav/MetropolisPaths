using LightGraphs, SimpleWeightedGraphs, StatsBase, MHPaths

# Test 1: Check that it is constructing properly
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
p_splice = 0.5

MHinstance(g, o, d, μ, p_splice)
