using MHPaths, LightGraphs, SkeletonCities2, SimpleWeightedGraphs, MetaGraphs, StatsBase

# αβ-network parameters
n_root = 15
α = 0.2
β = 1.3

# MH parameters
μ = 0
p_splice = 0.8

num_metropolis_steps = 1000

###
### Generate network and put into correct format
###

gg = αβ_network_meta(n_root, α, β)

weighted_adj = Float64.(copy(adjacency_matrix(gg)))
for edge in edges(gg)
    w = get_prop(gg, edge, :length)
    weighted_adj[edge.src, edge.dst] = w
    weighted_adj[edge.dst, edge.src] = w
end

g = SimpleWeightedDiGraph(weighted_adj)

# Uniform sampling density
p_ins = Weights(ones(nv(g)))

# Random origin and destination
#o = sample(1:nv(g), p_ins)
#d = sample(1:nv(g), p_ins)
#while d == o
#    # vertices(g) returns a range type
#    d = sample(vertices(g), p_ins)
#end

o = 1
d = 225

###
### Set up MH experiment
###

mh = MHInstance(g, o, d, μ, p_splice)
mh.p_insert = p_ins # Uniform probability of inserting nodes

# standard weight_func from 
mh_evolve!(mh, num_metropolis_steps, weight_func)

###
### Post-processing
###

function get_unique_paths(mh::MHInstance)
    paths = Array{Int64,1}[]
    for s in mh.history
        push!(paths, s.Γ)
    end
    unique_paths = unique(paths)
end

u_paths = get_unique_paths(mh)
num_u_paths = length(u_paths)
println("Number of unique paths: ", num_u_paths, "\n")

node_coords = get_node_pos(gg)
save_graph_tikz(g, node_coords, "skeleton_net.tex")
save_paths_tikz(g, u_paths, node_coords, "skeleton_metropolis_test.tex")

save_paths_tikz_exp(g,
                    u_paths,
                    node_coords,
                    "od_skeleton_metropolis_test.tex",
                    imp_nodes=[o,d])
