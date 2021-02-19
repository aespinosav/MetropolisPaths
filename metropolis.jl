import Base.show

"""
Object for running path sampling
"""
mutable struct MHInstance
    g # Graph
    o # Origin node
    d # Destination node
    sp_dist_mat # shortest path distance matrix
    μ # Inverse temp parameter
    p_splice # Probability of splice (constant)
    p_insert # Vector of probs of inserting a given node in splice
    history # Array of microstates
end
"""
Constructor for MHInstance that calculates the p_insert probabilities and the 
shortest-path distance matrix.

    MHinstance(g, o, d, μ, p_splice)
    
    g:          Graph (SimpleWeightedDiGraph?)
    o:          Origin node
    d:          Destination node
    μ:          Metropolis-Hastings parameter
    p_splice:   Probability of splicing (Constant)
"""
function MHInstance(g, o, d, μ, p_splice)

    N = length(vertices(g))

    fs = floyd_warshall_shortest_paths(g)

    geodesic_distance_matrix = geod_dist_mat_fs_in(g, fs)

    initial_path = enumerate_paths(fs, o, d)

    # Initial sampling weights uniform for the first microstate
    n_in_p = length(initial_path)
    
    weights_init = StatsBase.Weights(ones(n_in_p))
                 
    anchors = sample(1:n_in_p, 
                     weights_init,
                     3,
                     replace=false,
                     ordered=true)     
                                 
    # Initial state is shortest path with random anchors
    initial_state = MicroState(initial_path, anchors...)
    
    history = [initial_state]
    
    p_ins_func = make_p_insert_with_denom(o, d, g, μ, geodesic_distance_matrix)
    p_insert = zeros(N)
    for i in 1:N
        p_insert[i] = p_ins_func(i) 
    end
    p_insert = StatsBase.Weights(p_insert)
    
    MHInstance(g, o, d, geodesic_distance_matrix, μ, p_splice, p_insert, history)
end

function show(io::IO, mh::MHInstance)
    println(io, "Metropolis-Hastings instance:")
    showstr = "g:\t$(mh.g)\no:\t$(mh.o)\nd:\t$(mh.d)\nμ:\t$(mh.μ)\np_spl:\t$(mh.p_splice)"
    print(io, showstr)
end


"""
Get unique paths from the history of an MHInstance
    
    get_unique_paths(mh::MHInstance)
"""
function get_unique_paths(mh::MHInstance)
    paths = Array{Int64,1}[]
    for s in mh.history
        push!(paths, s.Γ)
    end
    unique_paths = unique(paths)
end
