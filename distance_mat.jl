"""
Calculates matrix of shortest path distances from a Floyd-Warshal state
"""
function geod_dist_mat_fs_in(g, fs)
    d_mat = zeros(nv(g), nv(g))
    
    for i in vertices(g)
        for j in vertices(g)
            if i != j
                p = enumerate_paths(fs, i, j)
                len = 0
                for k in 1:length(p)-1
                    # g.weights is indexed in reverse [dst,src]
                    src, dst = p[k], p[k+1]
                    len += g.weights[dst, src]
                end
                d_mat[i,j] = len
            end
        end
    end
    d_mat
end


"""
Calculates matrix of shortest path distances.
Uses Floyd Warshal algorithm and needs a WeightedGraph
"""
function geod_dist_mat(g)
    fs = floyd_warshall_shortest_paths(g)
    geod_dist_mat_fs_in(g, fs)
end


"""
Returns path lenght by adding length of edges in path
when given an array of nodes in the path (no multi-edges)

    path_length(path, g)
    
    path:   Array of nodes in path
    g:      WeightedGraph
    
Remember g.weights is indexed weirdly [dst, src]
"""
function path_length(path, g::SimpleWeightedDiGraph)
    num_edges = length(path)-1
    suma = 0
    for i in 1:num_edges
        s, d = path[i], path[i+1]
        suma += g.weights[d, s]
    end
    suma
end

path_length(path, mh::MHInstance) = path_length(path, mh.g)
path_length(state::MicroState, mh::MHInstance) = path_length(state.Γ, mh.g)


