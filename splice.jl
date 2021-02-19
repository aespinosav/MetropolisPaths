# TODO: Splicing seems to be returning paths with repeated edges. This must be urgently fixed. At the moment it seems to only happen when the protected  path segments disconnect the spliced vertex from the anchors. Also, for debugging purposes it seems smarter to have a splicing routine separate from the drawing routine.

"""
Generates new proposed state from a given microstate by doing the splicing
procedure in Flötteröd & Bierlaire 2013 (TR: Part B).

    splice(state, g, p_insert)
    
    state:      Micro state (made up of Γ, a, b, c)
    g:          Graph
    p_insert:   Probability distribution over nodes in g

Given a micro state (i.e. a Path and anchor points a, c and pivot b) generates
a proposed state (potentially the same) for a transition in the MH algorithm.
"""
function splice(state::MicroState, g, p_insert)

    insertion_set = setdiff(vertices(g),
                            state.Γ[1:state.a],
                            state.Γ[state.c:end])
    
    p_ν = p_insert[insertion_set]
    ν = sample(insertion_set, p_ν) # insertion node
    
   
    #TODO this should be checked very well...
     # Distance matrix with reduced node set (set Inf distance)
    dist_mat = copy(LightGraphs.weights(g))
    excluded_nodes = cat(state.Γ[1:state.a-1], state.Γ[state.c:end], dims=1)
    
    for i in excluded_nodes
        for j in outneighbors(g, i)
            dist_mat[i,j] = Inf
        end
        for j in inneighbors(g, i)
            dist_mat[j,i] = Inf
        end
    end
    
    ds = dijkstra_shortest_paths(g, state.Γ[state.a], dist_mat)
    paths1 = enumerate_paths(ds, ν)
    if length(paths1) == 0
        return state
    end
    Γ₁ = paths1
    
    ds = dijkstra_shortest_paths(g, ν, dist_mat)
    paths2 = enumerate_paths(ds, state.Γ[state.c])
    if length(paths2) == 0 
        return state
    end
    Γ₂ = paths2
    
    if length(intersect(Γ₁, Γ₂)) > 1 # Cycle in new path
        return state
    end
    
    new_Γ = cat(state.Γ[1:state.a-1], Γ₁[1:end-1], Γ₂, state.Γ[state.c+1:end], dims=1)
    new_a = state.a
    new_b = findfirst(isequal(ν), new_Γ)
    new_c = new_b + length(Γ₂) - 1
    
    MicroState(new_Γ, new_a, new_b, new_c)
end

splice(state::MicroState, mh::MHInstance) = splice(state, mh.g, mh.p_insert)


######################################################################

function draw_insertion_vertex(state::MicroState, g, p_insert)
    insertion_set = setdiff(vertices(g),
                            state.Γ[1:state.a],
                            state.Γ[state.c:end])
    
    p_ν = p_insert[insertion_set]
    ν = sample(insertion_set, p_ν)
end


function splice_sparse(state::MicroState, g::AbstractGraph, ν::Int)

    # 1 reduced node set for each splice segment
    dist_mat_1 = adjacency_matrix(g)
    excluded_nodes_1 = cat(state.Γ[1:state.a-1], state.Γ[state.c:end], dims=1)
    for i in excluded_nodes_1
        for j in outneighbors(g, i)
            dist_mat_1[i,j] = Inf
        end
        for j in inneighbors(g, i)
            dist_mat_1[j,i] = Inf
        end
    end
    
    # Having 2 matrices defined here seems very memory inefficient
    dist_mat_2 = adjacency_matrix(g)
    excluded_nodes_2 = cat(state.Γ[1:state.a], state.Γ[state.c+1:end], dims=1)
    for i in excluded_nodes_2
        for j in outneighbors(g, i)
            dist_mat_2[i,j] = Inf
        end
        for j in inneighbors(g, i)
            dist_mat_2[j,i] = Inf
        end
    end
    
    ds = dijkstra_shortest_paths(g, state.Γ[state.a], dist_mat_1)
    paths1 = enumerate_paths(ds, ν)
    Γ₁ = paths1
    
    # Check no disallowed edges are used
    suma = 0
    for i in 1:length(Γ₁)
        suma += dist_mat_1[i,i+1]
    end
    #@show suma
    if suma == Inf
        return state
    end
    
    ds = dijkstra_shortest_paths(g, ν, dist_mat_2)
    paths2 = enumerate_paths(ds, state.Γ[state.c])
    Γ₂ = paths2
    
    # Check no disallowed edges are used
    suma = 0
    for i in 1:length(Γ₂)
        suma += dist_mat_1[i,i+1]
    end
    #@show suma
    if suma == Inf
        return state
    end
    
    # Check no cycle has been created
    if length(intersect(Γ₁, Γ₂)) > 1 # only intersect at ν
        return state
    end
    
    new_Γ = cat(state.Γ[1:state.a-1], Γ₁[1:end-1], Γ₂, state.Γ[state.c+1:end], dims=1)
    new_a = state.a
    new_b = findfirst(isequal(ν), new_Γ)
    new_c = new_b + length(Γ₂) - 1
    MicroState(new_Γ, new_a, new_b, new_c)
end








"""
Same as splice but uses a sparse dist_mat for dijkstra shortest paths
"""
function splice_sp(state::MicroState, g, p_insert)

    insertion_set = setdiff(vertices(g),
                            state.Γ[1:state.a],
                            state.Γ[state.c:end])
    
    p_ν = p_insert[insertion_set]
    ν = sample(insertion_set, p_ν) # insertion node
    
   
    #TODO this should be checked very well...
     # Distance matrix with reduced node set (set Inf distance)
    #dist_mat = copy(LightGraphs.weights(g))
    dist_mat = adjacency_matrix(g)
    excluded_nodes = cat(state.Γ[1:state.a-1], state.Γ[state.c:end], dims=1)
    
    for i in excluded_nodes
        for j in outneighbors(g, i)
            dist_mat[i,j] = Inf
        end
        for j in inneighbors(g, i)
            dist_mat[j,i] = Inf
        end
    end
    
    ds = dijkstra_shortest_paths(g, state.Γ[state.a], dist_mat)
    paths1 = enumerate_paths(ds, ν)
    if length(paths1) == 0
        # I guess that this should actually return an error?
        return state
    end
    Γ₁ = paths1
    
    ds = dijkstra_shortest_paths(g, ν, dist_mat)
    paths2 = enumerate_paths(ds, state.Γ[state.c])
    if length(paths2) == 0 
        return state
    end
    Γ₂ = paths2
    
    if length(intersect(Γ₁, Γ₂)) > 1 # Cycle in new path
        return state
    end
    
    new_Γ = cat(state.Γ[1:state.a-1], Γ₁[1:end-1], Γ₂, state.Γ[state.c+1:end], dims=1)
    new_a = state.a
    new_b = findfirst(isequal(ν), new_Γ)
    new_c = new_b + length(Γ₂) - 1
    
    MicroState(new_Γ, new_a, new_b, new_c)
end


######################################################################


"""
Returns true or false depending if state is `spliceable` as defined in 
Flötteröd & Bierlaire 2013 (TR: Part B)

    is_spliceable(state, g, geodesic_dist_matrix)
    
    state:                  MicroState
    g:                      SimpleWeightedDiGraph
    geodesic_dist_matrix:   Matrix with shortest path distances
"""
function is_spliceable(state, g, geodesic_dist_matrix)

    a = state.a
    b = state.b
    c = state.c

    node_a = state.Γ[a]
    node_b = state.Γ[b]
    node_c = state.Γ[c]
    
    gd1 = geodesic_dist_matrix[node_a, node_b]
    gd2 = geodesic_dist_matrix[node_b, node_c]
    
    Γ₁ = state.Γ[a:b]
    Γ₂ = state.Γ[b:c]
    
    Γ₁_length = path_length(Γ₁, g)
    Γ₂_length = path_length(Γ₂, g)

    if Γ₁_length≈gd1 && Γ₂_length≈gd2
        return true
    else
        return false
    end
end

is_spliceable(state, mh::MHInstance) = is_spliceable(state, mh.g, mh.sp_dist_mat)
