"""
Generates new proposed state from a given microstate by doing the splicing
procedure in Flötteröd & Bierlaire 2013 (TR: Part B).

    splice(g, p_insert, state)
    
    g:          Graph
    p_insert:   Probability distribution over nodes in g
    state:      Micro state (made up of Γ, a, b, c)

Given a micro state (i.e. a Path and anchor points a, c and pivot b) generates
a proposed state (potentially the same) for a transition in the MH algorithm.
"""
function splice(g, p_insert, state)

    insertion_set = setdiff(vertices(g),
                            state.Γ[1:state.a],
                            state.Γ[state.c:end])
    
    p_ν = p_insert[insertion_set]
    ν = sample(insertion_set, p_ν) # insertion node
    
    dist_mat = copy(LightGraphs.weights(g))
    for i in append!(state.Γ[1:state.a-1], state.Γ[state.c:end])
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
    #println("Γ₁: ", Γ₁)
    
    ds = dijkstra_shortest_paths(g, ν, dist_mat)
    paths2 = enumerate_paths(ds, state.Γ[state.c])
    if length(paths2) == 0 
        return state
    end
    Γ₂ = paths2
    #println("Γ₂: ", Γ₂)
    
    if length(intersect(Γ₁, Γ₂)) > 1 # Cycle in new path
        return state
    end
    
    new_Γ = cat(state.Γ[1:state.a-1], Γ₁[1:end-1], Γ₂, state.Γ[state.c+1:end], dims=1)
    new_a = state.a
    new_b = findfirst(isequal(ν), new_Γ)
    new_c = new_b + length(Γ₂) - 1
    
    MicroState(new_Γ, new_a, new_b, new_c)
end

"""
Returns true or false depending if state is `spliceable` as defined in 
Flötteröd & Bierlaire 2013 (TR: Part B)

    is_spliceable(state, g, geodesic_dist_matrix)
    
    state:                  MicroState
    g:                      SimpleWeightedDiGraph
    geodesic_dist_matrix:   Matrix with shortest path distances
"""
function is_spliceable(state, g, geodesic_dist_matrix)

    node_a = state.Γ[state.a]
    node_b = state.Γ[state.b]
    node_c = state.Γ[state.c]
    
    gd1 = geodesic_dist_matrix[node_a, node_b]
    gd2 = geodesic_dist_matrix[node_b, node_c]
    
    Γ₁_length = 0
    for i in state.a:state.b-1
        s, d = state.Γ[i:i+1]
        Γ₁_length += g.weights[s,d]
    end
         
    Γ₂_length = 0    
    for i in state.b:state.c-1
        s, d = state.Γ[i:i+1]
        Γ₂_length += g.weights[s,d]
    end
    
    if Γ₁_length == gd1 && Γ₂_length ==gd2
        return true
    else
        return false
    end
end
