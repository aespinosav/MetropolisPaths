#Adjacency matrix for a braess network
using LightGraphs, SimpleWeightedGraphs, StatsBase, LinearAlgebra, MHPaths


function geod_dist_mat(g)
    d_mat = zeros(nv(g), nv(g))
    for i in vertices(g)
        d_mat[i,:] = gdistances(g, i)
    end
    d_mat
end


"""
Calculates matrix of shortest path distances.
Uses Floyd Warshal algorithm and needs a WeightedGraph
"""
function geod_dist_mat_fs(g)
    d_mat = zeros(nv(g), nv(g))
    
    fs = floyd_warshall_shortest_paths(g)
    
    for i in vertices(g)
        for j in vertices(g)
            if i != j
                p = enumerate_paths(fs, i, j)
                len = 0
                for k in 1:length(p)-1
                    # g.weights is indexed in reverse [dst,src]
                    len += g.weights[k+1,k]
                end
                d_mat[i,j] = len
            end
        end
    end
    d_mat
end

################################################################

"""
Function for calculating denominator for insert prob for a particular OD pair
for all vertices of a graph
"""
function p_insert_denom(o, d, μ, g, d_mat)
    suma = 0
    for w in vertices(g)
        suma += exp(-μ*(d_mat[o, w] + d_mat[w, d]))
    end
    suma
end

"""
Probability of inserting vertex, normalisation factor should be passed
externally for efficiency
"""
function p_insert(v, o, d, d_mat, μ, denominator)
    exp(-μ*(d_mat[o, v] +  d_mat[v, d]))/denominator
end

"""
Type for representing microstates (of paths).
The graph to which the paths refer is external to the type.
"""
mutable struct MicroState 
    Γ::Array{Int64,1} # Path
    a::Int64 # First acnchor
    b::Int64 # Splice point 
    c::Int64 # Second anchor
end

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


########################################################################
########################################################################


#Braess network
A_braess = [0 1 1 0 
            0 0 1 1
            0 0 0 1 
            0 0 0 0]
A_braess = Float64.(A_braess)
g = SimpleWeightedDiGraph(A_braess)

pp = StatsBase.Weights(ones(nv(g)))

state = MicroState([1,2,3,4], 1, 2, 4)

splice(g, pp, state)



#Geodesic distance matrix for braess network
d_mat = geod_dist_mat(g)

# Inverse tempr for metropolis hastings
μ = 0.5
o = 1
d = 4

# P insert probability
#Denominator of p
denominator = p_insert_denom(1, 4, μ, g, d_mat)

p_insert_dist = zeros(nv(g))
for v in vertices(g)
   global p_insert_dist[v] = p_insert(v, o, d, d_mat, μ, denominator)
end
# Use of global needed for REPL in v1.4 (careful with scoping)


############################################################################
############################################################################


#Braess network
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

p_insert = make_p_insert_with_denom(o, d, g, μ, d_mat)
p_dist = StatsBase.Weights([p_insert(v) for v in vertices(g)])

state = MicroState([1,2,3,4], 1, 2, 4)

splice(g, p_dist, state)


###########################################################################
###########################################################################
using LightGraphs, SimpleWeightedGraphs, StatsBase, MHPaths
#include("metropolis.jl")

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
