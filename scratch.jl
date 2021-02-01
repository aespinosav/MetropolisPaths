#Adjacency matrix for a braess network
using LightGraphs, SimpleWeightedGraphs, StatsBase

"""
Calculates matrix of shortest path distances
"""
function geod_dist_mat(g)
    d_mat = zeros(nv(g), nv(g))
    for i in vertices(g)
        d_mat[i,:] = gdistances(g, i)
    end
    d_mat
end

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

########################################################################
########################################################################


#Braess network
A_braess = [0 1 1 0 
            0 0 1 1
            0 0 0 1 
            0 0 0 0]
A_braess = Float64.(A_braess)
g = SimpleWeightedDiGraph(A_braess)


state = MicroState([1,2,3,4], 1, 2, 4)

insertion_set = setdiff(vertices(g),
                        state.Γ[1:state.a-1],
                        state.Γ[state.c:end])

pp = StatsBase.Weights(ones(nv(g)))

p_ν = pp[insertion_set]
ν = sample(insertion_set, p_ν) # insertion node

dist_mat_alt = LightGraphs.weights(g)
for i in append!(state.Γ[1:state.a-1], state.Γ[state.c:end])
    for j in outneighbors(g, i)
        global dist_mat_alt[i,j] = Inf
    end
    for j in inneighbors(g, i)
        global dist_mat_alt[j,i] = Inf
    end
end

ds = dijkstra_shortest_paths(g, state.Γ[state.a], dist_mat_alt)
paths1 = enumerate_paths(ds, ν)
Γ₁ = paths1

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


