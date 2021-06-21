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
function p_insert_explicit(v, o, d, μ, d_mat, denominator)
    exp(-μ*(d_mat[o, v] +  d_mat[v, d]))/denominator
end

# TODO: Make sure this is working properly!
"""
Returns function where denominator does not need to be passed explicitly.
Check that that this μ is different than the one for the weight function used to sample paths
"""
function make_p_insert_with_denom(o, d, g, μ, d_mat)
    den = p_insert_denom(o, d, μ, g, d_mat)
    v -> p_insert_explicit(v, o, d, μ, d_mat, den)
end

"""
Calculates p_insert (Fötteröd & Bierlaire) distribution with dijkstra's shortest paths algorithm

    make_p_insert(o, d, g, μ)

The shortest paths to destination ar calculated using d as the origin in the reversed graph
"""
function make_p_insert(o, d, g, μ)

    g_rev = reverse(DiGraph(g))
    distmx = LightGraphs.weights(g)
    distmx_t = distmx'

    ds_o = dijkstra_shortest_paths(g, o, distmx)
    ds_d = dijkstra_shortest_paths(g_rev, d, distmx_t)

    ps_from_o = enumerate_paths(ds_o)
    ps_to_d = enumerate_paths(ds_d)

    dists_from_o = map(x->path_length(x, g, distmx), ps_from_o)
    dists_to_d = map(x->path_length(x, g, distmx_t), ps_to_d)

    dev_dists = dists_from_o + dists_to_d
    exp_terms = exp.(-μ*dev_dists)
    denominator = sum(exp_terms)

    p_insert = exp_terms/denominator
    p_insert = StatsBase.Weights(p_insert)
end

function make_p_insert_and_ds_o(o, d, g, μ)

    g_rev = reverse(DiGraph(g))
    distmx = LightGraphs.weights(g)
    distmx_t = distmx' 

    ds_o = dijkstra_shortest_paths(g, o, distmx)
    ds_d = dijkstra_shortest_paths(g_rev, d, distmx_t)

    ps_from_o = enumerate_paths(ds_o)
    ps_to_d = enumerate_paths(ds_d)

    dists_from_o = map(x->path_length(x, g, distmx), ps_from_o)
    dists_to_d = map(x->path_length(x, g, distmx_t), ps_to_d)

    dev_dists = dists_from_o + dists_to_d
    exp_terms = exp.(-μ*dev_dists)
    denominator = sum(exp_terms)

    p_insert = exp_terms/denominator
    p_insert = StatsBase.Weights(p_insert)
    p_insert, ds_o
end
