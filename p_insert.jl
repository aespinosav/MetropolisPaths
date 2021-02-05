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

"""
Returns function where denominator does not need to be passed explicitly
"""
function make_p_insert_with_denom(o, d, g, μ, d_mat)
    den = p_insert_denom(o, d, μ, g, d_mat)
    v -> p_insert_explicit(v, o, d, μ, d_mat, den)
end
