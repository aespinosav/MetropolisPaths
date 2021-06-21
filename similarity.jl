"""
Similarity measure between paths in MHInstance history separated by d steps:
    
    similarity_measure(paths::Array{Array{Int,1},1}, d::Int; K)
    
Should check that it can run all the way to `K` for distance `d`
"""
function similarity_measure(path_array::Array{Array{Int,1},1}, d::Int, K::Int)
     
     M = length(path_array)
     #K = M - d
    
     suma = 0.0
     for i in 1:K
        path1 = path_array[i]
        path2 = path_array[i+d]
        
        l1 = length(path1)
        l2 = length(path2)
        
        suma += length(intersect(path1, path2))/(0.5*(l1 + l2))
     end
     suma/K
end

"""
Similarity measure between paths in MHInstance history separated by d steps:
    
    similarity_measure(mh::MHInstance, d::Int; K)
    
Should check that it can run all the way to `K` for distance `d`
"""
similarity_measure(mh::MHInstance, d::Int, K::Int) = similarity_measure(map(x->x.Γ, mh.history), d, K)


"""
Calculates similarity measures for an array of `d`s (d_range)

    similarity_curve(path_array::Array{Array{Int,1},1}, d_range, K)
"""
function similarity_curve(path_array::Array{Array{Int,1},1}, d_range, K)
    s_array = zeros(length(d_range))
    for (i,d) in enumerate(d_range)
        s_array[i] = similarity_measure(path_array, d, K)
        #println("finished calc number $i")
    end
    s_array
end



###
### Old definitions to be removed at some point I guess, although the one withouth explicitly defining K 
### might be useful, or the more natural definition for then defining the other methods
###

#########################################################################################################
function similarity_measure_dep(mh::MHInstance, d::Int, K::Int)
     
     M = length(mh.history)
     
     suma = 0.0
     for i in 1:K
        path1 = mh.history[i].Γ
        path2 = mh.history[i+d].Γ
        
        l1 = length(path1)
        l2 = length(path2)
        
        suma += length(intersect(path1, path2))/(0.5*(l1 + l2))
     end
     suma/K
end

"""
Runs over all the states in the history of the MHInstance mh
"""
similarity_measure(mh::MHInstance, d::Int) =  similarity_measure(mh, d, length(mh.history)-d)


function similarity_measure_old(mh::MHInstance, d::Int)
     
     M = length(mh.history)
     K = M - d
     
     suma = 0.0
     for i in 1:K
        path1 = mh.history[i].Γ
        path2 = mh.history[i+d].Γ
        
        l1 = length(path1)
        l2 = length(path2)
        
        suma += length(intersect(path1, path2))/(0.5*(l1 + l2))
     end
     suma/K
end

