"""
Similarity measure between paths in MHInstance history separated by d steps:
    
    similarity_measure(mh::MHInstance, d::Int)
"""
function similarity_measure(mh::MHInstance, d::Int)
     
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
