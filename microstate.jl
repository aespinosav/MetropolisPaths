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
