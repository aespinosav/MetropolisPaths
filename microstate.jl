"""
Type for representing microstates (of paths).
The graph to which the paths refer is external to the type.
"""
mutable struct MicroState 
    Γ::Array{Int64,1} # Path
    a::Int64 # First anchor
    b::Int64 # Splice point 
    c::Int64 # Second anchor
    MicroState(Γ, a, b, c) = a<b && b<c ? new(Γ, a, b, c) : error("Unordered anchors: a, b, c")
end
