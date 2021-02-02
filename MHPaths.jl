module MHPaths

using LightGraphs, SimpleWeightedGraphs, StatsBase

export 
    # From microstate.jl    
    MicroState,
    # From splice.jl
    splice,
    # From shuffle.jl
    shuffle!

include("microstate.jl")
include("splice.jl")
include("shuffle.jl")

end
