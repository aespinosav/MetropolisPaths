module MHPaths

using LightGraphs, SimpleWeightedGraphs, StatsBase

export 
    # From microstate.jl    
    MicroState,
    # From splice.jl
    splice

include("microstate.jl")
include("splice.jl")

end
