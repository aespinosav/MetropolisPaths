#Adjacency matrix for a braess network
using LightGraphs, SimpleWeightedGraphs

A_braess = [0 1 1 0 
            0 0 1 1
            0 0 0 1 
            0 0 0 0]

g = SimpleWeightedDiGraph(A_braess)
