# Test 1: Make sure it gives an answer
#A_braess = [0 1 1 0 
#            0 0 1 1
#            0 0 0 1 
#            0 0 0 0]
#A_braess = Float64.(A_braess)
#g = SimpleWeightedDiGraph(A_braess)
#d_mat = geod_dist_mat(g)

#p_ins = StatsBase.Weights(ones(nv(g)))
#p_spl = 0.3

#state1 = MicroState([1,2,3,4], 1, 2, 4)
#state2 = splice(state1, g, p_ins)

#state3 = MicroState([1,3,4], 1, 2, 3)
#state4 = MicroState([1,2,4], 1, 2, 3)

# state 1 is not slpiceable
#proposal_probability(state1, state2, p_spl, p_ins, g, d_mat)

# state 2 is spliceable
#proposal_probability(state2, state1, p_spl, p_ins, g, d_mat)

# state 3 and 4 are both spliceable
#proposal_probability(state3, state4, p_spl, p_ins, g, d_mat)

# Same state
#proposal_probability(state2, state3, p_spl, p_ins, g, d_mat)

# Should actually check this properly!!


###################################################################
###################################################################
#
#     Test network:
#
#
#  3  *------------*  4
#     |            | \
#     |            |   \
#     |            |    * 6
#     |            |   /
#     |            | /
#  2  *------------*  5
#      \          /
#        \       /
#          \   /
#            *
#            1


# Test 2: Larger Network

points = [3 0;
          2 2;
          2 4;
          5 4;
          5 2;
          8 3]
points = Float64.(points)          
          
edges_array = [(1, 2),
         (1, 5),
         (2, 3),
         (2, 5),
         (3, 4),
         (4, 5),
         (4, 6),
         (5, 6)]
         
edge_lens = [norm(points[edges_array[i][2],:] - points[edges_array[i][1],:]) for i in 1:length(edges_array)]

g = SimpleWeightedDiGraph(size(points)[1])
let g = g       
    for i in 1:length(edges_array)
        add_edge!(g, edges_array[i]..., edge_lens[i])
        add_edge!(g, reverse(edges_array[i])..., edge_lens[i])
    end
end

p_ins = StatsBase.Weights(ones(nv(g)))

d_mat = geod_dist_mat(g)

# Inverse tempr for metropolis hastings
μ = 0.5
o = 1
d = 6
p_splice = 0.5

mh = MHInstance(g, o, d, μ, p_splice)

s = mh.history[end]
s2 = splice(s, g, p_ins)

is_spliceable(s, g, mh.sp_dist_mat)
is_spliceable(s2, g, mh.sp_dist_mat)

# Testing all branches of flow chart:

s1 = MicroState([1, 5, 6], 1, 2, 3)
s2 = MicroState([1, 5, 4, 6], 1, 3, 4)
s3 = MicroState([1, 5, 4, 6], 1, 2, 4) # Not spliceable
s4 = MicroState([1, 5, 4, 6], 1, 2, 3)

s5 = MicroState([1, 2, 3, 4, 6], 1, 3, 4)
s6 = MicroState([1, 2, 3, 4, 6], 1, 2, 4)

s7 = MicroState([1, 2, 5, 4, 6], 1, 3, 5)
s8 = MicroState([1, 2, 5, 4, 6], 1, 4, 5)

# Test for prob of change in path through splice (Γ != Γ′)
@test proposal_probability(s1, s2, mh.p_splice, mh.p_insert, mh.g, mh.sp_dist_mat) == mh.p_splice*mh.p_insert[4]   

# Test for prob Γ to Γ′ through shuffle after not splicing even though poss.
@test proposal_probability(s2, s4, mh.p_splice, mh.p_insert, mh.g, mh.sp_dist_mat) ==  (1-mh.p_splice)*p_shuffle(s2.Γ, s4.a, s4.b, s4.c)

# Test for prob of Γ to Γ′ through shuffle only (because Γ′ not spliceable)
@test proposal_probability(s2, s3, mh.p_splice, mh.p_insert, mh.g, mh.sp_dist_mat) == p_shuffle(s3.Γ, s2.a, s2.b, s2.c)

# Test for prob of Γ to Γ′ through shuffle cause not spliceable (a, c unchanged)
@test proposal_probability(s7, s8, mh.p_splice, mh.p_insert, mh.g, mh.sp_dist_mat) == p_shuffle(s7.Γ, s8.a, s8.b, s8.c)

# Test for transition possible by both splicing and shuffling
@test proposal_probability(s5, s6, mh.p_splice, mh.p_insert, mh.g, mh.sp_dist_mat) == mh.p_splice*mh.p_insert[s6.Γ[s6.b]] + 
                   (1-mh.p_splice)*p_shuffle(s5.Γ, s6.a, s6.b, s6.c)

# Test for self loop detection (returns 1 to give multiplicative neutral)
@test proposal_probability(s2, s2, mh.p_splice, mh.p_insert, mh.g, mh.sp_dist_mat) == 1

##########################################
##########################################

#@test proposal_probability_verbatim(s3, s2, mh.p_splice, mh.p_insert, mh.g, mh.sp_dist_mat) == p_shuffle(s3.Γ, s2.a, s2.b, s2.c)

# Here s2 to s3 cannot happen through a splice (so something seems strange)
#@test proposal_probability_verbatim(s2, s3, mh.p_splice, mh.p_insert, mh.g, mh.sp_dist_mat) == mh.p_splice*mh.p_insert[s3.Γ[s3.b]] + (1 - mh.p_splice)*p_shuffle(s2.Γ, s3.a, s3.b, s3.c)

