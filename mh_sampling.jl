"""
Initial function to test if sample works
"""
function weight_func(state)
    length(state.Γ) 
end


"""
Samples path using the Metropolis Hastings procedure.

    mh_sample(mh, steps, b)
    
    mh:             MHInstance object
    steps:          Number of samples to draw
    weight_func:    Weight function associated to a state
    
At the moment it is coded to work with the cost function 
given in Flötteröd & Bierlaire that is passed as `weight_func`.

That is: 
    
    weight_func(state) =  exp(-μ*δ(Γ))/(|Γ|(|Γ|-1)(|Γ|-2)/6)

with δ(Γ) the path cost

... Have to check what modifications are needed if used for a 
different function.

I think weight_func should probably be inside MHInstance 
"""
function mh_sample(mh, weight_func)

    current_state = mh.history[end]
    
    draw = rand()
    if is_spliceable(current_state, mh.g, mh.sp_dist_mat) && draw<mh.p_splice
        candidate_state = splice(mh.g, mh.p_insert, current_state)
    else
        candidate_state = shuffle(current_state)
    end
    
    # Probability of drawing transition
    q_ij = proposal_probability(current_state,
                                candidate_state,
                                mh.p_splice,
                                mh.p_insert,
                                mh.g,
                                mh.sp_dist_mat)
    
    q_ji = proposal_probability(candidate_state,
                                current_state,
                                mh.p_splice,
                                mh.p_insert,
                                mh.g,
                                mh.sp_dist_mat)
    
    # Acceptance probability
    len_Γ = length(current_state.Γ)
    len_Γ′= length(candidate_state.Γ)
    
    b_i = weight_func(current_state)/(len_Γ*(len_Γ-1)*(len_Γ-2)/6)
    b_j = weight_func(candidate_state)/(len_Γ′*(len_Γ′-1)*(len_Γ′-2)/6)
    
    # Detailed balance?
    quotient =  b_j*q_ji/b_i*q_ij
    
    p_accept = minimum([quotient, 1])   
    
    # Accept or reject
    draw = rand()
    if draw<p_accept
        # Accept!
        return candidate_state
    else
        # Reject :(
        return deepcopy(current_state)
        # Is this bad for performance? I am just being careful
    end
end

"""
Draw samples with Metropolis Hastings algorithm and add to MHInstance history

    mh_evolve!(mh, num_steps, weight_func)
    
    mh:             MHInstance
    num_steps:      Number of samples
    weight_func:    Function for prob weights
"""
function mh_evolve!(mh, num_steps, weight_func)
    k = 0
    while k<num_steps
        k += 1
        new_state = mh_sample(mh, weight_func)
        push!(mh.history, new_state) 
    end
end
