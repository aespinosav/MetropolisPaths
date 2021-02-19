"""
Initial function to test if sample works
"""
function weight_func(state::MicroState, mh::MHInstance)
    g = mh.g
    μ = mh.μ
    exp(-μ*path_length(state.Γ, g))
end

function weight_func(path::Array{Int,1}, mh::MHInstance)
    g = mh.g
    μ = mh.μ
    exp(-μ*path_length(path, g))
end

function weight_func1(state, g)
    path_length(state.Γ, g) 
end


"""
Samples path using the Metropolis Hastings procedure.

    mh_sample(mh, weight_func)
    
    mh:             MHInstance object
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
        candidate_state = splice(current_state, mh.g, mh.p_insert)
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
    
    b_i = weight_func(current_state, mh)/(len_Γ*(len_Γ-1)*(len_Γ-2)/6)
    b_j = weight_func(candidate_state, mh)/(len_Γ′*(len_Γ′-1)*(len_Γ′-2)/6)
    
    # Detailed balance?
    quotient =  (b_j*q_ji) / (b_i*q_ij)
    
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

###
### This next function here is just vestigial, should probably get rid 
### of it at some point.
###

function mh_sample_simple(mh, weight_func)

    current_state = mh.history[end]
    println("Current state: ", current_state)

    draw = rand()
    println("draw: ", draw)
    if is_spliceable(current_state, mh) && draw<mh.p_splice
        println("We're splicing!")
        candidate_state = splice(current_state, mh)
    else
        println("We're shuffling!")
        candidate_state = shuffle(current_state)
    end
    println("candidate state: ", candidate_state)

    # Probability of drawing transition
    q_ij = proposal_probability(current_state,
                                candidate_state,
                                mh)
    println("q(i,j): ", q_ij)

    q_ji = proposal_probability(candidate_state,
                                current_state,
                                mh)
    println("q(j,i): ", q_ji)
    # Acceptance probability

    # number of nodes in path
    len_Γ = length(current_state.Γ)
    norm_1 = len_Γ*(len_Γ-1)*(len_Γ-2)/6

    len_Γ′= length(candidate_state.Γ)
    norm_2 = len_Γ′*(len_Γ′-1)*(len_Γ′-2)/6

    b_i = weight_func(current_state, mh) / norm_1
    b_j = weight_func(candidate_state, mh) / norm_2

    println("b(i): ", b_i)
    println("b(j): ", b_j)

    # TODO: change this to the logarithm comparison...
    # Detailed balance?
    quotient =  (b_j*q_ji)/(b_i*q_ij)
    println("t prob: ", quotient)

    p_accept = minimum([quotient, 1])   

    # Accept or reject
    draw = rand()
    println("draw: ", draw)
    if draw<p_accept
        # Accept!
        println("Transition accepted")
        return candidate_state
    else
        # Reject :(
        println("Transition rejected")
        return deepcopy(current_state)
        # Is this bad for performance? I am just being careful
    end
end


"""
Samples path using the Metropolis Hastings procedure.

    mh_sample(mh, weight_func)
    
    mh:             MHInstance object
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
function mh_sample_log(mh, weight_func)

    current_state = mh.history[end]
    #println("Current state: ", current_state)

    draw = rand()
    #println("draw: ", draw)
    if is_spliceable(current_state, mh) && draw<mh.p_splice
        #println("We're splicing!")
        candidate_state = splice(current_state, mh)
    else
        #println("We're shuffling!")
        candidate_state = shuffle(current_state)
    end
    #println("candidate state: ", candidate_state)

    # Probability of drawing transition
    q_ij = proposal_probability(current_state,
                                candidate_state,
                                mh)
    #println("q(i,j): ", q_ij)    

    q_ji = proposal_probability(candidate_state,
                                current_state,
                                mh)
    #println("q(j,i): ", q_ji)
    # Acceptance probability

    # number of nodes in path
    len_Γ = length(current_state.Γ)
    norm_1 = len_Γ*(len_Γ-1)*(len_Γ-2)/6

    len_Γ′= length(candidate_state.Γ)
    norm_2 = len_Γ′*(len_Γ′-1)*(len_Γ′-2)/6

    b_i = weight_func(current_state, mh) / norm_1
    b_j = weight_func(candidate_state, mh) / norm_2

    #println("b(i): ", b_i)
    #println("b(j): ", b_j)

    # Detailed balance?

    quotient =  (b_j*q_ji)/(b_i*q_ij)

    # Since probs will be small, compare log
    log_quotient = log(b_j) - log(b_i) + log(q_ji) - log(q_ij)

    # Accept or reject
    draw = rand()
    log_draw = log(draw)
    
    if log_draw < log_quotient
        # Accept!
        #println("Transition accepted")
        return candidate_state
    else
        # Reject :(
        #println("Transition rejected")
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
    
By default, uses the `mh_sample_log` function.
"""
function mh_evolve!(mh, num_steps, weight_func)
    k = 0
    while k<num_steps
        k += 1
        #new_state = mh_sample(mh, weight_func)
        new_state = mh_sample_log(mh, weight_func)
        push!(mh.history, new_state) 
    end
end
