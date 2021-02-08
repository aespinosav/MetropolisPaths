"""
Shuffles anchor points (a, c) and pivot of a microstate for generating 
proposed states for transition in MH algorithm as in the paper 
[Flötteröd & Bierlaire 2013 (TR: Part B) section 2.5].

    shuffle(state)

Returns new state
"""
function shuffle(state)
    N = length(state.Γ)
    weights = Weights(ones(N))
    new_anchors = sample(1:N, weights, 3, replace=false, ordered=true)
    new_state = MicroState(state.Γ, new_anchors...)
end

"""
Shuffles anchor points (a, c) and pivot of a microstate for generating 
proposed states for transition in MH algorithm as in the paper 
[Flötteröd & Bierlaire 2013 (TR: Part B) section 2.5].

    shuffle!(state)
        
Changes input `state`.
"""
function shuffle!(state)
    N = length(state.Γ)
    weights = Weights(ones(N))
    new_anchors = sample(1:N, weights, 3, replace=false, ordered=true)
    state.a, state.b, state.c = new_anchors
    state
end

"""
Probability of having chosen anchor points a, b and c. It is a function
of the length of the path.

    p_shuffle(path, a, b, c)
"""
function p_shuffle(path, a, b, c)
    N = length(path)
    if a<b && b<c
        return 1.0/(N*(N-1)*(N-2)/6)
    else
        return 0
    end
end
