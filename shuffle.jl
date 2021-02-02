"""
Shuffles anchor points (a, c) and pivot of a microstate for generating 
proposed states for transition in MH algorithm as in the paper 
[Flötteröd & Bierlaire 2013 (TR: Part B) section 2.5]
"""
function shuffle!(state)
    N = length(state.Γ)
    weights = Weights(ones(N))
    new_anchors = sample(1:N, weights, 3, replace=false, ordered=true)
    state.a, state.b, state.c = new_anchors
    state
end
