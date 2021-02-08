"""
Calculates proposal probability `q(i, i′)` according to flowchart from
Flötteröd & Bierlaire. 

    proposal_probability(old_state, new_state, p_splice, p_insert, g, M)

    old_state:  Original state (state `i`)
    new_state:  Proposed state (sampled) (state `j`)
    p_splice:   Probability of splicing (constant)
    p_insert:   Probability distribution over nodes (Array)
    g:          Graph
    M:          Shortest-path (geodesic) distance matrix

It is non-trivial because it has to account for
all possible ways of transitioning from i→i′.
"""
function proposal_probability(old_state, new_state, p_splice, p_insert, g, M)
    Γ, Γ′ = old_state.Γ, new_state.Γ
    a, b, c = old_state.a, old_state.b, old_state.c
    a′, b′, c′ = new_state.a, new_state.b, new_state.c
    
    ν = Γ′[b′]
    
    if Γ==Γ′
        if a==a′ && c==c′
            if b==b′
                q = 1 # Cancels out when calc. acceptance prob (self-loop)
            else
                if is_spliceable(new_state, g, M)
                    q = p_splice*p_insert[ν] + 
                        (1-p_splice)*p_shuffle(Γ, a′, b′, c′)
                else
                    q = p_shuffle(Γ, a′, b′, c′)
                end
            end
        else
            if is_spliceable(old_state, g, M)
                #This is different from flowchart but I suspect there is a typo
                #in the F&B paper
                q = (1-p_splice)*p_shuffle(Γ, a′, b′, c′)
            else
                q = p_shuffle(Γ, a′, b′, c′)
            end
        end
    else
        q = p_splice*p_insert[ν]
    end
end

###
### Note: Altenrate function (most lileley wrong)
### 

"""
This is the way of calculating the `proposal_probability` verbatim from the 
paper, however I suspec there is a typo in the flowchart, which should check 
whether the new state is spliceable to see if it was able to come from a
splice operation (the `if-then-else` block in the `else` of the a==a′ block)
"""
function proposal_probability_verbatim(old_state, new_state, p_splice, p_insert, g, M)
    Γ, Γ′ = old_state.Γ, new_state.Γ
    a, b, c = old_state.a, old_state.b, old_state.c
    a′, b′, c′ = new_state.a, new_state.b, new_state.c
    
    ν = Γ′[b′]
    
    if Γ==Γ′
        if a==a′ && c==c′
            if b==b′
                q = 1 # Cancels out when calc. acceptance prob (self-loop)
            else
                if is_spliceable(old_state, g, M)
                    q = p_splice*p_insert[ν] + 
                        (1-p_splice)*p_shuffle(Γ, a′, b′, c′)
                else
                    q = p_shuffle(Γ, a′, b′, c′)
                end
            end
        else
            if is_spliceable(old_state, g, M)
                q = (1-p_splice)*p_shuffle(Γ, a′, b′, c′)
            else
                q = p_shuffle(Γ, a′, b′, c′)
            end
        end
    else
        q = p_splice*p_insert[ν]
    end
end
