"""
Calculates proposal probability `q(i, i′)` according to flowchart from
Flötteröd & Bierlaire. It is non-trivial because it has to account for
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
                if is_spliceable(old_state, g, M)
                    q = p_splice*p_insert(ν) + 
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
