export close!, close, up!, up, reset!, reset, zone_graph


function close!(c::TAConstraint)
    n = size(c.D,1)
    for k ∈ 1:n
        for i ∈ 1:n
            for j ∈ 1:n
                c.D[i,j] = min(c.D[i,j], c.D[i,k]+c.D[k,j])
            end
        end
    end
    c
end
close(c::TAConstraint) = close!(copy(c))

function up!(c::TAConstraint)
    c.D[2:end,1] .= typemax(eltype(c.D))
    close!(c)
    c
end
up(c::TAConstraint) = up!(copy(c))

function reset!(c::TAConstraint, clocks::Set{Symbol})
    for clock ∈ clocks
        idx = findfirst(==(clock), c.clocks)+1
        c.D[idx,:] .= c.D[1,:]
        c.D[:,idx] .= c.D[:,1]
    end
    close!(c)
    c
end

reset(c::TAConstraint, clocks::Set{Symbol}) = reset!(copy(c), clocks)

function zone_graph(ta::TA{State, T}) where {State, T}
    states = Tuple{State,TAConstraint{T}}[(ta.initial_state, TAConstraint(ta.clocks, zeros(TABound{T}, length(ta.clocks)+1,length(ta.clocks)+1)))]
    transitions = Pair{Tuple{State,TAConstraint{T}},Tuple{TAMessage,Tuple{State,TAConstraint{T}}}}[]
    max_runs = 1e8

    function covers((location,constraint))
        ((l,c),) -> location == l && constraint ⊆ c
    end

    last = 0
    while last < length(states) && max_runs > 0
        max_runs -=1
        last += 1
        (l,D) = states[last]

        # Delay transition:
        # (l,D) ⇝ (l,D′), D′ = D↑ ∩ I(l)
        l′ = l
        D′ = up(D) ∩ ta.invariants[l]
        
        # if this state is consistent, keep it
        if !isempty(D′)

            idx = findfirst(covers((l′,D′)), states)
            if isnothing(idx)
                # if we just found a new state + clock region, save for later!
                push!(states, (l′,D′))
            else
                (l′,D′) = states[idx]
            end
            push!(transitions, (l,D) => (TAMessage() , (l′,D′)))
        end
        
        # Action transitions:
        # (l,D) ⇝ (l′, D′), D′ = r(D∧g) ∧ I(l′)  if  l →g,a,r→ l′
        for arc ∈ ta.arcs_by_state[l]
            # follow the arc
            l′ = arc.target
            
            # The arc constrains the zone
            D′′ = D ∩ arc.guard
            
            # if the arc is compatible with the current time zone, keep going
            if !isempty(D′′)
                # the transition resets some clocks and must satisfy the new state's guards
                D′ = TimedAutomata.reset(D′′, arc.resets) ∩ ta.invariants[l′]

                idx = findfirst(covers((l′,D′)), states)
                if isnothing(idx)
                    # if we just found a new state + clock region, save for later!
                    push!(states, (l′,D′))
                else
                    (l′,D′) = states[idx]
                end
                push!(transitions, (l,D) => (TAMessage() , (l′,D′)))
            end
        end
    end


    return (states=states, transitions=transitions)
end
