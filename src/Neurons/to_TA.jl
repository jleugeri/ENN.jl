import ..TimedAutomata

function TimedAutomata.TA(n::Neuron, name=:spike, τ_plateau::T=100, reset_only_necessary=false) where T
    N = length(n.all_segments)
    states = BitVector.(product(([false,true] for i in 1:N)...))[:]
    initial_state=states[1]
    symbols=Set(seg.input for seg in n.all_segments)
    clocks=[Symbol("clk_$(i)") for i in 1:N] # 🕑

    lookup = Dict(zip(n.all_segments, eachindex(n.all_segments)))

    function get_subtree(id)
        ids = Int[id]
        for child in n.all_segments[id].children
            append!(ids, get_subtree(lookup[child]))
        end
        return ids
    end

    # Regular transitions
    arcs = TAArc{BitVector, T}[]
    invariants = Dict{BitVector,TAConstraint{T}}()
    for (i,state) in enumerate(states)
        state_invariants = TADiagonalConstraint{T}[]
        # each segment could turn on, if > threshold of the children are on
        for (j,active) in enumerate(state)
            segment = n.all_segments[j]
            if active # allow segment to turn OFF when time runs out...
                # ... if its parent is not currently on (backprop prevents turning off)
                if !isa(segment.parent.x,Neuron) && !state[lookup[segment.parent[]]]
                    ids=get_subtree(j)
                    new_state = copy(state)
                    new_state[ids] .= false

                    # timeout happens exactly on the plateau-potential boundary
                    g = TAConstraint([
                            TADiagonalConstraint(clocks[j],nothing,TABound(false,τ_plateau)),
                            TADiagonalConstraint(nothing,clocks[j],TABound(false,-τ_plateau))
                        ], 
                        clocks
                    )
                    resets = reset_only_necessary ? Set{Symbol}() : Set(clocks[ids]) ∪ Set(clocks[.!state])
                    msg = TAMessage(MSG_SILENT, :↻)
                    push!(arcs, TAArc(state, new_state, g, msg, resets))

                    # set invariant accordingly
                    push!(state_invariants, TADiagonalConstraint(clocks[j],nothing,TABound(false,τ_plateau)))
                end
            else # allow segment to turn ON ...
                # ... if there is sufficient dendritic input from child-segments,
                active_children = 0
                for child in segment.children
                    k=lookup[child]
                    active_children += state[k]
                end
                
                # make an edge to trigger the next segment from synaptic input
                if active_children ≥ segment.dendritic_threshold
                    ids=get_subtree(j)
                    new_state = copy(state)
                    new_state[ids] .= true
                    # spikes can't arrive at exactly the same time
                    g = TAConstraint([TADiagonalConstraint(nothing,clocks[j],TABound(true,zero(T)))], clocks)
                    # reset all clocks of subtree
                    resets = reset_only_necessary ? Set([j]) : Set(clocks[ids]) ∪ Set(clocks[.!state])
                    push!(arcs, TAArc(state, new_state, g, TAMessage(MSG_IN, segment.input), resets))
                end
            end
        end
        invariants[state] = TAConstraint(state_invariants, clocks)
    end
        

    ## Handle the special cases for the spiking state
    spiking_state = states[findfirst(all,states)]
    
    # timeout happens immediately after the soma turned on
    g = TAConstraint([
            TADiagonalConstraint(clocks[lookup[n.dendrite]],nothing,TABound(false,0)),
            TADiagonalConstraint(nothing,clocks[lookup[n.dendrite]],TABound(false,0))], 
        clocks)
    
    resets = reset_only_necessary ? Set{Symbol}() : Set(clocks)
    # Add special transition for somatic spike!
    push!(arcs, TAArc(spiking_state, initial_state, g, TAMessage(MSG_OUT, name), resets))

    # Add special invariant for somatic AP
    invariants[spiking_state] = deepcopy(g)

    TA(states, initial_state, clocks, symbols, arcs, invariants)
end
