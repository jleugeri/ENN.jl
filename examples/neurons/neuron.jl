using TimedAutomata, Base.Iterators
using GLMakie, GraphMakie
##
abstract type SomaOrDendrite end

struct DendriteSegment <: SomaOrDendrite
    children::Vector{DendriteSegment}
    parent::Ref{SomaOrDendrite}
    input::Symbol
    dendritic_threshold::Int
    function DendriteSegment(children::Vector{DendriteSegment}, input::Symbol, dendritic_threshold::Int=isempty(children) ? 0 : 1)
        if dendritic_threshold > length(children)
            @warn "Threshold '$(dendritic_threshold)' is larger than number of children '$(length(children))'."
        end
        
        this = new(children,Ref{SomaOrDendrite}(), input, dendritic_threshold)
        for child in children
            child.parent[]=this
        end
        this
    end
end

function DendriteSegment(children, input::Symbol, dendritic_threshold::Int=isempty(children) ? 0 : 1)
    new_children = Vector{DendriteSegment}()
    for child in children
        push!(new_children, DendriteSegment(child...))
    end
    invoke(DendriteSegment, Tuple{Vector{DendriteSegment},Symbol,Int}, new_children, input, dendritic_threshold)
end

struct Neuron <: SomaOrDendrite
    all_segments::Vector{DendriteSegment}
    dendrite::DendriteSegment
    function Neuron(args...; kwargs...)
        dendrite = DendriteSegment(args...; kwargs...)

        # store all segments in BFS-order
        head=0
        all_segments = [dendrite]
        while head < length(all_segments)
            head+=1
            d=all_segments[head]
            append!(all_segments, d.children)
        end

        this = new(all_segments,dendrite)
        dendrite.parent[]=this
        this
    end
end

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

##
neuron_sequential = Neuron((((((),:A),),:B),),:C);
neuron_sequential_ta = TA(neuron_sequential)
neuron_sequential_p = prune(neuron_sequential_ta)

#zone = TAConstraint(neuron_sequential_p.clocks, zeros(TABound{Int},length(neuron_sequential_p.clocks)+1,length(neuron_sequential_p.clocks)+1))
zone = @constraint clk_1-clk_2==0 && clk_2-clk_3==0 && clk_3-clk_1==0 Int
TimedAutomata.set_clock_list!(zone, neuron_sequential_ta.clocks)
neuron_sequential_tts = TTS(neuron_sequential_p; initial_zone=zone)
@assert all(neuron_sequential_tts.states[end][1].components[1])
sequential_language_full=language(neuron_sequential_tts; final_states=neuron_sequential_tts.states[end:end], ignore_messages=(_)->false)
sequential_language_no_silent=language(neuron_sequential_tts; final_states=neuron_sequential_tts.states[end:end])
sequential_language_input=language(neuron_sequential_tts; final_states=neuron_sequential_tts.states[end:end], forbid_messages=(m)->m.symbol∈(:↻,), ignore_messages=(m)->m.direction!=MSG_IN)
sequential_language_no_reset_no_spike=language(neuron_sequential_tts; final_states=neuron_sequential_tts.states[end:end], forbid_messages=(m)->m.symbol∈(:spike,:↻))

##
neuron_parallel = Neuron((((),:A),((),:B)),:C);
neuron_parallel_ta = TA(neuron_parallel)
neuron_parallel_p = prune(neuron_parallel_ta)
neuron_parallel_tts = TTS(neuron_parallel_p; initial_zone=zone)
id = findfirst(((l,_),)->all(l.components[1]), neuron_parallel_tts.states)
parallel_language_full=language(neuron_parallel_tts; final_states=neuron_parallel_tts.states[id:id], ignore_messages=(_)->false)
parallel_language_no_silent=language(neuron_parallel_tts; final_states=neuron_parallel_tts.states[id:id])
parallel_language_input=language(neuron_parallel_tts; final_states=neuron_parallel_tts.states[id:id], forbid_messages=(m)->m.symbol∈(:↻,), ignore_messages=(m)->m.direction!=MSG_IN)
parallel_language_no_reset_no_spike=language(neuron_parallel_tts; final_states=neuron_parallel_tts.states[id:id], forbid_messages=(m)->m.symbol∈(:spike,:↻))

## Plotting

f,ax,_ = graphplot(neuron_sequential_ta);
hidedecorations!(ax)
display(f)
print("Press key to continue: ")
read(stdin, 1)
save(joinpath("examples","neurons", "figures","sequential_unpruned.png"), f)

f,ax,_ = graphplot(neuron_sequential_p);
hidedecorations!(ax)
display(f)
print("Press key to continue: ")
read(stdin, 1)
save(joinpath("examples","neurons", "figures","sequential_pruned.png"), f)

f,ax,_ = graphplot(neuron_parallel_p);
hidedecorations!(ax)
display(f)
print("Press key to continue: ")
read(stdin, 1)
save(joinpath("examples","neurons", "figures","parallel_pruned.png"), f)
