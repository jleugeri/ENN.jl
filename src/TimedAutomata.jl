module TimedAutomata
using LightGraphs, GraphPlot, Base.Iterators
export TimeSlice,TimePoint,TimeInterval,TimeLeftHalfLine,TimeRightHalfLine, ClockRegion, collect_constants, @Clk_str,
    MA,TMA,MATransition,TMATransition, to_graph, untime, language, ⊆, ==, PartialOrder, TotalOrder, OrderRelation, gplot
#, TMARun, RunState

##########################
#    Clock Conditions    #
##########################
include("clock_constraints.jl")


##########################
#     Timed Automata     #
##########################

struct TMATransition{TMAState}
    s::TMAState
    λ::Vector{Int}
    δ::ClockRegion
end

struct TMA{TMAState}
    Σ::Vector{TMAState}
    Σ′::Vector{TMAState}
    num_clocks::Int
    E::Dict{Tuple{TMAState,Symbol}, Vector{TMATransition{TMAState}}}
    I::Dict{TMAState, <:ClockRegion}
    s0::TMAState
end

#= 
mutable struct RunState{TMAState}
    state::TMAState
    clocks::Vector{Float64}
end

Base.:(==)(s1::RunState{TMAState}, s2::RunState{TMAState}) where {TMAState} = all(s1.state .== s2.state) && all(s1.clocks .≈ s2.clocks)
 =#

##########################
#   Nontimed Automata    #
##########################
struct MATransition{MAState}
    s::MAState
end

struct MA{MAState}
    Σ::Vector{MAState}
    Σ′::Vector{MAState}
    E::Dict{Tuple{MAState,Symbol}, Vector{MATransition{MAState}}}
    s0::MAState
end



function collect_constants(tma::TMA) 
    T = typeof(tma.I).parameters[2].parameters[2]
    res=[Set{T}(0) for i ∈ 1:tma.num_clocks]
    
    for i ∈ values(tma.I)
        collect_constants!(i, res)
    end

    for es ∈ values(tma.E)
        for e ∈ es
            collect_constants!(e.δ, res)
        end
    end
    
    return res
end

##########################
#  Untiming Operations   #
##########################
struct UMAState{TMAState, N, T}
    tma_state::TMAState
    clock_region::ClockRegion{N,T}
end
Base.:(==)(s1::UMAState{TMAState,N,T},s2::UMAState{TMAState,N,T}) where{TMAState,N,T} = (s1.tma_state == s2.tma_state) && (s1.clock_region == s2.clock_region)

function untime(m::TMA{TMAState}) where TMAState
    num_clocks = m.num_clocks
    
    # determine limits and scaling of the clocks

    clock_constants = collect_constants(m)
    clock_scaler = lcm(denominator.(reduce(union,clock_constants))...)//gcd(numerator.(reduce(union,clock_constants))...)
    if denominator(clock_scaler)==0
        clock_scaler = 0//1
    end
    c_max = Int.(maximum.(clock_constants) .* clock_scaler)

    """Utility function for updating clock regions in the setting with integer constraints"""
    function update_clock_region(r::ClockRegion, increase=ones(Bool, num_clocks), reset=zeros(Bool,num_clocks))
        new_slices = deepcopy(r.slices)
        new_order = deepcopy(r.order.order)
        for (i,s) ∈ enumerate(r.slices)
            new_slices[i] = if reset[i]
                TimePoint(0)
            elseif increase[i]
                if isa(s,TimePoint)
                    if s.val < c_max[i]
                        TimeInterval{Int,true,true}(s.val,s.val+1)
                    else
                        TimeRightHalfLine{Int,true}(s.val)
                    end
                elseif isa(s,TimeInterval)
                    TimePoint(s.val2)
                elseif isa(s,TimeRightHalfLine)
                    error("Can't go larger than $(s) for clock $(i)")
                else
                    error("$(s) should be TimePoint, TimeInterval or TimeRightHalfLine")
                end
            else
                new_slices[i]
            end
        end
        
        # if the clock region is a constant or out of bounds, treat it as exact and ignore order
        exact = isa.(new_slices, TimePoint)
        is_out_of_bound = isa.(new_slices, TimeRightHalfLine)
        unordered = exact .| is_out_of_bound

        # update the order of all of the clocks' fractional parts that are not exactly 0
        new_ranks = sort!(unique(new_order[.~(unordered)]))
        new_order[.~(unordered)] .= indexin(new_order[.~(unordered)], new_ranks)
        new_order[unordered] .= 0

        return ClockRegion{length(new_slices),Int}(new_slices, TotalOrder(new_order))
    end

    # construct region automaton
    states = UMAState{TMAState,num_clocks,Int}[]
    accepting_states = UMAState{TMAState,num_clocks,Int}[]
    edges = Dict{Tuple{UMAState{TMAState,num_clocks,Int},Symbol}, Vector{MATransition{UMAState{TMAState,num_clocks,Int}}}}()

    s0 = UMAState(m.s0,ClockRegion{num_clocks,Int}(fill(TimePoint(0), num_clocks), TotalOrder(zeros(Int,num_clocks))))

    queue = [s0]
    while ~isempty(queue)
        q = popfirst!(queue)
        
        # If we took care of the new state before, move on
        if q ∈ states
            continue
        end

        # determine properties of individual clocks
        is_exact = isa.(q.clock_region.slices,TimePoint)
        is_largest_rank = ismaximum(q.clock_region.order)
        is_out_of_bound = isa.(q.clock_region.slices,TimeRightHalfLine)

        # go through all symbol transitions
        for ((s′,a),ee) ∈ pairs(m.E)
            if s′ == q.tma_state    
                for e ∈ ee
                    # if there is an symbol transition, reset the corresponding clocks (both integral and fractional part)
                    # check if edge constraint holds within the region
                    if q.clock_region ⊆ e.δ
                        tmp = zeros(Bool, num_clocks)
                        tmp[e.λ] .= true
                        r′ = update_clock_region(q.clock_region, zeros(Bool,num_clocks), tmp)
                        q′ = UMAState(e.s, r′)
                        el = get(edges, (q,a), MATransition{UMAState{TMAState,num_clocks,Int}}[])
                        push!(el, MATransition(q′))
                        edges[(q,a)] = el
                        # If we didn't take care of the new state before, put it in the queue
                        if (q′ ∉ states) && (q′ ∉ queue)
                            push!(queue, q′)
                        end
                    end
                end
            end
        end
        
        # if all clocks are already out of bounds there's no time transitions
        # otherwise there is one if and only if the timing constraints are satisfied       
        if ~all(is_out_of_bound)
            # if any of the clocks are exactly specified, those clocks need to change
            # otherwise the clocks with (equal) largest rank will change first
            increase_clocks = if any(is_exact)
                is_exact
            else
                is_largest_rank
            end

            rsucc = update_clock_region(q.clock_region, increase_clocks)
            # time successor is only triggered if invariants allow it
            if q.tma_state ∉ keys(m.I) || rsucc ⊆ m.I[q.tma_state]
                q′=UMAState(q.tma_state,rsucc)
                edges[(q,:τ)] = [MATransition(q′)]
                push!(queue, q′)
            end
        end
        
        # add state to the list
        push!(states, q)
        # if the state is an accepting state, add it
        if q.tma_state ∈ m.Σ′
           push!(accepting_states, q)
        end
    end
 
    return MA(states, accepting_states, edges, s0)
end

function language(ma::MA, max_depth, prefix=Symbol[], s0=ma.s0, merge_transitions=[:τ])
    res = if max_depth == 0
        s0 ∈ ma.Σ′ ? [prefix] : []
    else
        if length(prefix)>0 && prefix[end] ∈ merge_transitions
            [[language(ma, max_depth-1, a==prefix[end] ? prefix : [prefix; a], e.s, merge_transitions) for ((s1,a),ee) ∈ ma.E for e ∈ ee if s1==s0]...;]
        else
            [[language(ma, max_depth-1, [prefix; a], e.s, merge_transitions) for ((s1,a),ee) ∈ ma.E for e ∈ ee if s1==s0]...;]
        end
    end

    return unique(res)
end

#=
struct TMARun{TMAState, DropRepetitions}
    tma::TMA{TMAState}
    inputs::Vector{Event}
    X::RunState{TMAState}
end


TMARun(tma::TMA{TMAState},inputs; drop_repetitions=false) where {TMAState} = TMARun{TMAState, drop_repetitions}(tma,inputs,RunState(deepcopy(tma.s0), zeros(Float64,tma.num_clocks)))

Base.IteratorSize(::Type{TMARun{S,true}}) where {S} = Base.SizeUnknown()
Base.IteratorSize(::Type{TMARun{S,false}}) where {S} = Base.HasLength()
Base.length(r::TMARun{S,false}) where {S} = length(r.inputs)+1

function Base.iterate(r::TMARun{S,DropRepetitions}, t0i=(0.0,0)) where {S,DropRepetitions}
    t0,i = t0i

    if i==0
        r.X.state = deepcopy(r.tma.s0)
        r.X.clocks = zeros(Float64,r.tma.num_clocks)

        return (t0,deepcopy(r.X)),(0.0, 1) 
    end
    
    t1 = t0

    keep_going = true
    while keep_going
        ret=iterate(r.inputs,i)
        
        if ret == nothing
            return nothing
        end

        (t1,a),i = ret

        dt=t1-t0
        t0=t1

        r.X.clocks .+= dt
        for e ∈ get(r.tma.E, (r.X.state, a), S[])
            if e.δ(r.X.clocks)
                r.X.clocks[e.λ] .= 0
                r.X.state = deepcopy(e.s)

                keep_going = false
                break
            end
        end

        if ~DropRepetitions
            keep_going = false
        end
    end
    return (t1,deepcopy(r.X)),(t1,i)
end

=#

function to_graph(ma::T) where T<:Union{MA,TMA}
    edges = Dict{Pair{Int,Int},Tuple{Symbol,Any}}()
    states_u = collect(ma.Σ)
    l = LightGraphs.SimpleDiGraph(length(states_u))

    for ((s,a),ee) ∈ pairs(ma.E)
        for e ∈ ee
            s_idx = findfirst(x->x==s,states_u)
            t_idx = findfirst(x->x==e.s,states_u)
            add_edge!(l,s_idx , t_idx)
            edges[s_idx=>t_idx] = (a,e)
        end
    end

    return (graph=l, edges=edges, vertices=states_u)
end

function GraphPlot.gplot(ma::MA{UMAState{TMAState,N,T}}) where {TMAState,N,T}
    g = to_graph(ma)

    tma_states = [v.tma_state for v ∈ g.vertices]
    clock_regions = [v.clock_region for v ∈ g.vertices]
    tma_states_u = unique(tma_states)
    clock_regions_u = unique(clock_regions)
    s_id = Vector{Int}(indexin(tma_states,tma_states_u))
    r_id = Vector{Int}(indexin(clock_regions,clock_regions_u))
    vertex_labels = ["S$(s)-R$(r)" for (s,r) ∈ zip(s_id, r_id)]
    edge_labels = ["$(first(g.edges[Pair(e)]))" for e ∈ edges(g.graph)]
    
    gplot(
        g.graph, Vector{Float64}(r_id), Vector{Float64}(s_id);
        nodelabel=vertex_labels,
        edgelabel=edge_labels
    )
end

function GraphPlot.gplot(m::Union{TMA,MA})
    g = to_graph(m)

    states_u = unique(g.vertices)
    s_id = indexin(g.vertices,states_u)
    
    vertex_labels = ["S$(s)" for s ∈ s_id]
    edge_labels = ["$(first(g.edges[Pair(e)]))" for e ∈ edges(g.graph)]
    
    gplot(
        g.graph;
        nodelabel=vertex_labels,
        edgelabel=edge_labels
    )
end



end