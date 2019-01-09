module TimedAutomata
using LightGraphs, Base.Iterators, IntervalArithmetic
export TimeSlice,TimePoint,TimeInterval,TimeLeftHalfLine,TimeRightHalfLine, ClockRegion, implies, collect_constants, @Clk_str,
    MA,TMA,MATransition,TMATransition, to_graph, untime, language
#, TMARun, RunState

@enum ORDER_RELATIONSHIP order_less=1 order_less_equal=2 order_equal=3 order_larger=4 order_larger_equal=5

#const TimeSlice{T} = Union{Nothing,T,Tuple{T,T},Tuple{T,Nothing}}

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


mutable struct RunState{TMAState}
    state::TMAState
    clocks::Vector{Float64}
end

Base.:(==)(s1::RunState{TMAState}, s2::RunState{TMAState}) where {TMAState} = all(s1.state .== s2.state) && all(s1.clocks .≈ s2.clocks)


##########################
#    Untimed Automata    #
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

function untime(m::TMA{TMAState}) where TMAState
    num_clocks = m.num_clocks
    
    clock_constants = collect_constants(m)

    clock_scaler = lcm(denominator.(reduce(union,clock_constants))...)//gcd(numerator.(reduce(union,clock_constants))...)
    if denominator(clock_scaler)==0
        clock_scaler = 0//1
    end

    c_max = Int.(maximum.(clock_constants) .* clock_scaler)

    r_max = 2 .* (c_max .+ 1)

    function condition_from_region((r,o))
        ts = Vector{TimeSlice}(undef, num_clocks)
        for clk ∈ eachindex(r)
            (c,i) = divrem(r[clk]-1, 2)
            ts[clk] = i==0 ? TimePoint{Int}(c) : (c<c_max[clk] ? TimeInterval{Int,true,true}(c,c+1) : TimeRightHalfLine{Int,true}(c))
        end
        return ClockRegion(Int,ts,o)
    end

    function update_clocks(state, increase=ones(Bool, num_clocks), reset=zeros(Bool,num_clocks))
        (r,o) = state
        r_new = copy(r)
        o_new = copy(o)

        # increase integral clock part
        r_new[increase] .+= 1

        # reset integral clock part
        r_new[reset] .= 1

        # clip clocks to unbounded region
        idx = r_new .> r_max
        r_new[idx] .= r_max[idx]

        # if the clock region is a constant or out of bounds, treat it as exact and ignore order
        exact = (r_new .% 2 .== 1)
        is_out_of_bound = r_new .== r_max
        unordered = exact .| is_out_of_bound

        # update the order of all of the clocks' fractional parts that are not exactly 0
        new_ranks = sort!(unique(o_new[.~(unordered)]))
        o_new[.~(unordered)] .= indexin(o_new[.~(unordered)], new_ranks)
        o_new[unordered] .= 0

        
        return r_new, o_new
        
    end

    # construct region automaton
    states = Tuple{TMAState, Tuple{Vector{Int},Vector{Int}}}[]
    accepting_states = Tuple{TMAState, Tuple{Vector{Int},Vector{Int}}}[]
    edges = Dict{Tuple{Tuple{TMAState, Tuple{Vector{Int},Vector{Int}}},Symbol}, Vector{MATransition{Tuple{TMAState, Tuple{Vector{Int},Vector{Int}}}}}}()

    s0 = (m.s0,(ones(Int, num_clocks),zeros(Int, num_clocks)))

    queue = [s0]
    while ~isempty(queue)
        (s,r) = q = popfirst!(queue)
        # If we took care of the new state before, move on
        if q ∈ states
            continue
        end

        is_exact = (r[1] .% 2) .== 1
        largest_rank = maximum(r[2])
        is_largest_rank = r[2] .== largest_rank
        is_out_of_bound = r[1] .>= r_max

        cnd = condition_from_region(r)
        # go through all symbol transitions
        for ((s′,a),ee) ∈ pairs(m.E)
            if s′ == s
                for e ∈ ee
                    # if there is an symbol transition, reset the corresponding clocks (both integral and fractional part)
                    # check if edge constraint holds within the region
                    if implies(cnd,e.δ)
                        tmp = zeros(Bool, num_clocks)
                        tmp[e.λ] .= true
                        r′ = update_clocks(r, zeros(Bool,num_clocks),tmp)
                        q1 = (e.s, r′)
                        el = get(edges, (q,a), MATransition{Tuple{TMAState, Tuple{Vector{Int},Vector{Int}}}}[])
                        push!(el, MATransition(q1))
                        edges[(q,a)] = el

                        # If we didn't take care of the new state before, put it in the queue
                        if q1 ∉ states
                            push!(queue, q1)
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

            rsucc = update_clocks(r, increase_clocks)
            # time successor is only triggered if invariants allow it
            if s ∉ keys(m.I) || implies(condition_from_region(rsucc),m.I[s])
                q1=(s,rsucc)
                edges[(q,:τ)] = [MATransition(q1)]
                push!(queue, q1)
            end
        end

        # add state to the list
        push!(states, q)
        # if the state is an accepting state, add it
        if s ∈ m.Σ′
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

function to_graph(ma::MA{MAState}) where MAState
    edges = Dict{Pair{Int,Int},Symbol}()
    states_u = collect(ma.Σ)
    l = LightGraphs.SimpleDiGraph(length(states_u))

    for ((s,a),ee) ∈ pairs(ma.E)
        for e ∈ ee
            s_idx = findfirst(x->x==s,states_u)
            t_idx = findfirst(x->x==e.s,states_u)
            add_edge!(l,s_idx , t_idx)
            edges[s_idx=>t_idx] = a
        end
    end

    return (graph=l, edges=edges, vertices=states_u)
end


function to_graph(tma::TMA{TMAState}) where TMAState
    edges = Dict{Pair{Int,Int},Tuple{Symbol,ClockRegion}}()
    states_u = collect(deepcopy.(tma.Σ))
    l = LightGraphs.SimpleDiGraph(length(states_u))

    for ((s,a),ee) ∈ pairs(tma.E)
        for e ∈ ee
            s_idx = findfirst(x->x==s,states_u)
            t_idx = findfirst(x->x==e.s,states_u)
            add_edge!(l,s_idx , t_idx)
            edges[s_idx=>t_idx] = (a,e.δ)
        end
    end

    T = typeof(tma.I).parameters[2].parameters[2]
    return (graph=l, edges=edges, vertices=(x->(x,get(tma.I, x, ClockRegion{tma.num_clocks,T}()))).(states_u))
end

end