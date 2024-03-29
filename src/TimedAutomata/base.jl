export TAMessage, TAArc, TA, TTS, MSG_DIRECTION, MSG_IN, MSG_OUT, MSG_SILENT

@enum MSG_DIRECTION MSG_IN=1 MSG_OUT=0 MSG_SILENT=-1

struct TAMessage
    direction::MSG_DIRECTION
    symbol::Symbol
end
TAMessage() = TAMessage(MSG_SILENT, Symbol("—"))
function Base.show(io::IO, m::TAMessage)
    dir = if m.direction==MSG_IN
        "?"
    elseif m.direction==MSG_OUT
        "!"
    else
        ""
    end
    print(io, "$(m.symbol)$(dir)")
end


Base.iterate(x::TAMessage) = (x, nothing)
Base.iterate(::TAMessage, ::Any) = nothing
Base.length(x::TAMessage) = 1
Base.Broadcast.broadcastable(c::TAMessage) = Ref(c)
Base.adjoint(c::TAMessage) = c
Base.copy(m::TAMessage) = TAMessage(m.direction,m.symbol)
Base.deepcopy(m::TAMessage) = TAMessage(deepcopy(m.direction),deepcopy(m.symbol))


mutable struct TAArc{State,T}
    source::State
    target::State
    guard::TAConstraint{T}
    message::TAMessage
    resets::Set{Symbol}
end

Base.iterate(x::TAArc) = (x, nothing)
Base.iterate(::TAArc, ::Any) = nothing
Base.length(x::TAArc) = 1
Base.Broadcast.broadcastable(c::TAArc) = Ref(c)
Base.adjoint(c::TAArc) = c
Base.copy(a::T) where T<:TAArc = T(a.source,a.target,a.guard,a.message,a.resets)
Base.deepcopy(a::T) where T<:TAArc = T(deepcopy(a.source),deepcopy(a.target),deepcopy(a.guard),deepcopy(a.message),deepcopy(a.resets))
function Base.show(io::IO, arc::TAArc)
    m=repr(arc.message)
    r=join(arc.resets, ",")
    r=isempty(r) ? "" : "{$(r)}"
    g=repr(arc.guard)
    i=join(filter(!isempty, (m,g,r)),",")
    i = isempty(i) ? i : "($(i))"
    print(io, "$(repr(arc.source))——$(i)—⟶$(repr(arc.target))")
end

function set_clock_list!(arc::TAArc{State,T}, clocks::Vector{Symbol}; mapping=identity) where {State,T}
    set_clock_list!(arc.guard, clocks; mapping)
    new_resets = mapping.(arc.resets)
    empty!(arc.resets)
    union!(arc.resets, new_resets)
    return arc
end

struct TAState{CombinedState}
    components::CombinedState
    function TAState(arg::T, args...) where T
        s = (arg,args...)
        new{typeof(s)}(s)
    end
end
Base.:(==)(s::TAState, o::TAState) = s.components==o.components
Base.hash(s::TAState) = hash(s.components)
Base.iterate(x::TAState) = iterate(x.components)
Base.iterate(x::TAState, s) = iterate(x.components, s)
Base.length(x::TAState) = length(x.components)
Base.lastindex(x::TAState) = Base.lastindex(x.components)
Base.getindex(x::TAState, s) = x.components[s]
Base.Broadcast.broadcastable(x::TAState) = x.components
Base.adjoint(x::TAState) = x
Base.copy(x::T) where T<:TAState = T(x.components)
Base.deepcopy(x::T) where T<:TAState = TAState(deepcopy.(x.components)...)
Base.:|(args::TAState...) = TAState(flatten(x.components for x in args)...)
function Base.show(io::IO, x::TAState)
    s=join(repr.(x.components),"|")
    print(io, "⟨$(s)⟩")
end

struct TA{State, T}
    states::Vector{State}                    # set of states
    initial_state::State                           # initial state
    clocks::Vector{Symbol}                   # set of clocks - the numeric index is used internally
    symbols::Set{Symbol}                         # alphabet of symbols
    arcs::Vector{TAArc{State,T}}       # arcs of the TA
    arcs_by_state::DefaultDict{State, Vector{TAArc{State,T}}}
    invariants::DefaultDict{State, TAConstraint{T}}   # invariants of the TA's states
    function TA(states::Vector{State},initial_state::State,clocks::Vector{Symbol},symbols::Set{Symbol},arcs::Vector{TAArc{State,T}},invariants::Dict{State, TAConstraint{T}}) where {State<:TAState, T}
        # make sure all arcs have the same clock-set
        arcs_by_state = DefaultDict{State, Vector{TAArc{State,T}}}(()->Vector{TAArc{State,T}}())

        for arc ∈ arcs
            source_id = findfirst(==(arc.source), states)
            @assert !isnothing(source_id) "Source state '$(arc.source)' of arc '$(arc)' not in set of states."
            target_id = findfirst(==(arc.target), states)
            @assert !isnothing(target_id) "Target state '$(arc.target)' of arc '$(arc)' not in set of states."

            # make sure the source/target actually refer to the exact object in the TA's state-vector
            arc.source = states[source_id]
            arc.target = states[target_id]

            set_clock_list!(arc.guard, clocks)
            for reset ∈ arc.resets
                @assert reset ∈ clocks "Clock '$(reset)' should be reset by arc '$(arc)', but no such clock is known! Clock list: '$(clocks)'."
            end
            push!(arcs_by_state[arc.source], arc)
        end

        # make sure all guards have the same clock-set, make sure there is an (empty?) invariance for each state
        default_guard = @constraint( true, T )
        set_clock_list!(default_guard, clocks)

        new_invariants = DefaultDict(default_guard, invariants)
        for (state,I) ∈ pairs(invariants)
            set_clock_list!(I, clocks)
            state_id = findfirst(==(state), states)
            @assert !isnothing(state_id) "State '$(state)' of guard expression '$(I)' not in the set states."
            new_invariants[states[state_id]] = I
        end

        new{State, T}(states,initial_state,clocks,symbols,arcs,arcs_by_state,new_invariants)
    end
end
TA(states::Vector,initial_state::State,clocks::Vector{Symbol},symbols::Set{Symbol},arcs::Vector{TAArc{State,T}},invariants::Dict{State, TAConstraint{T}}) where {State,T} =
    TA(TAState.(states),TAState(initial_state),clocks, symbols, [TAArc(TAState(a.source),TAState(a.target),a.guard,a.message,a.resets) for a in arcs], Dict(TAState(k)=>v for (k,v) in pairs(invariants)))

Base.copy(ta::T) where T<:TA = TA(ta.states,ta.initial_state,ta.clocks,ta.symbols,ta.arcs,Dict(ta.invariants))
Base.deepcopy(ta::T) where T<:TA = TA(deepcopy(ta.states),deepcopy(ta.initial_state),deepcopy(ta.clocks),deepcopy(ta.symbols),deepcopy(ta.arcs),Dict(deepcopy(ta.invariants)))
Base.iterate(x::TA) = (x, nothing)
Base.iterate(::TA, ::Any) = nothing
Base.length(x::TA) = 1
Base.Broadcast.broadcastable(c::TA) = Ref(c)
Base.adjoint(c::TA) = c

function set_clock_list!(ta::TA{State,T}, clocks::Vector{Symbol}; mapping=identity) where {State,T}
    clocks=copy(clocks)
    # update clock for each state invariant
    for state in ta.states
        set_clock_list!(ta.invariants[state], clocks; mapping)
    end

    # update clock for each arc guard & reset
    #should also apply changes to arcs_by_state, because they share the same objects
    for arc in ta.arcs
        set_clock_list!(arc, clocks; mapping)
    end

    empty!(ta.clocks)
    append!(ta.clocks,clocks)

    return ta
end

# Timed Transition System
struct TTS{State, T}
    states::Vector{Tuple{State,TAConstraint{T}}}
    transitions::Vector{Pair{Tuple{State,TAConstraint{T}}, Tuple{TAMessage, Tuple{State,TAConstraint{T}}}}}
    function TTS(ta::TA{State,T}, args...; kwargs...) where {State<:TAState, T} 
        new{State,T}(TimedAutomata.zone_graph_fwd(ta, args...; kwargs...)...)
    end
end

Base.iterate(x::TTS) = (x, nothing)
Base.iterate(::TTS, ::Any) = nothing
Base.length(x::TTS) = 1
Base.Broadcast.broadcastable(c::TTS) = Ref(c)
Base.adjoint(c::TTS) = c
Base.copy(tts::T) where T<:TTS = T(tts.states,tts.transitions)
Base.deepcopy(tts::T) where T<:TTS = T(deepcopy(tts.states),deepcopy(tts.transitions))
