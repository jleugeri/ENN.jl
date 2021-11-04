export TAMessage, TAArc, TA, TTS, MSG_DIRECTION, MSG_IN, MSG_OUT

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


struct TAArc{State,T}
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

struct TA{State, T}
    states::Vector{State}                    # set of states
    initial_state::State                           # initial state
    clocks::Vector{Symbol}                   # set of clocks - the numeric index is used internally
    symbols::Set{Symbol}                         # alphabet of symbols
    arcs::Vector{TAArc{State,T}}       # arcs of the TA
    arcs_by_state::DefaultDict{State, Vector{TAArc{State,T}}}
    invariants::DefaultDict{State, TAConstraint{T}}   # invariants of the TA's states
    function TA(states::Vector{State},initial_state::State,clocks::Vector{Symbol},symbols::Set{Symbol},arcs::Vector{TAArc{State,T}},invariants::Dict{State, TAConstraint{T}}) where {State, T}
        # make sure all arcs have the same clock-set
        arcs_by_state = DefaultDict{State, Vector{TAArc{State,T}}}(()->Vector{TAArc{State,T}}())

        for arc ∈ arcs
            set_clock_list!(arc.guard, clocks)
            for reset ∈ arc.resets
                @assert reset ∈ clocks "Clock '$(reset)' should be reset by arc '$(arc)', but no such clock is known! Clock list: '$(clocks)'."
            end
            push!(arcs_by_state[arc.source], arc)
        end

        # make sure all guards have the same clock-set, make sure there is an (empty?) invariance for each state
        default_guard = @constraint( true, T )
        set_clock_list!(default_guard, clocks)

        invariants = DefaultDict(default_guard, invariants)
        for (state,I) ∈ pairs(invariants)
            set_clock_list!(I, clocks)
            @assert state ∈ states "State '$(state)' of guard expression '$(I)' is not in the list of the TA's states ('$(states)')."
        end

        new{State, T}(states,initial_state,clocks,symbols,arcs,arcs_by_state,invariants)
    end
end

Base.copy(ta::T) where T<:TA = TA(ta.states,ta.initial_state,ta.clocks,ta.symbols,ta.arcs,Dict(ta.invariants))
Base.deepcopy(ta::T) where T<:TA = TA(deepcopy(ta.states),deepcopy(ta.initial_state),deepcopy(ta.clocks),deepcopy(ta.symbols),deepcopy(ta.arcs),Dict(deepcopy(ta.invariants)))
Base.iterate(x::TA) = (x, nothing)
Base.iterate(::TA, ::Any) = nothing
Base.length(x::TA) = 1
Base.Broadcast.broadcastable(c::TA) = Ref(c)
Base.adjoint(c::TA) = c

# Timed Transition System
struct TTS{State, T}
    states::Vector{Tuple{State,TAConstraint{T}}}
    transitions::Vector{Pair{Tuple{State,TAConstraint{T}}, Tuple{TAMessage, Tuple{State,TAConstraint{T}}}}}
    function TTS(ta::TA{State,T}; kwargs...) where {State, T} 
        new{State,T}(TimedAutomata.zone_graph_fwd(ta; kwargs...)...)
    end
end

Base.iterate(x::TTS) = (x, nothing)
Base.iterate(::TTS, ::Any) = nothing
Base.length(x::TTS) = 1
Base.Broadcast.broadcastable(c::TTS) = Ref(c)
Base.adjoint(c::TTS) = c
Base.copy(tts::T) where T<:TTS = T(tts.states,tts.transitions)
Base.deepcopy(tts::T) where T<:TTS = T(deepcopy(tts.states),deepcopy(tts.transitions))
