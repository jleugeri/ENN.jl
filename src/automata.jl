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

struct TAArc{State,T}
    source::State
    target::State
    guard::TAConstraint{T}
    message::TAMessage
    resets::Set{Symbol}
end

struct TA{State, T}
    states::Vector{State}                    # set of states
    initial_state::State                           # initial state
    clocks::Vector{Symbol}                   # set of clocks - the numeric index is used internally
    symbols::Set{Symbol}                         # alphabet of symbols
    arcs::Set{TAArc{State,T}}       # arcs of the TA
    arcs_by_state::DefaultDict{State, Set{TAArc{State,T}}}
    invariants::DefaultDict{State, TAConstraint{T}}   # invariants of the TA's states
    function TA(states::Vector{State},initial_state::State,clocks::Vector{Symbol},symbols::Set{Symbol},arcs::Set{TAArc{State,T}},invariants::Dict{State, TAConstraint{T}}) where {State, T}
        # make sure all arcs have the same clock-set
        arcs_by_state = DefaultDict{State, Set{TAArc{State,T}}}(()->Set{TAArc{State,T}}())

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

# Timed Transition System
struct TTS{State, T}
    states::Vector{Tuple{State,TAConstraint{T}}}
    transitions::Vector{Pair{Tuple{State,TAConstraint{T}}, Tuple{TAMessage, Tuple{State,TAConstraint{T}}}}}
    function TTS(ta::TA{State,T}; kwargs...) where {State, T} 
        new{State,T}(TimedAutomata.zone_graph(ta; kwargs...)...)
    end
end
