abstract type Automaton{ST,ET} end

struct UntimedAutomaton{ST,ET,TT} <: Automaton{ST,ET}
    X::Set{ST}
    E::Set{ET}
    T::Dict{Tuple{S,E}, TT}
    x0::ST
end

const TimedEdge{ST} = NamedTuple{(:guard,:reset,:new_state), Tuple{ClockConstraint,Set{Int},ST}}

struct TimedAutomaton{ST,ET,TT} <: Automaton{ST,ET}
    X::Set{ST}
    E::Set{ET}
    T::Dict{Tuple{S,E}, TimedEdge}
    I::ClockConstraint
    x0::ST
    num_clocks::Int
end