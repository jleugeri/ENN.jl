module TimedAutomata

export TMAState, ClockAndCondition, ClockClockRefCondition, ClockCondition, ClockConstRefCondition, ClockNegCondition, ClockNoCondition, TMA, TMARun, TMAState, TMATransition, RunState

const Event = NamedTuple{(:time,:symbol), Tuple{Float64, Symbol}}

Base.convert(::Type{Event}, ev::Tuple) = (time=ev[1], symbol=ev[2])

const TMAState = Vector
mutable struct RunState
    state::TMAState
    clocks::Vector{Float64}
end

Base.:(==)(s1::RunState, s2::RunState) = s1.state==s2.state && s1.clocks==s2.clocks

abstract type ClockCondition end


struct ClockConstRefCondition{OP} <: ClockCondition
    clock_id::Int
    ref::Float64
end

struct ClockClockRefCondition{OP} <: ClockCondition
    clock_id::Int
    ref_id::Int
end

struct ClockNegCondition <: ClockCondition
    cond::ClockCondition
end

struct ClockAndCondition <: ClockCondition
    cond1::ClockCondition
    cond2::ClockCondition
end

struct ClockNoCondition <: ClockCondition end

(c::ClockNoCondition)(x) = true

@generated (c::ClockConstRefCondition{OP})(clocks) where OP = quote
        return $(OP)(clocks[c.clock_id],c.ref)
end

@generated (c::ClockClockRefCondition{OP})(clocks) where OP = quote
    return $(OP)(clocks[c.clock_id],clocks[c.ref_id])
end

(c::ClockNegCondition)(clocks) = ~c.cond(clocks)

(c::ClockAndCondition)(clocks) = c.cond1(clocks) && c.cond2(clocks)


struct TMATransition
    s::TMAState
    s′::TMAState
    a::Symbol
    λ::Vector{Int}
    δ::ClockCondition
end

struct TMA
    Σ::Vector{TMAState}
    Σ′::Vector{TMAState}
    num_clocks::Int
    E::Vector{TMATransition}
    s0::TMAState
end

struct TMARun
    tma::TMA
    inputs::Vector{Event}
    X::RunState
end
TMARun(tma,inputs)=TMARun(tma,inputs,RunState(tma.s0, zeros(Float64,tma.num_clocks)))

function Base.iterate(r::TMARun, t0i=(0.0,0))
    t0,i = t0i

    if i==0
        return r.X,(0.0, 1) 
    end

    ret=iterate(r.inputs,i)
    
    if ret == nothing
        return nothing
    end

    inp,i = ret
    a=inp.symbol
    t1=inp.time

    dt=t1-t0
    r.X.clocks .+= dt
    for e in r.tma.E
        if a==e.a && r.X.state == e.s && e.δ(r.X.clocks)
            r.X.clocks[e.λ] .= 0
            r.X.state = e.s′
            break
        end
    end
    return RunState(r.X.state,r.X.clocks),(t1,i)
end

Base.length(r::TMARun) = length(r.inputs)+1

end