module TimedAutomata

export ClockAndCondition, ClockClockRefCondition, ClockCondition, ClockConstRefCondition, ClockNegCondition, ClockNoCondition, TMA, TMARun, TMATransition, RunState, parse_condition, @condition

const Event = Tuple{Float64, Symbol}

mutable struct RunState{TMAState}
    state::TMAState
    clocks::Vector{Float64}
end

Base.:(==)(s1::RunState{TMAState}, s2::RunState{TMAState}) where {TMAState} = all(s1.state .== s2.state) && all(s1.clocks .≈ s2.clocks)

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


struct TMATransition{TMAState}
    s::TMAState
    s′::TMAState
    a::Symbol
    λ::Vector{Int}
    δ::ClockCondition
end

struct TMA{TMAState}
    Σ::Vector{TMAState}
    Σ′::Vector{TMAState}
    num_clocks::Int
    E::Vector{TMATransition{TMAState}}
    s0::TMAState
end

struct TMARun{TMAState}
    tma::TMA{TMAState}
    inputs::Vector{Event}
    X::RunState{TMAState}
end
TMARun(tma,inputs)=TMARun(tma,inputs,RunState(deepcopy(tma.s0), zeros(Float64,tma.num_clocks)))


function parse_condition(expr)
    if isa(expr, Bool) && expr
        ClockNoCondition()
    elseif isa(expr, Expr) && expr.head == :call
        if expr.args[1] == :~
            ClockNegCondition(parse_condition(expr.args[2]))
        elseif expr.args[1] == :<
            if isa(expr.args[2],Expr) && isa(expr.args[3],Expr) && expr.args[2].head == :ref && expr.args[3].head == :ref && expr.args[2].args[1] == :c &&  expr.args[3].args[1] == :c
                ClockClockRefCondition{:<}(expr.args[2].args[2], expr.args[3].args[2])
            elseif isa(expr.args[2],Expr) && expr.args[2].head == :ref && isa(expr.args[3], Real) && expr.args[2].args[1] == :c
                ClockConstRefCondition{:<}(expr.args[2].args[2], expr.args[3])
            elseif isa(expr.args[3],Expr) && expr.args[3].head == :ref && isa(expr.args[2], Real) && expr.args[3].args[1] == :c
                ClockConstRefCondition{:>}(expr.args[3].args[2], expr.args[2])
            else
                raise("Either LHS or RHS must be reference to clock")
            end
        elseif expr.args[1] == :>
            if isa(expr.args[2],Expr) && isa(expr.args[3],Expr) && expr.args[2].head == :ref && expr.args[3].head == :ref && expr.args[2].args[1] == :c &&  expr.args[3].args[1] == :c
                ClockClockRefCondition{:>}(expr.args[2].args[2], expr.args[3].args[2])
            elseif isa(expr.args[2],Expr) && expr.args[2].head == :ref && isa(expr.args[3], Real) && expr.args[2].args[1] == :c
                ClockConstRefCondition{:>}(expr.args[2].args[2], expr.args[3])
            elseif isa(expr.args[3],Expr) && expr.args[3].head == :ref && isa(expr.args[2], Real) && expr.args[3].args[1] == :c
                ClockConstRefCondition{:<}(expr.args[3].args[2], expr.args[2])
            else
                raise("Either LHS or RHS must be reference to clock")
            end
        #else
        end
    elseif isa(expr, Expr) && expr.head == :(&&)
        c1 = parse_condition(expr.args[1])
        c2 = parse_condition(expr.args[2])

        ClockAndCondition(c1, c2)
    #else
    end
end


macro condition(expr)
    return parse_condition(expr)
end

function Base.iterate(r::TMARun, t0i=(0.0,0))
    t0,i = t0i

    if i==0
        r.X.state = deepcopy(r.tma.s0)
        r.X.clocks = zeros(Float64,r.tma.num_clocks)

        return (t0,deepcopy(r.X)),(0.0, 1) 
    end

    ret=iterate(r.inputs,i)
    
    if ret == nothing
        return nothing
    end

    (t1,a),i = ret

    dt=t1-t0
    r.X.clocks .+= dt
    for e in r.tma.E
        if a==e.a && r.X.state == e.s && e.δ(r.X.clocks)
            r.X.clocks[e.λ] .= 0
            r.X.state = deepcopy(e.s′)
            break
        end
    end
    return (t1,deepcopy(r.X)),(t1,i)
end

Base.length(r::TMARun) = length(r.inputs)+1

end