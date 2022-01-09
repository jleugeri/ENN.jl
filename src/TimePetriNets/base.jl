using SparseArrays

export TPNState, TPN, isenabled, isready, isvalidstate, update_clock!, max_delay, elapse!, fire!, fire_all!, fire_necessary!, run!


"""
    TPNState

Class to hold the state (=place and transition markings) of the Time Petri Net.
"""
struct TPNState{M,H}
    m::Vector{M}                # p-marking
    h::Vector{Union{Nothing,H}} # t-marking
end

# forward declare function so it can be used in constructor
function isenabled end
function update! end

"""
    TPN{M,H}

Class to hold the entire Time Petri Net.
This implementation uses an unusual parameterization inspired by [*Time and Petri Nets*, 2013][1]:

# Attributes:
- `P::Vector{Symbol}`: names of the places
- `T::Vector{Symbol}`: names of the transitions
- `R::SparseMatrixCSC{UInt}`: additional read arcs of the transitions
- `ΔF::SparseMatrixCSC{Int}`: token-changes caused by the transitions (one column per transition)
- `eft::Vector{H}`: earliest firing time of transition
- `lft::Vector{H}`: last firing time of transition
- `x₀::TPNState{M,H}`: initial state vector

[1]: https://link.springer.com/book/10.1007/978-3-642-41115-1
"""
struct TPN{M,H}
    P::Vector{Symbol}           # names of the places
    T::Vector{Symbol}           # names of the transitions
    C::SparseMatrixCSC{M}       # all dependencies of the transitions (read arcs + regular arcs)
    ΔF::SparseMatrixCSC{M}      # token-changes caused by the transitions (one column per transition)
    eft::Vector{H}              # earliest firing time of transition
    lft::Vector{H}              # last firing time of transition
    x₀::TPNState{M,H}           # initial state vector
    # internal helper variables
    may_disable::SparseMatrixCSC{Bool}
    may_enable::SparseMatrixCSC{Bool}

    """
        TPN(P,T,C,ΔF,eft,lft,m₀::Union{TPNState,Vector{<:Integer}})

    Construct the `TPN`.

    # Arguments:
    - `P`: names of the places
    - `T`: names of the transitions
    - `R`: additional read arcs of the transitions
    - `ΔF`: token-changes caused by the transitions (one column per transition)
    - `eft`: earliest firing time of transition
    - `lft`: last firing time of transition
    - `m₀::Union{TPNState,Vector{<:Integer}}`: initial state vector (if type `TPNState`) or initial marking
    """
    function TPN(P,T,R,ΔF,eft,lft,m₀::Union{TPNState,Vector{<:Integer}})
        M = eltype(ΔF)
        H = eltype(eft)

        C = sparse(R - ΔF.*(ΔF .< 0))
        nT = length(T)
        # initial state or initial marking provided
        x₀ = isa(m₀, TPNState) ? m₀ : TPNState(Vector{M}(m₀), Vector{Union{Nothing,H}}(nothing, nT))
    
        may_disable = spzeros(Bool, nT, nT)
        may_enable = spzeros(Bool, nT, nT)
        
        pn = new{M,H}(P,T,C,ΔF,eft,lft,x₀, may_disable, may_enable)
        
        update!(pn)
        pn
    end
end

"""
    isenabled(pn::TPN, x, t_id)

Returns `true` if the transition with ID `t_id` of the TPN `pn` is enabled in state `x` and false otherwise.

!!! info
    The notion of a transition being `enabled` in a Time Petri Net only refers to 
    the place marking (i.e. tokens) and does **not** check the transition markings (i.e. clocks)!

    For checking if the transition can be actually fired, see [`isready`](@ref)
"""
Base.@propagate_inbounds function isenabled(pn::TPN, x, t_id)
    @boundscheck @assert 1 ≤ t_id ≤ length(pn.T) "Transition $(t_id) out of bounds!";
    rows = rowvals(pn.C)
    vals = nonzeros(pn.C)
    @inbounds for i in nzrange(pn.C, t_id)
        if x.m[rows[i]] < vals[i]
            return false
        end
    end
    return true
end

"""
    isready(pn::TPN, x, t_id)

Returns `true` if the transition with ID `t_id` of the TPN `pn` is ready in state `x` and false otherwise.

!!! info
    The notion of a transition being `ready` in a Time Petri Net requires both
    the place marking (i.e. tokens) and the transition markings (i.e. clocks) to allow firing!

    For only checking the place markings, see [`isenabled`](@ref)
"""
Base.@propagate_inbounds function isready(pn::TPN, x, t_id)
    return !isnothing(x.h[t_id]) && pn.eft[t_id] ≤ x.h[t_id] ≤ pn.lft[t_id]
end

"""
    isvalidstate(pn::TPN, x)

Returns `true` if the place and transition markings of the state `x` in TPN `pn` are consistent, and `false` otherwise.
"""
function isvalidstate(pn::TPN, x)
    for (t_id,h_t) in enumerate(x.h)
        en = isenabled(pn, x, t_id)
        @inbounds if (en && h_t > lft[t_id]) || (!en && !isnothing(h_t))
            return false
        end
    end
    return true
end

"""
    update_clock!(pn::TPN, x::TPNState{M,H}, t_id) where {M,H}

Sets the clock of the TPN `pn`'s transition `t_id` in state `x` to
`zero(H)` if the transition was just enabled 
or to `nothing` if the transition was just disabled

!!! warning
    This modifies the state `x`!
"""
Base.@propagate_inbounds function update_clock!(pn::TPN, x::TPNState{M,H}, t_id) where {M,H}
    en = isenabled(pn, x, t_id)
    @inbounds val = x.h[t_id]
    @inbounds x.h[t_id] = if en && isnothing(val)
        zero(H)
    elseif en
        val
    else
        nothing
    end
end

"""
    update!(pn::TPN)

Updates the TPN after structural changes have been made (e.g. adding arcs).
"""
function update!(pn::TPN)
    # check which transitions can disable each other
    ΔF_rows = rowvals(pn.ΔF)
    ΔF_vals = nonzeros(pn.ΔF)
    C_rows = rowvals(pn.C)
    @inbounds for t_id in eachindex(pn.T)
        # if the number of tokens is decreased, may disable other transition
        decreased = [ΔF_rows[i] for i in nzrange(pn.ΔF, t_id) if ΔF_vals[i] < 0]
        # if the number of tokens is increased, may enable other transition
        increased = [ΔF_rows[i] for i in nzrange(pn.ΔF, t_id) if ΔF_vals[i] > 0]

        for t_id′ in eachindex(pn.T)
            checked = C_rows[nzrange(pn.C, t_id′)]
            if !isdisjoint(decreased, checked) 
                pn.may_disable[t_id′,t_id] = true 
            end
            if !isdisjoint(increased, checked) 
                pn.may_enable[t_id′,t_id] = true 
            end
        end
    end

    # set clocks of all enabled transitions in initial state to zero
    @inbounds for t_id in eachindex(pn.T)
        if isenabled(pn, pn.x₀, t_id)
            pn.x₀.h[t_id] = zero(H)
        end
    end
    return nothing
end


"""
    fire!(pn::TPN, x::TPNState, t_id; check_ready=true)

Fires the transition `x` of TPN `pn` from state `x`.

If `check_ready` is `true` (the default), then ensure that the transition is ready before firing. 

!!! warning
    This modifies the state `x`!
"""
Base.@propagate_inbounds function fire!(pn::TPN, x::TPNState, t_id; check_ready=true)
    if !check_ready || isready(pn, x, t_id)
        rows = rowvals(pn.ΔF)
        vals = nonzeros(pn.ΔF)
        @inbounds for i in nzrange(pn.ΔF, t_id)
            @inbounds x.m[rows[i]] += vals[i]
        end
 
        @inbounds for t_id′ in rowvals(pn.may_disable)[nzrange(pn.may_disable,t_id)] ∪ rowvals(pn.may_enable)[nzrange(pn.may_enable,t_id)]
            update_clock!(pn,x,t_id′)
        end 
    else
        error("Transition #$(t_id) not ready!")
    end
    nothing
end

"""
    fire_all!(pn::TPN, x; must_fire=isready)

Fires all transitions `t` of TPN `pn` in state `x` for which `must_fire(pn, x, t)` returns `true`.

Returns a vector of the fired transition labels.

!!! warning
    This modifies the state `x`!
"""
function fire_all!(pn::TPN, x; must_fire=isready)
    options = Set(t_id for t_id in eachindex(pn.T) if must_fire(pn, x, t_id))
    transitions = Symbol[]
    while !isempty(options)
        t_id = pop!(options, rand(options))
        fire!(pn, x, t_id; check_ready=false)
        push!(transitions, pn.T[t_id])

        # check if transition must fire again
        if must_fire(pn, x, t_id)
            push!(options, t_id)
        end

        # check transitions that might have been disabled
        @inbounds for t_id′ in rowvals(pn.may_disable)[nzrange(pn.may_disable,t_id)]
            if !isready(pn, x, t_id′)
                delete!(options, t_id′)
            end
        end

        # check transitions that might have been enabled
        @inbounds for t_id′ in rowvals(pn.may_enable)[nzrange(pn.may_enable,t_id)]
            if must_fire(pn, x, t_id′)
                push!(options, t_id′)
            end
        end
    end
    return transitions
end

"""
    fire_necessary!(pn::TPN, x::TPNState) 

Fires all transitions `t` of TPN `pn` in state `x` that are ready and cannot be delayed any further.

Returns a vector of the fired transition labels.

!!! warning
    This modifies the state `x`!
"""
function fire_necessary!(pn::TPN, x::TPNState) 
    isnecessary(pn, x, t_id) = 
        @inbounds isready(pn, x, t_id) && x.h[t_id] == pn.lft[t_id]
    fire_all!(pn, x; must_fire=isnecessary)
end

"""
    max_delay(pn, x::TPNState{TT}) where TT

Computes the maximum delay possible of TPN `pn` in state `x` that does not invalidate the state.
"""
function max_delay(pn, x::TPNState{TT}) where TT
    m = typemax(TT)
    for (tᵢ,lftᵢ) in zip(x.h,pn.lft)
        if !isnothing(tᵢ)
            m = min(m, lftᵢ-tᵢ)
        end
    end
    return m
end

"""
    elapse!(pn::TPN, x::TPNState, τ; check_valid=true)

Elapses the time `τ`, i.e. advances all clocks of TPN `pn` in state `x` by `τ`.
If `check_valid` is `true` (the default), then ensure that `τ` does not exceed the maximum allowed delay.

!!! warning
    This modifies the state `x`!
"""
Base.@propagate_inbounds function elapse!(pn::TPN, x::TPNState, τ; check_valid=true)
    if !check_valid || τ ≤ max_delay(pn, x)
        for (t_id,tᵢ) in enumerate(x.h)
            if !isnothing(tᵢ)
                @inbounds x.h[t_id] += τ
            end
        end
    else
        error("Time delay exceeds maximum delay ($(τ) > $(max_delay(pn,x)))!")
    end
end

"""

Runs the TPN `pn` from state `x`

!!! warning
    This modifies the state `x`!
"""
function run!(pn::TPN, x, events; tmax=first(last(events)))
    @assert issorted(events, by=first) "Events must be a sequence of (time,transition) pairs with increasing time."

    t_total = 0
    for (time, ext_t_id) in events
        while t_total < time
            delay = max_delay(pn, x)
            # elapse time until the next (external or internal event)
            if time ≤ t_total + delay
                # squeeze in our externally triggered transition
                t_total = time
                elapse!(pn, x, time-t_total)
                fire!(pn, x, ext_t_id)

                # fire all other events (to be sure)
                fire_necessary!(pn, x, ext_t_id)
            else
                # proceed to next regular event
                elapse!(pn, x, delay)
                fire_necessary!(pn, x, ext_t_id)
            end

            t_new = min(time, t_total+_max_delay)
            elapse!(pn, x, t_new-t_total)
            t_total = t_new
        end
    end

    x=copy(x)
    trace = []
    states = [copy(x)]
    for t_id in t_ids
        fire!(pn, x, t_id)
        push!(trace, pn.T[t_id])
        push!(states, copy(x))
    end
    return (;trace, states)
end

"""
    |(tpns::TPN{M,H}...) where {M,H}

Combines multiple TPNs together.

This procedure:
- merges transitions by adding arcs, 
- merges earliest and last firing times by taking the more restrictive bounds,
- merges initial tokens by adding them together

"""
function Base.:|(tpns::TPN{M,H}...) where {M,H}
    all_places = union((tpn.P for tpn in tpns)...)
    all_transitions = union((tpn.T for tpn in tpns)...)
    all_R = spzeros(M, length(all_places), length(all_transitions))
    all_ΔF = spzeros(M, length(all_places), length(all_transitions))
    all_eft = zeros(H, length(all_transitions))
    all_lft = fill(typemax(H), length(all_transitions))
    all_m₀ = zeros(M, length(all_places))
    
    # merge transitions by adding arcs, 
    # merge earliest and last firing times by taking the more restrictive bounds,
    # merge initial tokens by adding them together
    for tpn in tpns
        # map the TPN's places and transitions to the new places and transitions
        p_idx = indexin(tpn.P, all_places)
        t_idx = indexin(tpn.T, all_transitions)
        R = tpn.C + tpn.ΔF .* (tpn.ΔF .< 0)
        all_R[p_idx,t_idx] += R
        all_ΔF[p_idx,t_idx] += tpn.ΔF
        all_eft[t_idx] .= max(all_eft[t_idx], tpn.eft)
        all_lft[t_idx] .= min(all_eft[t_idx], tpn.lft)
        all_m₀[p_idx] .+= tpn.x₀.m
    end

    TPN(all_places, all_transitions, all_R, all_ΔF, all_eft, all_lft, all_m₀)
end
