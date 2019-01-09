abstract type TimeSlice{T} end

struct TimeUnconstrained{T} <: TimeSlice{T} 
end

struct TimePoint{T} <: TimeSlice{T}
    val::T
end

struct TimeInterval{T,OL,OR} <: TimeSlice{T}
    val1::T
    val2::T
end

struct TimeRightHalfLine{T,OL} <: TimeSlice{T}
    val1::T
end

struct TimeLeftHalfLine{T,OR} <: TimeSlice{T}
    val2::T
end

@inline comp(::Val{true},v1,v2) = v1 <= v2
@inline comp(::Val{false},v1,v2) = v1 < v2

Base.:⊆(r1::TimeUnconstrained,r2::T) where T = false
Base.:⊆(r1::T,r2::TimeUnconstrained) where T = true

Base.:⊆(r1::TimePoint, r2::TimePoint) = r1==r2
Base.:⊆(r1::TimePoint, r2::TimeInterval{T,OL,OR}) where {T,OL,OR} = comp(Val(~OL),r2.val1,r1.val) && comp(Val(~OR),r2.val2,r1.val)
Base.:⊆(r1::TimePoint, r2::TimeLeftHalfLine{T,OR}) where {T,OR} = comp(Val(~OR),r2.val2,r1.val)
Base.:⊆(r1::TimePoint, r2::TimeRightHalfLine{T,OL}) where {T,OL} = comp(Val(~OL),r2.val1,r1.val)
Base.:⊆(r1::TimeInterval, r2::TimePoint) = false
Base.:⊆(r1::TimeLeftHalfLine, r2::TimePoint) = false
Base.:⊆(r1::TimeRightHalfLine, r2::TimePoint) = false

Base.:⊆(r1::TimeInterval{T1,OL1,OR1},r2::TimeInterval{T2,OL2,OR2}) where {T1,T2,OL1,OR1,OL2,OR2} = comp(Val(OL1 || ~OL2),r2.val1,r1.val1) && comp(Val(OR1 || ~OR2),r1.val2,r2.val2)
Base.:⊆(r1::TimeInterval{T1,OL1,OR1},r2::TimeLeftHalfLine{T2,OR2}) where {T1,T2,OL1,OR1,OR2} = comp(Val(OR1 || ~OR2),r1.val2,r2.val2)
Base.:⊆(r1::TimeInterval{T1,OL1,OR1},r2::TimeRightHalfLine{T2,OL2}) where {T1,T2,OL1,OR1,OL2} = comp(Val(OL1 || ~OL2),r2.val1,r1.val1)
Base.:⊆(r1::TimeRightHalfLine,r2::TimeInterval) = false
Base.:⊆(r1::TimeLeftHalfLine,r2::TimeInterval) = false

Base.:⊆(r1::TimeRightHalfLine,r2::TimeLeftHalfLine) = false
Base.:⊆(r1::TimeLeftHalfLine,r2::TimeRightHalfLine) = false
Base.:⊆(r1::TimeRightHalfLine{T1,OL1},r2::TimeRightHalfLine{T2,OL2}) where {T1,T2,OL1,OL2} = comp(Val(OL1 || ~OL2),r2.val1,r1.val1)
Base.:⊆(r1::TimeLeftHalfLine{T1,OR1},r2::TimeLeftHalfLine{T2,OR2}) where {T1,T2,OR1,OR2} = comp(Val(OR1 || ~OR2),r1.val2,r2.val2)


struct ClockRegion{N,T}
    slices::Vector{TimeSlice{T}}
    order_constraints::Dict{Tuple{Int,Int},ORDER_RELATIONSHIP}
end

function ClockRegion(T::Type, slices::Vector{TimeSlice},order::Vector{Int})
    num_clocks = length(order)
    
    _order(a,b) = if a<b
        order_less
    elseif a==b
        order_equal
    else
        order_larger
    end

    cnd = Dict(((i,j) => _order(order[i],order[j]) for i ∈ 1:num_clocks for j ∈ 1:num_clocks)...)
    return ClockRegion{num_clocks,T}(slices,cnd)
end

function ClockRegion{num_clocks,T}(v::Pair{Int,<:TimeSlice{T}}...;order=Dict{Tuple{Int,Int},ORDER_RELATIONSHIP}()) where {num_clocks,T<:Real}
    r = Vector{TimeSlice{T}}(fill(TimeUnconstrained{T}(), num_clocks))
    
    for (clk,val) ∈ v
        r[clk] = val
    end

    ClockRegion{num_clocks,T}(r, order)
end


"""Does the constraint required by region c2 always hold within region c1?"""
function implies(c1::ClockRegion{N},c2::ClockRegion{N}) where N
    for i ∈ 1:N
        v1 = c1.slices[i]
        v2 = c2.slices[i]
        
        if ~(v1⊆v2)
            return false
        end
    end

    for (key,val2) ∈ pairs(c2.order_constraints)
        if ~haskey(c1.order_constraints, key)
            return false
        end

        val1 = c1.order_constraints[key]    
    
        if ~(((val1==order_less) && (val2==order_less || val2==order_less_equal)) ||
            ((val1==order_less_equal) && (val2==order_less_equal)) ||
            ((val1==order_equal) && (val2==order_equal)) ||
            ((val1==order_larger) && (val2==order_larger || val2==order_larger_equal)) ||
            ((val1==order_larger_equal) && (val2==order_larger_equal)))
            return false
        end
    end
    return true
end

macro Clk_str(str)
    T=Int
    r_point = r"^([0-9.]*)$"
    r_interval = r"^([[(])([0-9.∞]*),([0-9.∞]*)([])])$"

    return if str == "true"
    else
        m = match(r_point,str)
        if isa(m,Nothing)
            m = match(r_interval, str)
            if isa(m,Nothing)
                error("Can't parse expression '$(str)'.")
            else
                (lb,v1,v2,rb) = m.captures
                if v1 == v2 == "∞"
                    error("Can't parse expression '$(str)'.")
                elseif v1 == "∞"
                    :(TimeLeftHalfLine{$T,$(rb==")")}($(parse(T,v2))))
                elseif v2 == "∞"
                    :(TimeRightHalfLine{$T,$(lb=="(")}($(parse(T,v1))))
                else
                    :(TimeInterval{$T,$(lb=="("),$(rb==")")}($(parse(T,v1)),$(parse(T,v2))))
                end
            end
        else
            (v,) = m.captures
            :(TimePoint{$T}($(parse(T,v))))
        end
    end
end


function collect_constants!(c::ClockRegion{N,T}, counter::Vector) where {N,T}
    for (i,s) ∈ enumerate(c.slices)
        if isa(s,TimePoint)
            push!(counter[i], s.val)
        elseif isa(s,TimeInterval)
            push!(counter[i], s.val1, s.val2)
        elseif isa(s,TimeLeftHalfLine)
            push!(counter[i], s.val2)
        elseif isa(s,TimeRightHalfLine)
            push!(counter[i], s.val1)
        end
    end
    counter
end