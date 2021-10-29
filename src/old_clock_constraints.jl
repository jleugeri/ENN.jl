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

Base.hash(t::TimePoint{T}) where T= hash((hash(TimePoint),hash(T),hash(t.val)))
Base.hash(t::TimeUnconstrained{T}) where T = hash((hash(TimeUnconstrained),hash(T)))
Base.hash(t::TimeInterval{T,OL,OR}) where {T,OL,OR} = hash((hash(TimeInterval),hash(T),hash(OL),hash(OR),hash(t.val1),hash(t.val2)))
Base.hash(t::TimeRightHalfLine{T,OL}) where {T,OL} = hash((hash(TimeRightHalfLine),hash(T),hash(OL),hash(t.val1)))
Base.hash(t::TimeLeftHalfLine{T,OR}) where {T,OR} = hash((hash(TimeLeftHalfLine),hash(T),hash(OR),hash(t.val2)))
Base.:(==)(::TimeSlice,::TimeSlice) = false
Base.:(==)(r1::TimePoint,r2::TimePoint) = r1.val==r2.val
Base.:(==)(r1::TimeInterval{T,OL,OR},r2::TimeInterval{T,OL,OR}) where {T,OL,OR} = (r1.val1==r2.val1) && (r1.val2==r2.val2)
Base.:(==)(r1::TimeRightHalfLine{T,OL},r2::TimeRightHalfLine{T,OL}) where {T,OL} = r1.val1==r2.val1
Base.:(==)(r1::TimeLeftHalfLine{T,OR},r2::TimeLeftHalfLine{T,OR}) where {T,OR} = r1.val2==r2.val2

Base.:⊆(r1::TimeUnconstrained,r2::TimeUnconstrained) where T = true
Base.:⊆(r1::TimeUnconstrained,r2::T) where T = false
Base.:⊆(r1::T,r2::TimeUnconstrained) where T = true

Base.:⊆(r1::TimePoint, r2::TimePoint) = r1==r2
Base.:⊆(r1::TimePoint, r2::TimeInterval{T,OL,OR}) where {T,OL,OR} = comp(Val(~OL),r2.val1,r1.val) && comp(Val(~OR),r1.val,r2.val2)
Base.:⊆(r1::TimePoint, r2::TimeLeftHalfLine{T,OR}) where {T,OR} = comp(Val(~OR),r1.val,r2.val2)
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

abstract type OrderRelation{N} end

struct PartialOrder{N} <: OrderRelation{N}
    rel::Matrix{Symbol}
    function PartialOrder(n::Int, p::Pair{Tuple{Int,Int},Symbol}...)
        less = zeros(Bool, (n, n))
        larger = zeros(Bool, (n, n))
        equal = zeros(Bool, (n, n))

        for ((s1,s2),o) ∈ p
            if o == :(<=)
                less[s2,s1] = true
                larger[s1,s2] = true
                equal[s2,s1] = true
                equal[s1,s2] = true
            elseif o == :(<)
                less[s2,s1] = true
                larger[s1,s2] = true
            elseif o == :(==)
                equal[s2,s1] = true
                equal[s1,s2] = true
            elseif o == :(>)
                less[s1,s2] = true
                larger[s2,s1] = true
            elseif o == :(>=)
                less[s1,s2] = true
                larger[s2,s1] = true
                equal[s1,s2] = true
                equal[s2,s1] = true
            end
        end

        # a==a
        for i ∈ 1:n
            equal[i,i] = true
        end

        less′ = copy(less)
        larger′ = copy(larger)
        equal′ = copy(equal)
        # ensure transitivity
        for i ∈ 1:n
            less += less*less′
            larger += larger*larger′
            equal += equal*equal′
        end
        less = less .> 0
        larger = larger .> 0
        equal = equal .> 0

        rel = Matrix{Symbol}(undef, (n, n))

        #@assert (all(less .! larger) && all(less .! less') && all(larger .! larger') && all(equal .== equal')) "Not a valid (partial) order relation!"

        for i ∈ 1:n
            for j ∈ 1:n
                rel[i,j] = if less[i,j] && larger[i,j]
                    :?
                elseif less[i,j] && equal[i,j]
                    :(<=)
                elseif less[i,j]
                    :(<)
                elseif larger[i,j] && equal[i,j]
                    :(>=)
                elseif larger[i,j]
                    :(>)
                elseif equal[i,j]
                    :(==)
                else
                    :?
                end
            end
        end
        return new{n}(rel)
    end
end

(o::PartialOrder)(a::Int, b::Int) = o.rel[b,a]

struct TotalOrder{N} <: OrderRelation{N}
    order::Vector{Int}
end
TotalOrder(v::Vector{Int}) = TotalOrder{length(v)}(v)

ismaximum(o::TotalOrder) = o.order .== maximum(o.order)

(o::TotalOrder)(a::Int, b::Int) = if o.order[a]<o.order[b]
    :<
elseif o.order[a]==o.order[b]
    :(==)
elseif o.order[a]>o.order[b]
    :>
else
    :?
end


Base.hash(o::TotalOrder{N}) where N = hash((hash(TotalOrder),hash(N),hash(o.order)))
Base.:(==)(o1::TotalOrder,o2::TotalOrder) = all(o1.order .== o2.order)
#Base.:(==)(o1::OrderRelation,o2::OrderRelation) = (o1 ⊆ o2) && (o2 ⊆ o1)

function Base.:⊆(o1::OrderRelation{N},o2::OrderRelation{N}) where N
    # check if each ordering of pairs in o1 is *equally or more specific* than (and compatible with) the corresponding ordering in o2
    for i ∈ 1:N
        for j ∈ 1:N
            val1 = o1(i,j)
            val2 = o2(i,j)

            if val2 == :?
                # o2 doesn't care about this relation
                continue
            elseif val1 == :? && val2 != :?
                # o1 doesn't care about this relation -> can't be more specific than anything
                return false
            elseif ~(((val1==:<) && (val2==:< || val2==:(<=))) ||
                    ((val1==:(<=)) && (val2==:(<=))) ||
                    ((val1==:(==)) && (val2==:(==) || val2==:(<=) || val2==:(>=))) ||
                    ((val1==:>) && (val2==:> || val2==:(>=))) ||
                    ((val1==:(>=)) && (val2==:(>=))))
                # implication violated by this clock
                return false
            end
        end
    end

    # if we got this far, there're no contradictions in the two order relations -> o1 ⊆ o2
    return true
end    

struct ClockRegion{N,T}
    slices::Vector{TimeSlice{T}}
    order::OrderRelation
end
function ClockRegion{num_clocks,T}(v::Pair{Int,<:TimeSlice{T}}...;order=PartialOrder(num_clocks)) where {num_clocks,T<:Real}
    r = Vector{TimeSlice{T}}(fill(TimeUnconstrained{T}(), num_clocks))
    
    for (clk,val) ∈ v
        r[clk] = val
    end

    ClockRegion{num_clocks,T}(r, order)
end

Base.:(==)(c1::ClockRegion{N},c2::ClockRegion{N}) where N = all(c1.slices .== c2.slices) && (c1.order == c2.order)
Base.hash(r::ClockRegion{N,T}) where {N,T} = hash((hash(ClockRegion),hash(N),hash(T),hash(r.slices),hash(r.order)))

"""Check if the constraints of clock region 1 are compatible with, yet more specific than the constraints of clock region 2"""
function Base.:⊆(c1::ClockRegion{N},c2::ClockRegion{N}) where N
    # check if each clock from clock region 1 is *more constrained* than the corresponding one from clock region 2
    for i ∈ 1:N
        v1 = c1.slices[i]
        v2 = c2.slices[i]
        
        if ~(v1⊆v2)
            return false
        end
    end

    # if we got this far, it all depends on the order relationship
    # check if the order from clock region 1 is *more constrained and compatible with* the corresponding one from clock region 2
    return c1.order ⊆ c2.order
end

macro Clk_str(str)
    T=Int
    r_point = r"^([0-9.]*)$"
    r_interval = r"^([[(])([0-9.∞]*),([0-9.∞]*)([])])$"

    return if str == ""
        :(TimeUnconstrained{$T}())
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