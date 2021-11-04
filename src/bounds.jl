export TABound

mutable struct TABound{T}
    strict::Bool
    value::T
    infinite::Bool
    function TABound(s,v::T) where T
        b=new{T}(s,v,false)
        b.value=v
        b
    end
end
function Base.setproperty!(b::TABound{T}, name::Symbol, value) where T
    if name == :strict
        setfield!(b, :strict, value)
    elseif name == :value
        if isinf(value) || value==typemax(value) || value==typemin(value)
            if value > 0
                setfield!(b, :value, typemax(T))
            else
                setfield!(b, :value, typemin(T))
            end
            setfield!(b, :infinite, true)
        else
            setfield!(b, :value, value)
            setfield!(b, :infinite, false)
        end
    else
        throw(ArgumentError("Can't set property '$(name)' of '$(b)'!"))
    end
end

Base.isinf(b::TABound) = b.infinite
Base.iszero(b::TABound) = iszero(b.value)
Base.hash(b::TABound) = hash((hash(b.strict), hash(b.value)))
Base.:(==)(b1::TABound, b2::TABound) = ((b1.value == b2.value) && (b1.strict == b2.strict)) || (isinf(b1) && isinf(b2))
Base.isequal(b1::TABound, b2::TABound) = b1==b2
Base.isless(b1::TABound{T}, b2::TABound) where {T} = (!isinf(b1) && isinf(b2) && b2.value > zero(T)) || (b1.value < b2.value) || ((b1.value == b2.value) && b1.strict && !b2.strict )
(Base.zero(::Type{TABound{T}}) where T) = TABound(false, zero(T))
(Base.typemax(::Type{TABound{T}}) where T) = TABound(false, typemax(T))
(Base.typemin(::Type{TABound{T}}) where T) = TABound(false, typemin(T))
function Base.:+(b1::TABound{T}, b2::TABound) where T
    if isinf(b1) && isinf(b2) && sign(b1.value)!=sign(b2.value)
        throw(ArgumentError("Sum of bounds undefined! (∞ - ∞)"))
    elseif isinf(b1) && b1.value > 0 || isinf(b2) && b2.value > 0
        typemax(TABound{T})
    elseif isinf(b1) && b1.value < 0 || isinf(b2) && b2.value < 0
        typemin(TABound{T})
    else
        TABound(b1.strict || b2.strict, b1.value + b2.value)
    end
end
Base.:-(c::TABound) = TABound(c.strict, -c.value)
Base.:-(c1::TABound, c2::TABound) = c1 + (-c2)

Base.copy(b::TABound) = TABound(b.s,b.v)
Base.deepcopy(b::TABound) = TABound(deepcopy(b.s),deepcopy(b.v))
Base.iterate(x::TABound) = (x, nothing)
Base.iterate(::TABound, ::Any) = nothing
Base.length(x::TABound) = 1
Base.Broadcast.broadcastable(b::TABound) = Ref(b)
Base.adjoint(c::TABound) = c

function Base.show(io::IO, M::MIME"text/plain", b::TABound) 
    s = if isinf(b) && b.value > 0
        "∞"
    elseif isinf(b) && b.value ≤ 0
        "-∞"
    elseif b.strict
        "< "*repr(b.value)
    else
        "≤ "*repr(b.value)
    end
    print(io, s)
end
