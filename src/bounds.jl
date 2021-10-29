export TABound

mutable struct TABound{T}
    strict::Bool
    value::T
    infinite::Bool
    function TABound(s,v::T) where T
        new{T}(s,v,isinf(v))
    end
end
function Base.setproperty!(b::TABound, name::String, value)
    if name == :strict
        setfield!(b, :strict, value)
    elseif name == :value
        if isinf(value)
            setfield!(b, :value, typemax(T))
            setfield!(b, :infinite, true)
        else
            setfield!(b, :value, value)
            setfield!(b, :infinite, false)
        end
    else
        ArgumentError("Can't set property '$(name)' of '$(b)'!")
    end
end

Base.isinf(b::TABound) = b.infinite
Base.:(==)(b1::TABound, b2::TABound) = ((b1.value == b2.value) && (b1.strict == b2.strict)) || (isinf(b1) && isinf(b2))
Base.isless(b1::TABound, b2::TABound) = (!isinf(b1) && isinf(b2)) || (b1.value < b2.value) || ((b1.value == b2.value) && b1.strict && !b2.strict )
(Base.zero(::Type{TABound{T}}) where T) = TABound(false, zero(T))
(Base.typemax(::Type{TABound{T}}) where T) = TABound(false, typemax(T))
Base.:+(b1::TABound, b2::TABound) = isinf(b1) || isinf(b2) ? Base.typemax(typeof(b1)) : TABound(b1.strict || b2.strict, b1.value + b2.value)
Base.:-(c::TABound) = TABound(c.strict, -c.value)
Base.length(b::TABound) = 1
Base.size(b::TABound) = ()
Base.iterate(b::TABound, step=1) = step==1 ? (b,nothing) : nothing
Base.Broadcast.broadcastable(b::TABound) = Ref(b)
Base.similar(b::TABound, ::Type{T}, dims::Dims) where {T} = SparseArray(T, dims)
function Base.show(io::IO, b::TABound) 
    s = if isinf(b) && b.value > 0
        "∞"
    elseif isinf(b) && b.value ≤ 0
        "-∞"
    elseif b.strict
        "<"*repr(b.value)
    else
        "≤"*repr(b.value)
    end
    print(io, s)
end
Base.adjoint(c::TABound) = c

(Base.BroadcastStyle(::Type{TABound{T}}) where T) = Base.Broadcast.DefaultArrayStyle{0}()
