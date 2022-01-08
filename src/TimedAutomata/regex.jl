
abstract type TARegex end

struct TARegexAtom{T} <: TARegex
    symbol::T
end

struct TARegexRepetition <: TARegex
    content::TARegex
end

struct TARegexUnion <: TARegex
    terms::Set{TARegex}
end
TARegexUnion(r::TARegexUnion) = r
TARegexUnion(r::TARegex) = TARegexUnion(Set([r]))
Base.:+(r1::TARegexUnion, r2::TARegexUnion) = TARegexUnion(r1.terms ∪ r2.terms)
Base.:+(r1::TARegex, r2::TARegex) = TARegexUnion(r1)+TARegexUnion(r2)
#Base.zero(::Type{TARegexUnion}) = TARegexUnion(Set{TARegex}())

struct TARegexConcatenation <: TARegex
    terms::Vector{TARegex}
end
TARegexConcatenation(r::TARegexConcatenation)=r
TARegexConcatenation(r::TARegex) = TARegexConcatenation([r])
Base.zero(::Type{T}) where T <:TARegex = TARegexUnion(Set{TARegex}())
Base.one(::Type{T}) where T <:TARegex = TARegexConcatenation(TARegex[])
Base.zero(::T) where T <:TARegex = zero(T)
Base.one(::T) where T <:TARegex = one(T)
Base.iszero(x::TARegexUnion) = isempty(x.terms)
Base.iszero(x::TARegexConcatenation) = length(x.terms)==1 && iszero(first(x.terms))
Base.iszero(x::TARegexRepetition) = iszero(x.content)
Base.iszero(x::TARegex) = false

Base.:*(r1::TARegexConcatenation, r2::TARegexConcatenation) = TARegexConcatenation([r1.terms ; r2.terms])
Base.:*(r1::TARegexUnion, r2::TARegexUnion) = TARegexUnion(Set([t1*t2 for t1 in r1.terms for t2 in r2.terms]))
Base.:*(r1::TARegex, r2::TARegexUnion) = TARegexUnion(Set(r1 .* r2.terms))
Base.:*(r1::TARegexUnion, r2::TARegex) = TARegexUnion(Set(r1.terms .* r2))
Base.:*(r1::TARegex, r2::TARegex) = TARegexConcatenation(r1)*TARegexConcatenation(r2)

Base.hash(x::TARegexConcatenation) = hash(x.terms)
Base.hash(x::TARegexUnion) = hash(x.terms)
Base.hash(x::TARegexRepetition) = hash(x.content)
Base.hash(x::TARegexAtom) = hash(x.symbol)
Base.:(==)(x::TARegex, y::TARegex) = hash(x)==hash(y)

Base.broadcastable(x::TARegex) = Ref(x)
Base.show(io::IO, ::MIME"text/plain", x::TARegexAtom) = print(io, x.symbol)
Base.show(io::IO, ::MIME"text/plain", x::TARegexRepetition) = (s=repr("text/plain", x.content); print(io, "($(s))*"))
Base.show(io::IO, ::MIME"text/plain", x::TARegexUnion) = print(io, join(repr.(Ref("text/plain"),x.terms), "|"))
Base.show(io::IO, ::MIME"text/plain", x::TARegexConcatenation) = print(io, join(repr.(Ref("text/plain"),x.terms), "→"))

Base.promote_rule(::Type{T1},::Type{T2}) where {T1<:Union{TARegexAtom,TARegexConcatenation,TARegexRepetition,TARegexUnion},T2<:Union{TARegexAtom,TARegexConcatenation,TARegexRepetition,TARegexUnion}} = T1==T2 ? T1 : TARegex

"""
Brzozowski's method
"""
function language(ts::TTS{State,T}; forbid_messages=(_)->false, ignore_messages=(m)->m.direction==MSG_SILENT, final_states=ts.states, initial_state=ts.states[1]) where {State,T}
    X = zeros(TARegex,length(ts.states))
    M = zeros(TARegex, length(ts.states), length(ts.states))
    
    
    for state in final_states
        id = findfirst(==(state), ts.states)
        @assert !isnothing(id) "Final state $(state) not in list of states!"
        X[id] = one(TARegex)
    end
    
    lookup(s) = findfirst(==(s),ts.states)
    
    for (source,(message,target)) in ts.transitions
        M[lookup(source), lookup(target)] += if forbid_messages(message)
                zero(TARegex)
            elseif ignore_messages(message)
                one(TARegex)
            else
                TARegexAtom(message) 
            end
    end
    
    initial_id = findfirst(==(initial_state), ts.states)
    @assert !isnothing(initial_id) "Initial state $(initial_state) not in state list."

    # make sure initial state is at position 1
    permutation = circshift(collect(1:length(ts.states)),1-initial_id)
    X = X[permutation]
    M = M[permutation,permutation]

    for n in reverse(eachindex(ts.states))
        if !iszero(M[n,n])
            X[n] = TARegexRepetition(M[n,n])*X[n]
            M[n,1:n-1] .= TARegexRepetition(M[n,n]) .* M[n,1:n-1]
        end
        X[1:n-1] .+= M[1:n-1,n] .* X[n]
        M[1:n-1,1:n-1] .+= M[1:n-1,n].*reshape(M[n,1:n-1],(1,:))
    end
    return X[1]
end
