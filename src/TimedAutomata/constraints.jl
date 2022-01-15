export TAConstraint, TADiagonalConstraint, get_diagonal_constraints, set_clock_list!

struct TADiagonalConstraint{T}
    clock1::Union{Nothing, Symbol}
    clock2::Union{Nothing, Symbol}
    bound::TABound{T}
end
Base.hash(c::TADiagonalConstraint) = hash((hash(c.clock1), hash(c.clock2), hash(c.bound)))
Base.:(==)(c1::TADiagonalConstraint{T},c2::TADiagonalConstraint{T}) where T = c1.clock1==c2.clock1 && c1.clock2==c2.clock2 && c1.bound==c2.bound
Base.isequal(c1::TADiagonalConstraint,c2::TADiagonalConstraint) = c1==c2
function Base.:!(c::TADiagonalConstraint{T}) where T
    bound = TABound(!c.bound.strict,-c.bound.value)
    TADiagonalConstraint{T}(c.clock2, c.clock1, bound)
end

Base.iterate(x::TADiagonalConstraint) = (x, nothing)
Base.iterate(::TADiagonalConstraint, ::Any) = nothing
Base.length(x::TADiagonalConstraint) = 1
Base.Broadcast.broadcastable(c::TADiagonalConstraint) = Ref(c)
Base.adjoint(c::TADiagonalConstraint) = c

function Base.show(io::IO, M::MIME"text/plain", c::TADiagonalConstraint)   
    res = if isnothing(c.clock1) && isnothing(c.clock2) && c.bound.value ≥ 0
        "⊤"
    elseif isnothing(c.clock1) && isnothing(c.clock2) && c.bound.value < 0
        "⊥"
    elseif isnothing(c.clock1)
        v = if isinf(c.bound) && c.bound.value > 0 
            "-∞"
        elseif isinf(c.bound) && c.bound.value < 0
            "∞"
        else
            "$(-c.bound.value)"
        end
        "$(c.clock2) $(c.bound.strict ? :> : :≥) $(v)"
    elseif isnothing(c.clock2)
        v = if isinf(c.bound) && c.bound.value > 0 
            "∞"
        elseif isinf(c.bound) && c.bound.value < 0
            "-∞"
        else
            "$(c.bound.value)"
        end
        "$(c.clock1) $(c.bound.strict ? :< : :≤) $(v)"
    else
        v = if isinf(c.bound) && c.bound.value > 0 
            "∞"
        elseif isinf(c.bound) && c.bound.value < 0
            "-∞"
        else
            "$(c.bound.value)"
        end
        "$(c.clock1)-$(c.clock2) $(c.bound.strict ? :< : :≤) $(v)"
    end
    print(io,res)
end

mutable struct TAConstraint{T}
    clocks::Vector{Symbol}
    D::Matrix{TABound{T}}
end

function TAConstraint(terms::Vector{TADiagonalConstraint{T}}, clocks=Symbol[]) where T
    if isempty(clocks)
        for term ∈ terms
            if !isnothing(term.clock1) && !(term.clock1 in clocks)
                push!(clocks, term.clock1)
            end
            if !isnothing(term.clock2) && !(term.clock2 in clocks)
                push!(clocks, term.clock2)
            end
        end
    end
    D = make_dbm(terms, clocks)
    TAConstraint(clocks, D)
end

Base.copy(ta::TAConstraint) = TAConstraint(ta.clocks,ta.D)
Base.deepcopy(ta::TAConstraint) = TAConstraint(deepcopy(ta.clocks), deepcopy(ta.D))

function make_dbm(terms::Vector{TADiagonalConstraint{T}}, clocks) where T
    num_clocks = length(clocks)+1
    D = fill(typemax(TABound{T}), num_clocks, num_clocks)
    strict = zeros(Bool, num_clocks, num_clocks)
    lookup = Dict(zip(clocks,keys(clocks)))

    # Diagonal entries: all clocks are identical to themselves
    D[diagind(D)] .= zero(TABound{T})

    # first row: all clocks are positive
    D[1,:] .= zero(TABound{T})

    # Add specific constraints
    for term in terms
        key1 = isnothing(term.clock1) ? 1 : lookup[term.clock1]+1
        key2 = isnothing(term.clock2) ? 1 : lookup[term.clock2]+1

        # Refine constraints
        D[key1,key2] = min(term.bound,D[key1,key2])
    end
    return D
end


function Base.intersect!(c::TAConstraint{T}, g::TADiagonalConstraint) where T
    idx1 = isnothing(g.clock1) ? 1 : findfirst(==(g.clock1), c.clocks)+1
    idx2 = isnothing(g.clock2) ? 1 : findfirst(==(g.clock2), c.clocks)+1
    
    if c.D[idx2,idx1] < -g.bound || (isinf(g.bound) && g.bound.value < 0)
        c.D[1,1]=TABound(true,-one(T))
    elseif g.bound < c.D[idx1,idx2]
        c.D[idx1,idx2] = g.bound
        c.D .= min.(c.D, c.D[:,idx1] .+ c.D[idx1,:]', c.D[:,idx2] .+ c.D[idx2,:]')
    end
    c
end
Base.intersect(c::TAConstraint, g::TADiagonalConstraint) = intersect!(deepcopy(c), g)


function Base.intersect!(c1::TAConstraint, others::TAConstraint...)
    c1.D .= min.(c1.D, getfield.(others, :D)...)
    close!(c1)
    c1
end
Base.intersect(c1::TAConstraint, others::TAConstraint...) = intersect!(deepcopy(c1),others...)


Base.iterate(x::TAConstraint) = (x, nothing)
Base.iterate(::TAConstraint, ::Any) = nothing
Base.length(x::TAConstraint) = 1
Base.Broadcast.broadcastable(c::TAConstraint) = Ref(c)
Base.adjoint(c::TAConstraint) = c


function Base.issubset(c1::TAConstraint, c2::TAConstraint) 
    @assert size(c1.D) == size(c2.D) "Must be of same size, got $(size(c1.D)) and $(size(c2.D))"
    for idx ∈ eachindex(c1.D)
        if c1.D[idx] > c2.D[idx]
            return false
        end
    end
    return true
end

Base.hash(c::TAConstraint) = hash((hash(c.clocks), hash(c.D)))
Base.:(==)(c1::TAConstraint, c2::TAConstraint) = (c1.clocks == c2.clocks) && (c1.D==c2.D)

function Base.isempty(c::TAConstraint) 
    D  = c.D
    D′ = D'
    for idx ∈ eachindex(c.D)
        c1 = D[idx]
        c2 = D′[idx]
        # Situation: -m2 </≤ x-y </≤ m1 
        # if either of the inequalities is strict, we have a problem if m1 ≤ -m2
        # if neither of the inequalities is strict, we only have a problem if m1 < -m2
        if (
                (c1.strict || c2.strict) && (
                    c1.value ≤ -c2.value || 
                    (isinf(c1.value) && c1.value < 0) || 
                    (isinf(c2.value) && c2.value < 0)
                )
            ) || (c1.value < -c2.value)
            return true
        end
    end
    return false
end

function set_clock_list!(ta::TAConstraint{T}, clocks::Vector{Symbol}; mapping=identity) where T
    clocks = copy(clocks)
    D = make_dbm(TADiagonalConstraint{T}[], clocks)

    @assert mapping.(ta.clocks) ⊆ clocks "Not all clocks of the constraint set ($(mapping.(ta.clocks))) are in the given list of clocks ($(clocks))."

    lookup = Dict(zip(clocks,keys(clocks)))
    for (i,clock1) ∈ enumerate([nothing; ta.clocks])
        idx1 = isnothing(clock1) ? 1 : lookup[mapping(clock1)]+1
        for (j,clock2) ∈ enumerate([nothing; ta.clocks])
            idx2 = isnothing(clock2) ? 1 : lookup[mapping(clock2)]+1
            D[idx1,idx2] = ta.D[i,j]
        end
    end

    ta.clocks = clocks
    ta.D = D
    return ta
end

function get_diagonal_constraints(b::TAConstraint{T}) where T
    if length(b.clocks) > 0
        TADiagonalConstraint{T}.([nothing; b.clocks], reshape([nothing; b.clocks],(1,:)), b.D)[2:end]
    else
        TADiagonalConstraint{T}[]
    end
end

function Base.show(io::IO, ::MIME"text/plain", b::TAConstraint) 
    if isempty(b)
        print(io, "∅")
        return
    end
    constraints = String[]
    for j ∈ 1:length(b.clocks)+1
        for i ∈ j+1:length(b.clocks)+1
            c1 = b.D[i,j]
            c2 = b.D[j,i]
            clock = j==1 ? "$(b.clocks[i-1])" : "$(b.clocks[i-1])-$(b.clocks[j-1])"
            v1 = isinf(c1) ? (c1.value>0 ? "∞" : "-∞") : "$(c1.value)"
            v2 = isinf(c2) ? (c2.value>0 ? "-∞" : "∞") : "$(-c2.value)"
            if -c2 < c1
                if isinf(c1) && isinf(c2)
                    ""
                else
                    lhs = isinf(c2) || (j==1 && iszero(c2.value)) ? "" : "$(v2) $(c2.strict ? :< : :≤) "
                    rhs = isinf(c1) ? "" : " $(c1.strict ?  :< : :≤) $(v1)"
                    if !isempty(lhs) || !isempty(rhs)
                        push!(constraints, lhs*"$(clock)"*rhs)
                    end
                end
            elseif c1 < -c2
                push!(constraints, "⊥")
            else
                push!(constraints, "$(clock) == $(v1)")
            end
        end
    end
    print(io, join(constraints, "\n"))
end
