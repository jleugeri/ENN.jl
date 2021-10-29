export TAConstraint, TADiagonalConstraint

struct TADiagonalConstraint{T}
    clock1::Union{Nothing, Symbol}
    clock2::Union{Nothing, Symbol}
    bound::TABound{T}
end

mutable struct TAConstraint{T}
    clocks::Vector{Symbol}
    D::Matrix{TABound{T}}
end

function TAConstraint(terms::Vector{TADiagonalConstraint{T}}, clocks=Symbol[]) where T
    if isempty(clocks)
        for term ∈ terms
            if !isnothing(term.clock1)
                push!(clocks, term.clock1)
            end
            if !isnothing(term.clock2)
                push!(clocks, term.clock2)
            end
        end
    end
    D = make_dbm(terms, clocks)
    TAConstraint(clocks, D)
end

function Base.copy(ta::TAConstraint{T}) where T
    TAConstraint(copy(ta.clocks), copy(ta.D))
end

function make_dbm(terms::Vector{TADiagonalConstraint{T}}, clocks) where T
    num_clocks = length(clocks)+1
    D = fill(TABound(false,typemax(T)), num_clocks, num_clocks)
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


function Base.intersect!(c1::TAConstraint, c2::TAConstraint)
    c1.D .= min.(c1.D, c2.D)
    close!(c1)
    c1
end

Base.intersect(c1::TAConstraint, c2::TAConstraint) = intersect!(copy(c1),c2)

function Base.issubset(c1::TAConstraint, c2::TAConstraint) 
    @assert size(c1.D) == size(c2.D) "Must be of same size, got $(size(c1.D)) and $(size(c2.D))"
    for idx ∈ eachindex(c1.D)
        if c1.D[idx] > c2.D[idx]
            return false
        end
    end
    return true
end

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
                    (isinf(c2.value) && c2.value > 0)
                )
            ) || (c1.value < -c2.value)
            return true
        end
    end
    return false
end

function set_clock_list!(ta::TAConstraint{T}, clocks::Vector{Symbol}) where T
    D = make_dbm(TADiagonalConstraint{T}[], clocks)

    @assert ta.clocks ⊆ clocks "Not all clocks of the constraint set ($(ta.clocks)') are in the given list of clocks ($(clocks))."

    lookup = Dict(zip(clocks,keys(clocks)))
    for (i,clock1) ∈ enumerate([nothing; ta.clocks])
        idx1 = isnothing(clock1) ? 1 : lookup[clock1]+1
        for (j,clock2) ∈ enumerate([nothing; ta.clocks])
            idx2 = isnothing(clock2) ? 1 : lookup[clock2]+1
            D[idx1,idx2] = ta.D[i,j]
        end
    end

    ta.clocks = clocks
    ta.D = D
    return ta
end

function Base.show(io::IO, ::MIME"text/plain", b::TAConstraint) 
    constraints = String[]
    for j ∈ 1:length(b.clocks)+1
        for i ∈ j+1:length(b.clocks)+1
            c1 = b.D[i,j]
            c2 = b.D[j,i]
            
            
            clock = j==1 ? "$(b.clocks[i-1])" : "$(b.clocks[i-1])-$(b.clocks[j-1])"
            if -c2 < c1
                push!(constraints, "$(-c2.value) $(c2.strict ? :< : :≤) $(clock) $(c1.strict ?  :< : :≤) $(c1.value)")
            elseif c1 < -c2
                push!(constraints, "⊥")
            else
                push!(constraints, "$(clock) == $(c1.value)")
            end
        end
    end
    print(io, join(constraints, " ∧ "))
end
