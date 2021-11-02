export @arc, @constraint 

macro arc(s1, s2, constraints, messages, resets, T=:(Rational{Int}))
    cnst = parse_constraint_set(constraints, T)
    esc(:(TAArc{typeof($s1),$T}($s1, $s2, $cnst, $messages, Set($resets))))
end


function parse_constraint_set(expr, T)
    acc = Tuple{Union{Symbol,Nothing},Union{Symbol,Nothing},Bool,Any}[]
    
    function make_inequality(clock1, clock2, op, c)
        # Flip greater-than(-or-equal) operators
        if op ∈ (:≥,:>)
            op = (op == :≥) ? :≤ : :<
            clock1,clock2 = clock2,clock1
            c = :(-($c)) 
        end
        strict = op ∈ (:<,:>)
        (clock1, clock2, strict,c)
    end

    function parse_diagonal_constraints(e::Expr)
        error(e) = "Must be of the form: 'true' or 'x-y ∼ c' or 'x ~ c' or 'c ∼ x' or 'c ∼ x-y' where x and y are clocks, c is a constant and ∼ ∈ {<,≤,==,≥,>}. (Got: $(e) instead)."
        @assert e.head == :call && !isempty(e.args) error(e)
        op = e.args[1]
        @assert op ∈ (:<,:≤,:(==),:≥,:>) && length(e.args)==3 error(e) 
        arg1,arg2 = e.args[2:3]    

        flip = Dict(:< => :>,:≤ => :≥, :> => :<,:≥ => :≤)

        clock1,clock2,op,bound = if isa(arg1,Expr) && arg1.head == :call && length(arg1.args)==3 && arg1.args[1] == :-
            (arg1.args[2:3]...,op,arg2)
        elseif isa(arg2,Expr) && arg2.head == :call && length(arg2.args)==3 && arg2.args[1] == :-
            (arg2.args[2:3]...,flip[op],arg1)
        elseif isa(arg1,Symbol)
            (arg1,nothing,op,arg2)
        elseif isa(arg2,Symbol)
            (arg2,nothing,flip[op],arg1)
        else
            throw(AssertionError(error(e)))
        end
    
        return if   op == :(==)
            [make_inequality(clock1,clock2,:≤,bound),make_inequality(clock1,clock2,:≥,bound)]
        else
            [make_inequality(clock1,clock2,op,bound)]
        end
    end
    
    function recurse_constraint_set!(e, acc)
        if isa(e, Bool)
            @assert e "A constant false constraint is not supported!"
            return
        elseif e.head == :(&&)
            for arg in e.args
                recurse_constraint_set!(arg, acc)
            end
        else
            append!(acc, parse_diagonal_constraints(e))
        end
    end

    recurse_constraint_set!(expr, acc)

    # Infer proper type for time-constants
    if isnothing(T)
        T = promote_type(getindex.(getfield.(typeof.(acc),:parameters), 4)...)
    end

    vals = [:(TADiagonalConstraint{$T}($(isa(c1,Symbol) ? QuoteNode(c1) : c1),$(isa(c2,Symbol) ? QuoteNode(c2) : c2), TABound($(s),$(m)))) for (c1,c2,s,m) in acc]
    return :(TAConstraint(TADiagonalConstraint{$T}[$(vals...)]))
end

macro constraint(expr, T=:(Rational{Int}))
    res = parse_constraint_set(expr, T)    
    esc(:($res))
end
