export close!, close, up!, up, reset!, reset, zone_graph, get_clock_ceilings, get_diagonal_constraints, get_diagonal_guards, normalize_zone, normalize_zone_k


function close!(c::TAConstraint)
    n = size(c.D,1)
    for k ∈ 1:n
        for i ∈ 1:n
            for j ∈ 1:n
                c.D[i,j] = min(c.D[i,j], c.D[i,k]+c.D[k,j])
            end
        end
    end
    c
end
close(c::TAConstraint) = close!(copy(c))

function up!(c::TAConstraint)
    c.D[2:end,1] .= typemax(eltype(c.D))
    close!(c)
    c
end
up(c::TAConstraint) = up!(copy(c))

function reset!(c::TAConstraint, clocks::Set{Symbol})
    for clock ∈ clocks
        idx = findfirst(==(clock), c.clocks)+1
        c.D[idx,:] .= c.D[1,:]
        c.D[:,idx] .= c.D[:,1]
    end
    close!(c)
    c
end

reset(c::TAConstraint, clocks::Set{Symbol}) = reset!(copy(c), clocks)

function zone_graph(ta::TA{State, T}; k_fun=kk->fill(maximum(kk),size(kk))) where {State, T}
    lD_initial = (
        ta.initial_state, 
        ta.invariants[ta.initial_state]
    )
    states = Tuple{State,TAConstraint{T}}[lD_initial]
    transitions = Pair{Tuple{State,TAConstraint{T}},Tuple{TAMessage,Tuple{State,TAConstraint{T}}}}[]
    max_runs = 1e3

    Gd = get_diagonal_guards(ta)
    k  = k_fun(get_clock_ceilings(ta))

    function add_transition!((l,D),msg,(l′, D′))
        idx = findfirst(==((l′, D′)), states)
        if isnothing(idx)
            # if we just found a new state + clock region, save for later!
            push!(states, (l′, D′))
        end
        push!(transitions, (l,D) => (msg , (l′, D′)))
    end

    last = 0
    while last < length(states) && max_runs > 0
        max_runs -=1
        last += 1
        (l,D) = states[last]

        # Delay transition:
        # (l,D) ⇝ (l,D′), D′ = D↑ ∩ I(l)
        l′ = l
        D′ = up(D) ∩ ta.invariants[l]

        # we don't need loops within the same zone
        if D′!=D 
            # split & normalize the zone
            for D′′ in normalize_zone(D′, Gd, k)
                # if this state is consistent, keep it
                if !isempty(D′′)
                    add_transition!((l,D),TAMessage(),(l′,D′′))
                end
            end
        end
        
        # Action transitions:
        # (l,D) ⇝ (l′, D′), D′ = r(D∧g) ∧ I(l′)  if  l →g,a,r→ l′
        for arc ∈ ta.arcs_by_state[l]
            # follow the arc
            l′ = arc.target
            
            # The arc constrains the zone
            # if the arc is compatible with the current time zone, keep going
            D′ = D ∩ arc.guard
            if isempty(D′)
                continue
            end

            # the transition resets some clocks and must satisfy the new state's guards
            D′ = TimedAutomata.reset(D′, arc.resets) ∩ ta.invariants[l′]
            if isempty(D′)
                continue
            end
            
            # split & normalize the zone
            for D′′ in normalize_zone(D′, Gd, k)
                if !isempty(D′′)
                    add_transition!((l,D),arc.message,(l′,D′′))
                end
            end

        end
    end


    return (states=states, transitions=transitions)
end


function get_diagonal_guards(ta::TA{State,T}) where {State,T}
    Gd = Vector{TADiagonalConstraint{T}}()
    if !isempty(ta.clocks)
        for arc in ta.arcs
            tmp = get_diagonal_constraints(arc.guard)[2:end]
            filter!(tmp) do arc
                !isnothing(arc.clock1) && !isnothing(arc.clock2) && arc.clock1 != arc.clock2
            end

            append!(Gd, tmp)
        end
    end
    unique!(Gd)
    Gd
end

function split(D::TAConstraint{T},Gd) where {T}
    Q = Set{TAConstraint{T}}([D])
    Q′= Set{TAConstraint{T}}()

    for g in Gd
        for D′ in Q
            D1 = D′ ∩ g
            D2 = D′ ∩ !g
            if !isempty(D1) && !isempty(D2)
                # if g cuts the region in two, keep both parts
                push!(Q′,D1)
                push!(Q′,D2)
            else
                # otherwise keep the region
                push!(Q′,copy(D′))
            end
        end
        # swap buffers and empty Q′
        (Q,Q′) = (Q′,Q)
        empty!(Q′)
    end
    Q
end

function normalize_zone_k!(c::TAConstraint{T},k::Vector{T}) where T
    for i in 1:length(c.clocks)+1
        for j in 1:length(c.clocks)+1
            if !isinf(c.D[i,j]) && TABound(false,k[i]) < c.D[i,j]
                c.D[i,j] = typemax(TABound{T})
            elseif !isinf(c.D[i,j]) && c.D[i,j] < TABound(true,-k[j])
                c.D[i,j] = TABound(true,-k[j])
            end
        end
    end
    close!(c)
    c
end

function normalize_zone_k(c::TAConstraint,k)
    _c=copy(c)
    normalize_zone_k!(_c,k)
    _c
end

function normalize_zone(D::TAConstraint{T},Gd,k) where T
    Q = Vector{TAConstraint{T}}()
    for D′ in split(D,Gd)
        G_unsat = Set{TADiagonalConstraint{T}}()
        for g in Gd
            ng = !g
            if isempty(D′ ∩ g)
                push!(G_unsat, ng)
            end
            if isempty(D′ ∩ ng)
                push!(G_unsat, g)
            end
        end

        normalize_zone_k!(D′,k)
        for g in G_unsat
            intersect!(D′,g)
        end
        push!(Q,D′)
    end
    Q
end

function get_clock_ceilings(c::TAConstraint{T}) where T
    max.(
        maximum(x->isinf(x) ? -1 : x.value, c.D, dims=2)[:],
        maximum(x->isinf(x) ? -1 : -x.value, c.D, dims=1)[:]
    )
end

function get_clock_ceilings(ta::TA{State,T}) where {State, T}
    k = zeros(T, length(ta.clocks)+1)
    for arc in ta.arcs
        k .= max.(k, get_clock_ceilings(arc.guard))
    end

    for inv in values(ta.invariants)
        k .= max.(k, get_clock_ceilings(inv))
    end
    k
end

