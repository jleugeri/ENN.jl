using Base.Iterators
export close!, close, up!, up, reset!, reset, zone_graph_fwd, get_clock_ceilings, get_diagonal_constraints, get_diagonal_guards, normalize_zone, normalize_zone_k, reachable_states, unreachable_states, prune, prune!


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
close(c::TAConstraint) = close!(deepcopy(c))

function up!(c::TAConstraint)
    c.D[2:end,1] .= typemax(eltype(c.D))
    close!(c)
    c
end
up(c::TAConstraint) = up!(deepcopy(c))

function reset!(c::TAConstraint, clocks::Set{Symbol})
    for clock ∈ clocks
        idx = findfirst(==(clock), c.clocks)+1
        c.D[idx,:] .= c.D[1,:]
        c.D[:,idx] .= c.D[:,1]
    end
    close!(c)
    c
end

reset(c::TAConstraint, clocks::Set{Symbol}) = reset!(deepcopy(c), clocks)


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
                push!(Q′,deepcopy(D′))
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
    _c=deepcopy(c)
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


function zone_graph_fwd(ta::TA{State, T}; k_fun=kk->fill(maximum(kk),size(kk)), ignore_implicit_delays=true, max_iter = 1e4) where {State, T}
    lD_initial = (
        ta.initial_state, 
        ta.invariants[ta.initial_state]
    )
    states = Tuple{State,TAConstraint{T}}[lD_initial]
    transitions = Pair{Tuple{State,TAConstraint{T}},Tuple{TAMessage,Tuple{State,TAConstraint{T}}}}[]
    

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
    while last < length(states) && max_iter > 0
        max_iter -=1
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
            if ignore_implicit_delays
                intersect!(up!(TimedAutomata.reset!(D′, arc.resets)), ta.invariants[l′])
            else
                intersect!(TimedAutomata.reset!(D′, arc.resets), ta.invariants[l′])
            end
            
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

    if max_iter==0
        println("Warning: reached maximum number of iterations!")
    end

    return (states=states, transitions=transitions)
end


function language(ts::TTS{State,T}) where {State,T}
    CombinedState = Tuple{State, TAConstraint{T}}
    
    function recurse(visited, state)
        if state in visited
            return [("",state)]
        end

        push!(visited, state)
        sequences = Tuple{String,Union{CombinedState,Nothing}}[]
        for trans in ts.transitions
            if trans[1]==state
                msg,target=trans[2]
                paths = recurse(deepcopy(visited), target)
                seqs=map(paths) do (tail,loop)
                    if loop==state
                        "($(msg)$(tail))*",nothing
                    else
                        "$(msg)$(tail)",loop
                    end
                end
                append!(sequences,seqs)
            end
        end

        return sequences
    end

    res = recurse(Vector{CombinedState}(), ts.states[1])
    return first.(res)
end

function Base.:|(tas::(TA{S, T} where S)...; names=[], separate_clocks=true) where {T}
    # Construct parallel state type
    _TT = Tuple{first.(getfield.(typeof(tas).parameters, :parameters))...}
    N = length(tas)
    TT = if isempty(names)
        names = 1:N
        _TT
    else
        names = Tuple(Symbol.(names))
        @assert length(names) == N
        NamedTuple{names,_TT}
    end

    clock_mappers = if separate_clocks
        [(s::Symbol)->Symbol("$(s)_$(i)") for i in 1:N]
    else
        fill(identity, N)
    end

    all_states = collect(product((ta.states for ta in tas)...))[:]
    all_edges = Vector{TAArc{TT,T}}()
    all_invariants = Dict{TT, TAConstraint{T}}()

    # compute clocks
    all_clocks = union((clock_mappers[i].(tas[i].clocks) for i in 1:N)...)
    
    # compute internal and external alphabet
    alphabet=Dict{Symbol, Symbol}()
    for tasᵢ in tas
        for arc in tasᵢ.arcs
            if haskey(alphabet, arc.message.symbol)
                # if we use the symbol as both in- and output, it must be an internal symbol (no mixed usage!)
                if alphabet[arc.message.symbol]==:input && arc.message.direction==MSG_OUT || alphabet[arc.message.symbol]==:output && arc.message.direction==MSG_IN
                    alphabet[arc.message.symbol] = :internal
                end
            else
                if arc.message.direction==MSG_IN
                    alphabet[arc.message.symbol] = :input
                elseif arc.message.direction == MSG_OUT
                    alphabet[arc.message.symbol] = :output
                end
            end
        end
    end
    
    # set the clock list for the combined automaton
    extend_clock_list(c::TAConstraint, automaton_id::Int) = set_clock_list!(deepcopy(c),all_clocks;mapping=clock_mappers[automaton_id])
    # compute the combined invariant for a given state from multiple automata
    compute_invariant(state) = intersect((extend_clock_list(tas[i].invariants[state[i]],i) for i in 1:N)...)
    
    initial_state = TT(first.(getfield.(tas,:states)))

    for state in all_states
        invariant = compute_invariant(state)
        all_invariants[state] = invariant

        @assert !isempty(invariant) "Invariant empty!"

        # step 1: make list of all available edges for this symbol
        # step 2: take ONE sender out of edge list
        # step 3: build every possible subset of the remaining edges (receivers and senders alike, one edge per process each)
            # step 3.a: compute combined guard = union of guards
            # step 3.b: ignore edge set if combined invariant of current state, guard & reset is incompatible with invariant of target state
        # step 4: for every valid edge set, check every superset
            # step 4.a: if the guard of the superset includes the guard of this set, remove this set 
            # step 4.b: else cut out the superset's guard region from this guard region <- this requires splitting this zone!!! -> multiple edges

        # step 1: make list of all available edges for this symbol
        options = DefaultDict{Symbol,NamedTuple{(:sender,:receiver,:silent),NTuple{3,Vector{Vector{Int}}}}}(
            ()->(sender=[Vector{Int}() for i in 1:N],receiver=[Vector{Int}() for i in 1:N],silent=[Vector{Int}() for i in 1:N]))
        for (i,(lᵢ,tasᵢ)) in enumerate(zip(state, tas))
            #check all processes' transitions starting in this state
            for (j,arc) in enumerate(tasᵢ.arcs)
                if arc.source==lᵢ
                    if arc.message.direction == MSG_IN
                        push!(options[arc.message.symbol].receiver[i],j)
                    elseif arc.message.direction == MSG_OUT
                        push!(options[arc.message.symbol].sender[i],j)
                    else
                        push!(options[arc.message.symbol].silent[i],j)
                    end
                end
            end
        end

        #see which transitions work
        for (symbol,transition_options) in pairs(options)
            potential_edges = TAArc{TT,T}[]

            edge_sets = if haskey(alphabet,symbol) && alphabet[symbol] in (:internal, :output)
                # step 2: fix exactly ONE process as sender
                tmp = []
                for (i,transition_ids) in enumerate(transition_options.sender)
                    if isempty(transition_ids)
                        continue
                    end
                    push!(tmp,product((j==i ? transition_ids : [nothing; transition_options.receiver[j]] for j=1:N)...))
                end
                flatten(tmp)
            else
                product(([nothing; transition_options.receiver[j]] for j=1:N)...)
            end

                
            # step 3: build every possible subset of the remaining edges (receivers and senders alike, one edge per process each)
            for edge_set in edge_sets
                # step 3.a: compute combined guard = union of guards
                guard = TAConstraint(TADiagonalConstraint{T}[],all_clocks)
                invalid = false
                tmp = Any[]
                resets = Set{Symbol}()
                for (k,edge_id) in enumerate(edge_set)
                    if isnothing(edge_id)
                        push!(tmp, state[k])
                    else
                        edge = tas[k].arcs[edge_id]
                        push!(tmp, edge.target)
                        union!(resets, clock_mappers[k].(edge.resets))
                        intersect!(guard, extend_clock_list(edge.guard,k))
                        if isempty(guard)
                            invalid=true
                            break
                        end
                    end
                end
                if invalid
                    continue
                end
                target_state = TT(tmp)

                # step 3.b: ignore edge if combined invariant of current state, guard & reset is incompatible with invariant of target state
                c = guard ∩ invariant
                if isempty(c)
                    continue
                end
                intersect!(TimedAutomata.reset!(c,resets),compute_invariant(target_state))
                if isempty(c)
                    continue
                end
                
                push!(potential_edges, TAArc(state, target_state, guard, TAMessage(MSG_OUT,symbol), resets))
            end
                
            # step 4: for every valid edge set, check every superset
            to_remove = zeros(Bool,length(potential_edges))
            for i in 1:length(potential_edges)
                if to_remove[i]
                    continue
                end
                for j in i+1:length(potential_edges)
                    if to_remove[j]
                        continue
                    end
                    edgeᵢ = potential_edges[i]
                    edgeⱼ = potential_edges[j]

                    # step 4.a: if the guard of the superset includes the guard of this set, remove this set 
                    if all((edgeᵢ.target .== state) .| (edgeᵢ.target .== edgeⱼ.target)) &&
                        (edgeⱼ.guard ⊆ edgeᵢ.guard)
                        to_remove[i] = true
                    elseif all((edgeⱼ.target .== state) .| (edgeᵢ.target .== edgeⱼ.target)) &&
                        (edgeᵢ.guard ⊆ edgeⱼ.guard)
                        to_remove[j] = true
                    elseif all((edgeᵢ.target .== state) .| (edgeᵢ.target .== edgeⱼ.target)) || 
                        all((edgeⱼ.target .== state) .| (edgeᵢ.target .== edgeⱼ.target))
                        # step 4.b: else cut out the superset's guard region from this guard region <- this requires splitting this zone!!! -> multiple edges
                        throw(ErrorException("Broadcast with non-trivially intersecting guards is not implemented yet"))
                    else
                        # sets don't have a super/subset relation at all
                    end
                end
            end
            deleteat!(potential_edges, to_remove)
            append!(all_edges, potential_edges)
        end
    end
    TA(all_states, initial_state, all_clocks, Set(keys(alphabet)), all_edges, all_invariants)
end


reachable_states(ta::TA) = first.(TTS(ta).states)
unreachable_states(ta::TA) = setdiff(ta.states, reachable_states(ta))

function prune!(ta::TA, prune_states=unreachable_states(ta))
    tts=TTS(ta)

    prune_at = ta.states .∈ Ref(prune_states)
    deleteat!(ta.states, prune_at)
    delete!.(Ref(ta.arcs_by_state), prune_states)

    # prune edges
    prune_at = zeros(Bool, length(ta.arcs))
    for (i,arc) in enumerate(ta.arcs)
        if arc.source ∈ prune_states || arc.target ∈ prune_states
            prune_at[i] = true
        end
    end
    deleteat!(ta.arcs, prune_at)

    for (_, arcs) in pairs(ta.arcs_by_state)
        prune_at = zeros(Bool, length(arcs))
        for (i,arc) in enumerate(arcs)
            if arc.source ∈ prune_states || arc.target ∈ prune_states
                prune_at[i] = true
            end
        end
        deleteat!(arcs, prune_at)
    end
    ta
end
prune(ta::TA) = prune!(deepcopy(ta))
