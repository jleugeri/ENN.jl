using Graphs, CairoMakie, GraphMakie

function GraphMakie.graphplot(ta::TA, args...; kwargs...)
    g = SimpleDiGraph(length(ta.states))
    
    nlabels = map(ta.states) do state
        "State: "*repr("text/plain", state)#*"\Guard: "*repr("text/plain", ta.invariants[state])
    end

    elabels = map(ta.arcs) do arc
        id1 = findfirst(==(arc.source), ta.states)
        id2 = findfirst(==(arc.target), ta.states)
        add_edge!(g, id1, id2)

        msg = repr("text/plain", arc.message)
        clks = join(",", ta.resets)
        guard = repr("text/plain", arc.guard)

        "$(msg), {$(clks)}\n$(guard)"
    end

    f,ax,p= graphplot(g, args...; nlabels, elabels, kwargs...)

    deregister_interaction!(ax, :rectanglezoom)
    register_interaction!(ax, :ndrag, NodeDrag(p))
    return f,ax,p
end

function GraphMakie.graphplot(ts::TTS, args...; kwargs...)
    g = SimpleDiGraph(length(ts.states))
    
    nlabels = map(ts.states) do (state, zone)
        "State: "*repr("text/plain", state)*"\nZone: "*repr("text/plain", zone)
    end

    elabels = map(ts.transitions) do (from, (msg, to))
        id1 = findfirst(==(from), ts.states)
        id2 = findfirst(==(to), ts.states)
        add_edge!(g, id1, id2)
        repr("text/plain", msg)
    end

    f,ax,p= graphplot(g, args...; nlabels, elabels, kwargs...)

    deregister_interaction!(ax, :rectanglezoom)
    register_interaction!(ax, :ndrag, NodeDrag(p))
    return f,ax,p
end