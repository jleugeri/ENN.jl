using Graphs, CairoMakie, GraphMakie

function GraphMakie.graphplot(ta::TA, args...; show_guards=true, show_invariants=false, show_resets=true, show_messages=(MSG_IN,MSG_OUT), nlabels_align=(:center,:bottom), nlabels_distance=25, elabels_distance=5, kwargs...)
    g = SimpleDiGraph(length(ta.states))
    
    nlabels = map(ta.states) do state
        inv = show_invariants ? "\n"*repr("text/plain", ta.invariants[state]) : ""
        repr("text/plain", state)*inv
    end

    elabels=Dict{Graphs.SimpleGraphs.SimpleEdge,String}()
    foreach(ta.arcs) do arc
        id1 = findfirst(==(arc.source), ta.states)
        id2 = findfirst(==(arc.target), ta.states)
        e = Graphs.SimpleGraphs.SimpleEdge(id1,id2)
        add_edge!(g, e)

        msg = arc.message.direction in show_messages ? repr("text/plain", arc.message) : ""
        clks = show_resets ? join(String.(arc.resets) .* ":=0", ", ") : ""
        guard = show_guards ? repr("text/plain", arc.guard) : ""

        elabels[e]=join(filter(!isempty,(msg,clks,guard)), "\n")
    end

    f,ax,p= graphplot(g, args...; nlabels, elabels=getindex.(Ref(elabels),edges(g)), nlabels_align, nlabels_distance, elabels_distance, kwargs...)

    deregister_interaction!(ax, :rectanglezoom)
    register_interaction!(ax, :ndrag, NodeDrag(p))
    return f,ax,p
end

function GraphMakie.graphplot(ts::TTS, args...; kwargs...)
    g = SimpleDiGraph(length(ts.states))
    
    nlabels = map(ts.states) do (state, zone)
        "State: "*repr("text/plain", state)*"\nZone: "*repr("text/plain", zone)
    end

    elabels=Dict{Graphs.SimpleGraphs.SimpleEdge,String}()
    foreach(ts.transitions) do (from, (msg, to))
        id1 = findfirst(==(from), ts.states)
        id2 = findfirst(==(to), ts.states)
        e = Graphs.SimpleGraphs.SimpleEdge(id1,id2)
        add_edge!(g, e)
        elabels[e]=repr("text/plain", msg)
    end

    f,ax,p= graphplot(g, args...; nlabels, elabels=getindex.(Ref(elabels),edges(g)), kwargs...)

    deregister_interaction!(ax, :rectanglezoom)
    register_interaction!(ax, :ndrag, NodeDrag(p))
    return f,ax,p
end