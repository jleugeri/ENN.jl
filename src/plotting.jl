using Graphs, CairoMakie, GraphMakie

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

    return graphplot(g, args...; nlabels, elabels, kwargs...)
end