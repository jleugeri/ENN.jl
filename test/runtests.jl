
using TimedAutomata, Test, LightGraphs, GraphPlot, Compose

tma = TMA(
    [
        [0,0],
        [0,1],
        [1,0],
        [1,1]
    ], 
    Vector{Int}[], 
    2, 
    [
        TMATransition([0,0],[0,1],:a,[1],@condition true),
        TMATransition([0,1],[1,0],:b,[2],@condition true),
        TMATransition([1,0],[1,1],:c,Int[],@condition c[1] < 1.0),
        TMATransition([1,1],[0,0],:d,Int[],@condition c[2] > 2.0),
    ], 
    [0,0]
)

testrun = TMARun(
    tma, 
    [
        (2.0,:a),
        (2.7,:b),
        (2.8,:b),
        (2.8,:c),
        (5.0,:d)
    ];
    drop_repetitions=true)


@test all(collect(last.(testrun)) .== [
    RunState([0, 0], [0.0, 0.0]),
    RunState([0, 1], [0.0, 2.0]),
    RunState([1, 0], [0.7, 0.0]),
    RunState([1, 1], [0.8, 0.1]),
    RunState([0, 0], [3.0, 2.3]),
])

g = to_graph(tma)

draw(SVG("test,svg", 16cm, 16cm), gplot(g.graph; nodelabel=string.(g.vertices), edgelabel=[string(g.edges[Pair(e)].a) for e ∈ edges(g.graph)]))