
using TimedAutomata
using Test, LightGraphs, GraphPlot, Compose, IntervalArithmetic

#= 
tma = TMA(
    [
        [0,0],
        [0,1],
        [1,0],
        [1,1]
    ], 
    Vector{Int}[], 
    2, 
    Dict(
        ([0,0],:a)=>[TMATransition([0,1],[1],@condition true)],
        ([0,1],:b)=>[TMATransition([1,0],[2],@condition true)],
        ([1,0],:c)=>[TMATransition([1,1],Int[],@condition c[1] < 1.0), TMATransition([0,0],Int[],@condition c[1] == 1.0)],
        ([1,1],:d)=>[TMATransition([0,0],Int[],@condition c[2] > 2.0)],
    ),
    Dict([1,0]=>@condition c[1] < 1.0),
    [0,0]
)
@test collect_constants(tma) == [Set((1//1,)),Set((2//1,))]
=#

tma = TMA(
    [
        0,
        1,
        2,
        3
    ],
    [3],
    1, 
    Dict(
        (0,:msg)=>[TMATransition(1,[1],ClockRegion(Int,[nothing],[1]))],
        (1,:msg)=>[TMATransition(1,[1],ClockRegion(Int,[1],[1])),TMATransition(1,[1],ClockRegion(Int,[(1,nothing)],[1])),TMATransition(2,[1],ClockRegion(Int,[(0,1)],[1]))],
        (2,:alarm)=>[TMATransition(3,Int[],ClockRegion(Int,[nothing],[1]))]
    ),
    #Dict{Int,ClockRegion{1,Int}}(),
    Dict(2=>ClockRegion(Int,[(0,1)],[1])),
    0
)


@test collect_constants(tma) == [Set([0,1])]


tg = to_graph(tma)

draw(SVG("test_timed.svg", 16cm, 16cm), gplot(tg.graph; nodelabel=[string(tg.vertices[i][1]) for i ∈ vertices(tg.graph)], edgelabel=[string(tg.edges[Pair(e)][1]) for e ∈ edges(tg.graph)]))

ma = untime(tma)

ug = to_graph(ma)

draw(SVG("test_untimed.svg", 16cm, 16cm), gplot(ug.graph; nodelabel=[string(ug.vertices[i]) for i ∈ vertices(ug.graph)], edgelabel=[string(ug.edges[Pair(e)]) for e ∈ edges(ug.graph)]))

#=
for s ∈ ma.Σ
    println("State $(s):")
    for ((s′,a),ee) ∈ pairs(ma.E)
        if s′==s
            for e ∈ ee
                println("\t$(a)=>$(e.s)")
            end
        end
    end
end
=#

println.(language(ma, 20))
#= 
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

r = collect(last.(testrun))

@test length(r) == 5 && all(r .== [
    RunState([0, 0], [0.0, 0.0]),
    RunState([0, 1], [0.0, 2.0]),
    RunState([1, 0], [0.7, 0.0]),
    RunState([1, 1], [0.8, 0.1]),
    RunState([0, 0], [3.0, 2.3]),
])
 =#