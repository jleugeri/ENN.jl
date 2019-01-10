
using TimedAutomata
using Test, LightGraphs, GraphPlot, Compose, IntervalArithmetic

# Clock related test

@testset "All Tests" begin
    op1 = PartialOrder(3, (1,2)=>:<, (2,3)=>:(<=))
    op2 = PartialOrder(3, (2,1)=>:(>=), (1,3)=>:<)
    op3 = PartialOrder(3, (2,3)=>:<)
    ot1 = TotalOrder([1,2,3])
    ot2 = TotalOrder([1,1,2])
    ot3 = TotalOrder([1,2,2])
    
    t1 = Clk"0"
    t2 = Clk"1"
    t3 = Clk"3"
    t4 = Clk"(0,1)"
    t5 = Clk"[0,2)"
    t6 = Clk"(0,∞)"
    t7 = Clk"(∞,1]"
    t8 = Clk"[0,∞)"
    t9 = Clk""

    c1 = ClockRegion{3,Int}(1=>Clk"(0,1)", 2=>Clk"[0,1]", 3=>Clk"1"; order=op1)
    c2 = ClockRegion{3,Int}(1=>Clk"(0,1)", 2=>Clk"[0,1]", 3=>Clk"1"; order=op2)
    c3 = ClockRegion{3,Int}(1=>Clk"(0,1)", 2=>Clk"[0,1]"; order=op2)
    c4 = ClockRegion{3,Int}(1=>Clk"(0,1)", 2=>Clk"[0,1)"; order=op2)
    c5 = ClockRegion{3,Int}(1=>Clk"(0,1)", 2=>Clk"1"; order=op2)

    @testset "Equality tests" begin
        @test Clk"1" == Clk"1"
        @test Clk"1" != Clk"0"
        @test Clk"(0,1]" == Clk"(0,1]"
        @test Clk"(0,1]" != Clk"(0,1)"

        @test ClockRegion{1,Int64}(TimeSlice{Int64}[TimePoint{Int64}(0)], TotalOrder{1}([0])) == ClockRegion{1,Int64}(TimeSlice{Int64}[TimePoint{Int64}(0)], TotalOrder{1}([0]))
    end

    @testset "Order Relation tests" begin
        expected_results = [
            true  true  false false false false;
            false true  false false false false;
            false false true  false false false;
            true  true  true  true  false false;
            false true  true  false true  false;
            true  true  false false false true;
        ]
        results = zeros(Bool, size(expected_results))

        for (i,r1) ∈ enumerate((op1,op2,op3,ot1,ot2,ot3))
            for (j,r2) ∈ enumerate((op1,op2,op3,ot1,ot2,ot3))
                results[i,j] = (r1 ⊆ r2)
                @test results[i,j] == expected_results[i,j]
            end
        end

        failures = findall(results .!= expected_results)
        
        if ~isempty(failures)
            println("Failed for: $(failures)")
        end
    end

    @testset "Time slice tests" begin
        expected_results = [
            true  false false false true  false true  true  true;
            false true  false false true  true  true  true  true;
            false false true  false false true  false true  true;
            false false false true  true  true  true  true  true;
            false false false false true  false false true  true;
            false false false false false true  false true  true;
            false false false false false false true  false true;
            false false false false false false false true  true;
            false false false false false false false false true;
        ]
        results = zeros(Bool, size(expected_results))

        for (i,r1) ∈ enumerate((t1,t2,t3,t4,t5,t6,t7,t8,t9))
            for (j,r2) ∈ enumerate((t1,t2,t3,t4,t5,t6,t7,t8,t9))
                results[i,j] = (r1 ⊆ r2)
                @test results[i,j] == expected_results[i,j]
            end
        end

        failures = findall(results .!= expected_results)
        
        if ~isempty(failures)
            println("Failed for: $(failures)")
        end
    end

    @testset "Clock Region tests" begin
        expected_results = [
            true  true  true false false;
            false true  true false false;
            false false true false false;
            false false true true  false;
            false false true false true
        ]
        results = zeros(Bool, size(expected_results))


        for (i,r1) ∈ enumerate((c1,c2,c3,c4,c5))
            for (j,r2) ∈ enumerate((c1,c2,c3,c4,c5))
                results[i,j] = (r1 ⊆ r2)
                @test results[i,j] == expected_results[i,j]
            end
        end

        failures = findall(results .!= expected_results)

        if ~isempty(failures)
            println("Failed for: $(failures)")
        end
    end

end


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
        (0,:msg)=>[TMATransition(1,[1],ClockRegion{1,Int}())],
        (1,:msg)=>[TMATransition(1,[1],ClockRegion{1,Int}(1=>Clk"1")),TMATransition(1,[1],ClockRegion{1,Int}(1=>Clk"[1,∞)")),TMATransition(2,[1],ClockRegion{1,Int}(1=>Clk"(0,1)"))],
        (2,:alarm)=>[TMATransition(3,Int[],ClockRegion{1,Int}())]
    ),
    #Dict{Int,ClockRegion{1,Int}}(),
    Dict(2=>ClockRegion{1,Int}(1=>Clk"(0,1)")),
    0
)

@test collect_constants(tma) == [Set([0,1])]

ma = untime(tma)

L = language(ma, 20)
println.(L)

tg = to_graph(tma)
#ug = to_graph(ma)

draw(SVG("test_timed.svg", 16cm, 16cm), gplot(tma))
draw(SVG("test_untimed.svg", 16cm, 16cm), gplot(ma))

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