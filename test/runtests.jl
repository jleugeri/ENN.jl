using TimedAutomata, Test

@testset "All tests" begin
@testset "Bounds tests" begin
    for T in (Rational{Int},Int,Float64)
        @testset "Time-type $(T)" begin
            b1 = TABound(true,typemin(T))
            b2 = TABound(true,-one(T))
            b3 = TABound(false,-one(T))
            b4 = TABound(true,zero(T))
            b5 = TABound(true,one(T))
            b6 = TABound(false,one(T))
            b7 = TABound(false,typemax(T))

            @test b1≤b2≤b3≤b4≤b5≤b6≤b7
            @test b1<b2<b3<b4<b5<b6<b7
            @test !(b6<b5)

            @test isinf.([b1,b2,b3,b4,b5,b6,b7]) == [true,false,false,false,false,false,true]

            b2.value = typemin(T)
            @test isinf(b2)

            @test b1==b2 && hash(b1)==hash(b2) && isequal(b1,b2)

            @test_throws ArgumentError b1+b7

            @test b1+b2 == typemin(TABound{T})
            @test b6+b7 == typemax(TABound{T})

            @test b3+b4 == TABound(true,-one(T))
            @test b3+b6 == TABound(false,zero(T))
            @test b5+b6 == TABound(true,one(T)+one(T))

            @test b3-b4 == TABound(true,-one(T))
            @test b3-b6 == TABound(false,-one(T)-one(T))
            @test b5-b6 == TABound(true,zero(T))

            @test repr.(Ref("text/plain"),[b1,b2,b3,b4,b5,b6,b7]) == ["-∞","-∞","≤ -$(one(T))","< $(zero(T))","< $(one(T))","≤ $(one(T))","∞"]
        end
    end
end

@testset "Constraint tests" begin
    c1 = TADiagonalConstraint(:x, :y, TABound(true, 5))
    c2 = TADiagonalConstraint(:y, :x, TABound(false, -5))
    c3 = TADiagonalConstraint(nothing, :x, TABound(false, -2))
    c4 = TADiagonalConstraint(:x, :y, TABound(false, -3))
    c5 = TADiagonalConstraint(:x, :y, TABound(false, 7))
    c6 = TADiagonalConstraint(:x, nothing, TABound(false, 5))

    @test repr.(Ref("text/plain"),[c1,c2,c3]) == ["x-y < 5","y-x ≤ -5","x ≥ 2"]
    @test c1 != c2
    @test (c1 == !c2) && isequal(c1,!c2) && hash(c1) == hash(!c2)

    tac = TAConstraint([c2,c3],[:x,:y])
    @test repr("text/plain", tac) == "2 ≤ x\ny-x ≤ -5"

    @test !isempty(tac)
    @test TimedAutomata.close(tac).D[1,2].value==-5
    @test isempty(tac ∩ c4)
    @test repr("text/plain", tac ∩ c4) == "∅"
    @test !isempty(tac ∩ c5 ∩ c6)
    @test repr("text/plain", tac ∩ c5 ∩ c6) == "x == 5\ny == 0\ny-x == -5"

    @test tac ∩ c5 ∩ c6 == tac ∩ TAConstraint([c5,c6],[:x,:y])
    @test tac ∩ c5 ∩ c6 ⊆ tac ∩ c5 ⊆ tac
    @test !(tac ∩ c5 ⊆ tac ∩ c6)

    tac2 = copy(tac)
    set_clock_list!(tac2, [:y,:x,:z])
    @test tac2.clocks == [:y,:x,:z] && tac2.D[[3,2],[3,2]] == tac.D[[2,3],[2,3]]

    cs = get_diagonal_constraints(tac)
    @test [c2,c3] ⊆ cs
    tac3 = TAConstraint(cs, [:x,:y])
    @test tac3 == tac && tac ⊆ tac3 && tac3 ⊆ tac
    
    tac4 = @constraint (y-x ≤ -5 && y ≥ 0 && 2 ≤ x) Int64
    set_clock_list!(tac4, [:x,:y])
    @test tac4 == tac3
end

@testset "Automata tests" begin
    states = ["S0", "S1", "S2", "S3"]
    initial_state = "S0"
    clocks = [:x,:y,:z]
    symbols = Set(Symbol[])
    arcs = Vector([
        @arc("S0","S1",true,TAMessage(),[:z],Int64),
        @arc("S1","S2",y>2,TAMessage(),[:y],Int64),
        @arc("S2","S3",x-z<1 && z-y<1,TAMessage(),[],Int64)
    ])
    invariants = Dict("S0"=>@constraint x-y==0 && y-z==0 && z-x==0 Int64)

    automaton1 = TA(
        states,
        initial_state,
        clocks,
        symbols,
        arcs,
        invariants
    )

    zg1 = TTS(automaton1)
end

end