using TimedAutomata, GLMakie, GraphMakie
##

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

##
states = ["start", "loop", "end"]
initial_state = "start"
clocks = [:x,:y]
symbols = Set(Symbol[])
arcs = Vector([
    @arc("start","loop",true,TAMessage(),[:x,:y],Int64),
    @arc("loop","loop",x==10,TAMessage(),[:x],Int64),
    @arc("loop","end",y≥20,TAMessage(),[:x,:y],Int64)
])
invariants = Dict("start"=>@constraint(x-y==0,Int64), "loop"=>@constraint(x≤10,Int64))

automaton2 = TA(
    states,
    initial_state,
    clocks,
    symbols,
    arcs,
    invariants
)

##
τ = 1//1
automaton3 = TA(
    [(0,0), (0,1), (1,0), (1,1)],
    (0,0),
    [:x, :y],
    Set([:A; :B]),
    Vector([
        @arc((0,0), (1,0), true, TAMessage(MSG_IN, :A), [:x] ), 
        @arc((1,0), (0,0), x==τ, TAMessage(), [:x,:y]),
        @arc((1,0), (1,1), x<τ, TAMessage(MSG_IN, :B), [:y]),
        @arc((1,1), (0,0), y==τ, TAMessage(), [:x,:y]),
    ]),
    Dict(
        (1,0) => @constraint(x≤τ),
        (1,1) => @constraint(y≤τ)
    )
)
##
automaton4=TA(
    ["off","dim","bright"],
    "off",
    [:x],
    Set([:press]),
    Vector([
        @arc("off", "dim", true, TAMessage(MSG_IN, :press), [:x]),
        @arc("dim", "off", x>10//1, TAMessage(MSG_IN, :press), []),
        @arc("dim", "bright", x≤10//1, TAMessage(MSG_IN, :press), []),
        @arc("bright", "off", true, TAMessage(MSG_IN, :press), [])
    ]),
    Dict{String,TAConstraint{Rational{Int}}}()
)
zg4 = TTS(automaton4)


automaton5=TA(
    ["idle","t","study", "relax"],
    "idle",
    [:y],
    Set([:press]),
    Vector([
        @arc("t", "study", true, TAMessage(MSG_OUT, :press), []),
        @arc("study", "study", true, TAMessage(), []),
        @arc("study", "idle", true, TAMessage(MSG_OUT, :press), []),
        @arc("idle", "t", true, TAMessage(MSG_OUT, :press), [:y]),
        @arc("idle", "relax", true, TAMessage(MSG_OUT, :press), [:y]),
        @arc("relax", "idle", y>10//1, TAMessage(MSG_OUT, :press), [])
    ]),
    Dict(
        "t" => @constraint y<5//1
    )
)

automaton45 = automaton4 | automaton5
automaton44 = (automaton4 | automaton4)
automaton54 = automaton5 | automaton4

find_isomorphism(automaton45,automaton54)

automaton445 = |(automaton4, automaton4, automaton5)

##

automata = ("Example 1"=>automaton1, "Example 2"=>automaton2, "2 Neuron-Segments"=>automaton3, "Lamp"=>automaton4, "Student"=>automaton5, "Lamp + Student"=>automaton45, "2 Lamps + Student"=>automaton445);
for (name,automaton) in automata
    for condition in ("automaton", "pruned automaton", "zone graph")
        obj=if condition=="automaton"
            automaton
        elseif condition=="pruned automaton"
            prune(automaton)
        elseif condition=="zone graph"
            TTS(automaton)
        end

        f,ax,_ = graphplot(obj)
        hidedecorations!(ax)
        ax.title[] = "$(name) ($(condition))"
        save(joinpath("examples","figures","$(name)_$(condition).png"), f)
    end
end
