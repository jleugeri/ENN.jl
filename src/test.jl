using TimedAutomata, GLMakie, GraphMakie

states = ["S0", "S1", "S2", "S3"]
initial_state = "S0"
clocks = [:x,:y,:z]
symbols = Set(Symbol[])
arcs = Set([
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

states = ["start", "loop", "end"]
initial_state = "start"
clocks = [:x,:y]
symbols = Set(Symbol[])
arcs = Set([
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

zg2 = TTS(automaton2)

