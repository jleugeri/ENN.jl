using TimedAutomata, CairoMakie, GraphMakie

τ = 1//1
Neuron = TA(
    [(0,0), (0,1), (1,0), (1,1)],
    (0,0),
    [:x, :y],
    Set([:A; :B]),
    Set([
        @arc((0,0), (1,0), true, TAMessage(MSG_IN, :A), [:x] ), 
        @arc((1,0), (0,0), x==τ, TAMessage(), []),
        @arc((1,0), (1,1), x<τ, TAMessage(MSG_IN, :B), [:y]),
        @arc((1,1), (0,0), y==τ, TAMessage(), []),
    ]),
    Dict(
        (1,0) => @constraint(x≤τ),
        (1,1) => @constraint(y≤τ)
    )
)

ts = TTS(Neuron)

f, ax, p = graphplot(ts, layout=GraphMakie.NetworkLayout.Stress())
save("test.svg", f)