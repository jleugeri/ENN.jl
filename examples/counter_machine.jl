using ENN.Neurons, ENN.TimePetriNets, GLMakie, GraphMakie
GLMakie.activate()

counter_neurons = Dict((Symbol("c_$(i)")=>Neuron((((),[:input],:dendrite),),[:trigger],:soma) for i in 1:10)...)
counter_ff_connections = Dict((Symbol("c_$(i)")=>[(Symbol("c_$(i+1)"),:dendrite,:input)] for i in 1:9)...)
counter_trig_connections = Dict()

net = NeuralNetwork(
    counter_neurons,    
    counter_ff_connections
)
net_tpn = TPN(net)


f,ax,p = graphplot(net_tpn); f
