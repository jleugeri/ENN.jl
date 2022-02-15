using ENN.Neurons, ENN.TimePetriNets, GLMakie, GraphMakie, DataStructures

counter_neurons = OrderedDict((Symbol("c_$(i)")=>Neuron((((),[:input],:dendrite),),[:trigger],:soma) for i in 1:10)...)
counter_ff_connections = OrderedDict((Symbol("c_$(i)")=>[(Symbol("c_$(i+1)"),:input, :excitatory)] for i in 1:9)...)
inputs = OrderedDict{Symbol,Vector{Tuple{Symbol,Symbol,Symbol}}}()
outputs = OrderedDict{Symbol,Symbol}()

net = NeuralNetwork(
    counter_neurons,    
    inputs,
    outputs,
    counter_ff_connections
)
#net_tpn = TPN(net)

f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
res=netplot!(ax, net, linewidth=2)
display(f)
