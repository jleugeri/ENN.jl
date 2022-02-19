using ENN.Neurons, ENN.TimePetriNets, GLMakie, GraphMakie, DataStructures

counter_neurons = OrderedDict((Symbol("c_$(i)")=>Neuron{Int}((((),[:input],:dendrite),),[:trigger],:soma; axon_delay=2) for i in 1:10)...)
counter_ff_connections = OrderedDict(
    (Symbol("c_$(i)")=>[(Symbol("c_$(i+1)"),:input)] for i in 1:9)...,
    :init => [(:c_1, :input)],
    :trigger => [(key, :trigger) for key in keys(counter_neurons)]
)

net = NeuralNetwork(
    counter_neurons,    
    [:init,:trigger],
    Symbol[],
    counter_ff_connections
)

##
f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
res=netplot!(ax, net)
net_tpn = res[:tpn][]
state = res[:state]
t_init_id = res[:tpn][].T_index[:input_init_start]
t_trigger_id = res[:tpn][].T_index[:input_trigger_start]
display(f)

##
fire!(net_tpn, state[], t_init_id)
notify(state)

##
τ=max_delay(net_tpn,state[])
elapse!(net_tpn, state[], τ)
fire_necessary!(net_tpn, state[])
τ=max_delay(net_tpn,state[])
elapse!(net_tpn, state[], τ)
fire!(net_tpn, state[], t_trigger_id)
notify(state)
