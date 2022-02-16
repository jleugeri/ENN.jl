using ENN.Neurons, ENN.TimePetriNets, GLMakie, GraphMakie, DataStructures

counter_neurons = OrderedDict((Symbol("c_$(i)")=>Neuron((((),[:input],:dendrite),),[:trigger],:soma) for i in 1:10)...)
counter_ff_connections = OrderedDict((Symbol("c_$(i)")=>[(Symbol("c_$(i+1)"),:input, excitatory::SynapseType)] for i in 1:9)...)
inputs = OrderedDict(
    :init => [(:c_1, :input, excitatory::SynapseType)],
    :trigger => [(key, :trigger, excitatory::SynapseType) for key in keys(counter_neurons)]
)
outputs = OrderedDict{Symbol,Symbol}()

net = NeuralNetwork(
    counter_neurons,    
    inputs,
    outputs,
    counter_ff_connections
)

f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
res=netplot!(ax, net, linewidth=2)
display(f)
##

function set_colors_to_state!(netplot, state::Observable{<:TPNState}; input_color=(active=:red,passive=:black), axon_color=input_color, dendrite_color=(active=:red,passive=:silver), spine_color=dendrite_color)
    on(state) do state
        net = netplot[:net][]
        
        # Process input axons
        for (name,synapses) in pairs(net.inputs)
            netplot[:input_colors][][name] = input_color.passive
            for synapse in synapses
                for dendrite in synapse.target_dendrites
                    (neuron,input) = synapse.target
                    neuron_name = neuron[].name[]
                    tgt_idx = findfirst(==(Symbol("syn_$(neuron_name)_$(dendrite.name)_$(input)_trigger")), net_tpn.P)
                    if state.m[tgt_idx] > 0
                        netplot[:input_colors][][name] = input_color.active
                        break
                    end
                end
            end
        end
        
        # Process internal axons
        for neuron_name in keys(netplot[:axon_colors][])
            soma_name = net.neurons[neuron_name].dendrite.name
            tgt_idx = findfirst(==(Symbol("seg_$(neuron_name)_$(soma_name)_on")), net_tpn.P)
            netplot[:axon_colors][][neuron_name] = 
                state.m[tgt_idx] > 0 ? axon_color.active : axon_color.passive
        end

        # Process spines
        for (neuron_name, dendrite_name, input_name) in keys(netplot[:spine_colors][])
            tgt_idx = findfirst(==(Symbol("syn_$(neuron_name)_$(dendrite_name)_$(input_name)_on")), net_tpn.P)
            netplot[:spine_colors][][(neuron_name, dendrite_name, input_name)] = 
                state.m[tgt_idx] > 0 ? spine_color.active : spine_color.passive
        end

        # Process dendrites
        for (neuron_name, dendrite_name) in keys(netplot[:dendrite_colors][])
            tgt_idx = findfirst(==(Symbol("seg_$(neuron_name)_$(dendrite_name)_on")), net_tpn.P)
            netplot[:dendrite_colors][][(neuron_name, dendrite_name)] = 
                state.m[tgt_idx] > 0 ? dendrite_color.active : dendrite_color.passive
        end

        notify(netplot[:spine_colors])
        notify(netplot[:axon_colors])
        notify(netplot[:input_colors])
        notify(netplot[:dendrite_colors])
    end

    notify(state)
end

##

net_tpn = TPN(net)
state = Observable(net_tpn.x₀)
set_colors_to_state!(res, state)

##
fire!(net_tpn, state[], 1)
notify(state)
##
τ=max_delay(net_tpn,state[])
elapse!(net_tpn, state[], τ)
fire_necessary!(net_tpn, state[])
notify(state)

##

τ=max_delay(net_tpn,state[])
elapse!(net_tpn, state[], τ)
fire_necessary!(net_tpn, state[])
notify(state)
fire!(net_tpn, state[], 2)
notify(state)
