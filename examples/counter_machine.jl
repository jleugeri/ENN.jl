##
using ENN.Neurons, ENN.TimePetriNets, GLMakie, GraphMakie, DataStructures

function combine(all_neurons, all_connections)
    neurons = merge(all_neurons...)
    connections = all_connections[1]
    for cons in all_connections[2:end]
        for (key, con) in pairs(cons)
            if key in keys(connections)
                append!(connections[key], con)
            else
                connections[key] = con
            end
        end
    end
    (; neurons, connections)
end

function make_up_counter(name, size)
    neurons = OrderedDict((Symbol("$(name)_$(i)") => Neuron{Int}((((), [:input1, :input2], [:inh], :dendrite),), [:inc], [], :soma; axon_delay=2) for i in 1:size)...)
    connections = OrderedDict(
        (Symbol("$(name)_$(i)") => [(Symbol("$(name)_$(i+1)"), :input1)] for i in 1:size-1)...,
        Symbol("$(name)_init") => [(Symbol("$(name)_1"), :input1)],
        Symbol("$(name)_inc") => [(key, :inc) for key in keys(neurons)]
    )
    (; neurons, connections)
end

function make_loop(name, loop_label, loop_position, loop_origin, condition)
    neurons = OrderedDict(Symbol("$(name)_loop_$(loop_label)") => Neuron{Int}((((), [:condition], [], :dendrite),), [:trigger], [], :soma))
    connections = OrderedDict(
        Symbol("$(name)_loop_$(loop_label)") => [(Symbol("$(name)_$(loop_origin)"), :input2), (Symbol("$(name)_$(loop_position+1)"), :inh)],
        condition => [(Symbol("$(name)_loop_$(loop_label)"), :condition)]
    )
    (; neurons, connections)
end

function make_up_down_counter(name, size)
    neurons = OrderedDict((Symbol("$(name)_$(i)") => Neuron{Int}(
        (
            # left dendrite branch
            (
                (
                # left top branch
                    ((), [:input_12], [], :dendrite_12),
                ),
                [:input_11],
                :dendrite_11
            ),
            # right dendrite branch
            (
                (
                # right top branch
                    ((), [:input_22], [], :dendrite_22),
                ),
                [:input_21],
                :dendrite_21
            )
            # soma
        ), [:input_01, :input_02], [], :soma; axon_delay=2)
                           for i in 1:size)...
    )

    connections = OrderedDict(
        # from i to i+1, left top dendrite branch
        Symbol("$(name)_1") => [(Symbol("$(name)_2"), :input_12), (Symbol("$(name)_1"), :input_22)],
        (Symbol("$(name)_$(i)") => [(Symbol("$(name)_$(i+1)"), Symbol("$(name)_input_12")), (Symbol("$(name)_$(i-1)"), :input_22)] for i in 2:size-1)...,
        Symbol("$(name)_$(size)") => [(Symbol("$(name)_$(size-1)"), :input_22)],
        (Symbol("$(name)_dec") => [
            [(Symbol("$(name)_$(i)"), :input_02) for i in 1:size-1] ; 
            [(Symbol("$(name)_$(i)"), :input_21) for i in 1:size-1]
        ]),
        (Symbol("$(name)_inc") => [
            [(Symbol("$(name)_$(i)"), :input_01) for i in 2:size] ;
            [(Symbol("$(name)_$(i)"), :input_11) for i in 2:size]
            (Symbol("$(name)_1"), :input_01)
        ]),
        Symbol("$(name)_init") => [(Symbol("$(name)_1"), :input_11), (Symbol("$(name)_1"), :input_12)]
    )
    (; neurons, connections)
end
    

##

prog_counter=make_up_counter("program_counter", 5);
prog_loop = make_loop("program_counter", "loop_1", 3, 1, :condition);
prog = combine([prog_counter.neurons, prog_loop.neurons], [prog_counter.connections, prog_loop.connections]);

net = NeuralNetwork(
    prog.neurons,    
    [:program_counter_init,:program_counter_inc, :condition],
    Symbol[],
    prog.connections
);
f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
res=netplot!(ax, net)
display(f)


##
net = NeuralNetwork(
    res.neurons,    
    [:up1_init,:up1_inc,:up1_dec],
    Symbol[],
    res.connections
);

##

net = NeuralNetwork(
    program_counter_neurons,    
    [:init,:trigger],
    Symbol[],
    counter_ff_connections
)
##

macro program(exp)
    dump(exp)
end

@program begin
    INC(1)
    DEC(2)
    JZ(1,5)
end


##
f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
res=netplot!(ax, net)
net_tpn = res[:tpn][]
state = res[:state]
t_init_id = res[:tpn][].T_index[:input_init_start]
t_trigger_id = res[:tpn][].T_index[:input_trigger_start]
display(f)

trigger_period = 50
last_trigger = 0
sim_time = 0
##
function trigger_init()
    global sim_time = 0
    global last_trigger = 0
    fire!(net_tpn, state[], t_init_id)
    notify(state)
    fire_necessary!(net_tpn, state[])
    notify(state)
    println("Time: $(sim_time)")
end

function advance_trigger()
    next_trigger = last_trigger + trigger_period
    # execute all transitions until next trigger time
    while true
        τ = max_delay(net_tpn, state[])
        if sim_time + τ > next_trigger
            break
        end
        global sim_time += τ
        elapse!(net_tpn, state[], τ)
        fire_necessary!(net_tpn, state[])
        notify(state)
    end

    # advance to next trigger time
    elapse!(net_tpn, state[], next_trigger - sim_time)
    fire!(net_tpn, state[], t_trigger_id)
    notify(state)

    global sim_time = next_trigger
    global last_trigger = next_trigger
    # execute all transitions until next trigger time
    while true
        τ = max_delay(net_tpn, state[])
        if sim_time + τ > next_trigger
            break
        end
        global sim_time += τ
        elapse!(net_tpn, state[], τ)
        fire_necessary!(net_tpn, state[])
        notify(state)
    end

    println("Time: $(sim_time)")
end


##
trigger_init()
advance_trigger()
##
τ=max_delay(net_tpn,state[])
elapse!(net_tpn, state[], τ)
fire_necessary!(net_tpn, state[])
fire!(net_tpn, state[], t_trigger_id)
notify(state)
