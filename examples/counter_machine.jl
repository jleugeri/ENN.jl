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

function prune!(neurons, connections)
    all_targets = Tuple{Symbol,Symbol}[]
    for connection in values(connections)
        append!(all_targets, connection)
    end
    unique!(all_targets)

    for (key,neuron) in pairs(neurons)
        for segment in neuron.all_segments
            to_delete = [i for (i,inp) in enumerate(segment.exc_inputs) if (key, inp) ∉ all_targets]
            deleteat!(segment.exc_inputs, to_delete)
            
            to_delete = [i for (i,inp) in enumerate(segment.inh_inputs) if (key, inp) ∉ all_targets]
            deleteat!(segment.inh_inputs, to_delete)            
        end
    end
end

function make_up_counter(name, size, init_name=Symbol("$(name)_init"), inc_name=Symbol("$(name)_inc"))
    neurons = OrderedDict((Symbol("$(name)_$(i)") => Neuron{Int}((((), [:input1, :input2], [:inh], :dendrite),), [:inc], [], :soma; axon_delay=2) for i in 1:size)...)
    connections = OrderedDict(
        (Symbol("$(name)_$(i)") => [(Symbol("$(name)_$(i+1)"), :input1)] for i in 1:size-1)...,
        init_name => [(Symbol("$(name)_1"), :input1)],
        inc_name => [(key, :inc) for key in keys(neurons)]
    )
    (; neurons, connections)
end

function make_loop(name, loop_label, loop_position, loop_origin, condition_name=Symbol("$(name)_condition"))
    neurons = OrderedDict(Symbol("$(name)_loop_$(loop_label)") => Neuron{Int}((((), [:condition], [], :dendrite),), [:trigger], [], :soma; spike_duration=1))
    connections = OrderedDict(
        Symbol("$(name)_loop_$(loop_label)") => [(Symbol("$(name)_$(loop_origin)"), :input2), (Symbol("$(name)_$(loop_position+1)"), :inh)],
        Symbol("$(name)_$(loop_position)") => [(Symbol("$(name)_loop_$(loop_label)"),:trigger)],
        condition_name => [(Symbol("$(name)_loop_$(loop_label)"), :condition)]
    )
    (; neurons, connections)
end

function make_up_down_counter(name, size, init_name=Symbol("$(name)_init"), inc_name=Symbol("$(name)_inc"), dec_name=Symbol("$(name)_dec"), refresh_name=Symbol("$(name)_refresh"))
    neurons = OrderedDict(
        Symbol("$(name)_$(:interneuron)")=>Neuron{Int}((),[:input_dec, :input_inc, :input_refresh, :input_init],[],:soma; axon_delay=1),
        (Symbol("$(name)_$(i)") => Neuron{Int}(
        (
            # left dendrite branch
            (
                (
                # left top branch
                    ((), [:input_121], [:input_122], :dendrite_12),
                ),
                [:input_111], [:input_112],
                :dendrite_11
            ),
            # center dendrite branch
            (
                (
                # left top branch
                    ((), [:input_321], [:input_322], :dendrite_32),
                ),
                [:input_311], [:input_312],
                :dendrite_31
            ),
            # right dendrite branch
            (
                (
                # right top branch
                    ((), i==1 ? [:input_221,:input_223] : [:input_221], [:input_222], :dendrite_22),
                ),
                [:input_211], [:input_212],
                :dendrite_21
            )
            # soma
        ), i==1 ? [:input_01, :input_02, :input_03, :input_04] : [:input_01, :input_02, :input_03], [], :soma; axon_delay=10)
                           for i in 1:size)...
    )

    connections = OrderedDict(
        # from i to i+1, left top dendrite branch
        
        (Symbol("$(name)_$(i)") => [
            (Symbol("$(name)_$(i+1)"), :input_121), 
            (Symbol("$(name)_$(i)"), :input_321), 
            (Symbol("$(name)_$(i-1)"), :input_221)
        ] for i in 2:size-1)...,
        dec_name => [
            [(Symbol("$(name)_$(i)"), :input_02) for i in 1:size-1]; 
            [(Symbol("$(name)_$(i)"), :input_211) for i in 1:size-1];
             (Symbol("$(name)_$(:interneuron)"), :input_dec)
        ],
        inc_name => [
            [(Symbol("$(name)_$(i)"), :input_01) for i in 2:size];
            [(Symbol("$(name)_$(i)"), :input_111) for i in 2:size]; 
             (Symbol("$(name)_1"), :input_01);
             (Symbol("$(name)_$(:interneuron)"), :input_inc)
        ],
        refresh_name => [
            [(Symbol("$(name)_$(i)"), :input_03) for i in 1:size]; 
            [(Symbol("$(name)_$(i)"), :input_311) for i in 1:size];
             (Symbol("$(name)_$(:interneuron)"), :input_refresh) 
        ],
        Symbol("$(name)_$(:interneuron)") => [
            [(Symbol("$(name)_$(i)"), :input_222) for i in 1:size];
            [(Symbol("$(name)_$(i)"), :input_212) for i in 1:size];
            [(Symbol("$(name)_$(i)"), :input_122) for i in 1:size];
            [(Symbol("$(name)_$(i)"), :input_112) for i in 1:size];
            [(Symbol("$(name)_$(i)"), :input_322) for i in 1:size];
            [(Symbol("$(name)_$(i)"), :input_312) for i in 1:size]
        ],
        init_name => [
            (Symbol("$(name)_1"), :input_04), 
            (Symbol("$(name)_1"), :input_111), 
            (Symbol("$(name)_1"), :input_121),
            (Symbol("$(name)_$(:interneuron)"), :input_init) 
        ],
    )
    if size>1
        push!(connections,
            Symbol("$(name)_1") => [
                (Symbol("$(name)_2"), :input_121), 
                (Symbol("$(name)_1"), :input_321), 
                (Symbol("$(name)_1"), :input_223)
            ],
            Symbol("$(name)_$(size)") => [
                (Symbol("$(name)_$(size)"), :input_321), 
                (Symbol("$(name)_$(size-1)"), :input_221)
            ]
        )
    else
        push!(connections,
            Symbol("$(name)_1") => [
                (Symbol("$(name)_1"), :input_321), 
                (Symbol("$(name)_1"), :input_223)
            ]
        )
    end
    (; neurons, connections)
end
    
function make_relay_neuron(name, sources)
    neurons = OrderedDict(name => Neuron{Int}((), [Symbol("input_$(src)") for src in sources], [], :soma))

    connections = OrderedDict(
        (src => [(name, Symbol("input_$(src)"))] for src in sources)...
    )
    (; neurons, connections)
end

function make_input_trigger(netplot, input_name)
    net_tpn = netplot[:tpn][]
    state = netplot[:state]
    t_id = net_tpn.T_index[Symbol("input_$(input_name)_start")]

    function trigger()
        fire!(net_tpn, state[], t_id)
        fire_necessary!(net_tpn, state[])
        notify(state)
    end
end

function make_time_advance(netplot, trigger_period)
    net_tpn = netplot[:tpn][]
    state = netplot[:state]

    function advance_trigger()
        sim_time = 0
        # execute all transitions until next trigger time
        while true
            τ = max_delay(net_tpn, state[])
            if sim_time + τ > trigger_period
                break
            end
            sim_time += τ
            elapse!(net_tpn, state[], τ)
            fire_necessary!(net_tpn, state[])
        end
    
        # advance to next trigger time
        elapse!(net_tpn, state[], trigger_period - sim_time)
        notify(state)
    end
end


##

prog_counter=make_up_counter(:program_counter, 5, :init);
prog_loop = make_loop(:program_counter, "loop_1", 3, 1, :condition);
counter1 = make_up_down_counter(:counter1, 3, :init);
relay_inc1 = make_relay_neuron(:counter1_inc, [:program_counter_2, :program_counter_3])
relay_dec1 = make_relay_neuron(:counter1_dec, Symbol[:program_counter_4,:program_counter_5])
relay_refresh1 = make_relay_neuron(:counter1_refresh, Symbol[:program_counter_1])
prog = combine([prog_counter.neurons, prog_loop.neurons, relay_inc1.neurons, relay_dec1.neurons, relay_refresh1.neurons, counter1.neurons], [prog_counter.connections, counter1.connections, prog_loop.connections, relay_inc1.connections, relay_dec1.connections, relay_refresh1.connections]);
prune!(prog.neurons, prog.connections)

net = NeuralNetwork(
    prog.neurons,    
    [:init,:program_counter_inc, :condition],
    Symbol[],
    prog.connections
);


f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
res=netplot!(ax, net)
display(f)
trigger_init = make_input_trigger(res, :init)
trigger_condition = make_input_trigger(res, :condition)
trigger_inc = make_input_trigger(res, :program_counter_inc)
advance_cycle = make_time_advance(res, 75)

## Run through without looping
trigger_init()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()

## Run through with looping
trigger_init()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
trigger_condition()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()

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
function parse_program(expr::Expr, counters...)
    neurons = []
    connections = []
    groups = Vector{Symbol}[]

    counter_sizes = OrderedDict(counters...)

    counter_ops = OrderedDict()
    for key in keys(counter_sizes)
        counter_ops[key] = (inc=Int[],dec=Int[])
    end
    # go through program and create loop neurons as we go
    i=0
    for instruction in expr.args
        if (instruction isa Expr) && (instruction.head == :call)
            i+=1
            if(length(instruction.args)==3 && instruction.args[1] == :JZ)
                (counter, loopback) = instruction.args[2:3]
                (n,c) = make_loop(:program_counter,Symbol("inst_$(i)_loop"),i+1,loopback,Symbol("counter_$(counter)_1"))
                if !haskey(counter_ops, counter)
                    counter_ops[counter] = (inc=Int[],dec=Int[])
                end
                push!(neurons,n)
                push!(connections,c)
            elseif (length(instruction.args)==2 && instruction.args[1] == :INC)
                counter = instruction.args[2]
                if !haskey(counter_ops, counter)
                    counter_ops[counter] = (inc=Int[],dec=Int[])
                end
                push!(counter_ops[counter].inc, i)
            elseif (length(instruction.args)==2 && instruction.args[1] == :DEC)
                counter = instruction.args[2]
                if !haskey(counter_ops, counter)
                    counter_ops[counter] = (inc=Int[],dec=Int[])
                end
                push!(counter_ops[counter].dec, i)
            else
                error("Cannot parse instruction '$(instruction)'.")
            end
        end
    end
    num_instructions = i+2 # +1 for the start instruction, +1 to have a buffer if the program ends in a JZ
    # make program counter
    (n,c) = make_up_counter(:program_counter, num_instructions, :init, :step)
    pushfirst!(neurons, n)
    pushfirst!(connections, c)

    prog = combine(neurons, connections)
    push!(groups, collect(keys(prog.neurons)))

    # make up-down counters & relay neurons
    for (counter,ops) in pairs(counter_ops)
        sz = get(counter_sizes, counter, length(ops.inc))+1
        cnt = make_up_down_counter(Symbol("counter_$(counter)"), sz, :init)
        rly_inc = make_relay_neuron(Symbol("counter_$(counter)_inc"),[Symbol("program_counter_$(i)") for i in ops.inc])
        rly_dec = make_relay_neuron(Symbol("counter_$(counter)_dec"),[Symbol("program_counter_$(i)") for i in ops.dec])
        rly_refresh = make_relay_neuron(Symbol("counter_$(counter)_refresh"),[Symbol("program_counter_$(i)") for i in 1:num_instructions if ((i ∉ ops.inc) && (i∉ ops.dec))])

        push!(neurons, cnt.neurons, rly_inc.neurons, rly_dec.neurons, rly_refresh.neurons)
        push!(connections, cnt.connections, rly_inc.connections, rly_dec.connections, rly_refresh.connections)
        push!(groups, collect(keys(cnt.neurons)))
        append!(groups[1],collect(keys(rly_inc.neurons)),collect(keys(rly_dec.neurons)),collect(keys(rly_refresh.neurons)))
    end
    
    prog = combine(neurons, connections)
    prune!(prog.neurons, prog.connections)

    net=NeuralNetwork(
            (prog.neurons),
            [:init,:step],
            Symbol[],
            (prog.connections)
    )
    return net,groups
end

macro program(expr, counters...)
    net,groups = parse_program(expr, eval.(counters)...)
    return :($net,$groups)
end

net,groups = @program begin
    INC(1)
    INC(1)
    INC(1)
    INC(1)
    INC(1)
    INC(2)
    INC(2)
    INC(2)
    INC(2)
    INC(2)
    DEC(2)
    DEC(2)
    DEC(2)
    DEC(2)
    DEC(2)
    DEC(1)
    DEC(1)
    DEC(1)
    DEC(1)
    DEC(1)
end 1=>5 2=>5












net,groups = @program begin
    INC(1)
    INC(1)
    INC(1)
    INC(1)
    INC(1)
    INC(2)
    INC(2)
    INC(2)
    JZ(2,13)
    DEC(2)
    DEC(1)
    JZ(3,9)
    DEC(3)
end 1=>5 2=>5 3=>1











begin
f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
res=netplot!(ax, net, groups)
hidedecorations!(ax)
display(f)
trigger_init = make_input_trigger(res, :init)
trigger_inc = make_input_trigger(res, :step)
advance_cycle = make_time_advance(res, 75)
advance_minimal = make_time_advance(res, 1)
end

## Run through without looping
begin
sleep(2)
trigger_init()
for i in 1:75
    sleep(0.001)
    advance_minimal()
end
for j in 1:15
    trigger_inc()
    for i in 1:75
        sleep(0.001)
        advance_minimal()
    end
end
end
advance_cycle()
trigger_inc()
## Run through with looping
trigger_init()
advance_cycle()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
trigger_condition()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()
trigger_inc()
advance_cycle()


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
