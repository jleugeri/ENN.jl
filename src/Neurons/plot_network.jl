using Base.Iterators, LinearAlgebra
using DataStructures, IntervalSets
using ENN.Neurons, ENN.TimePetriNets

export set_colors_to_state!

# Optional dependencies from Requires.jl:
using .GLMakie

include("plot_neuron.jl")
include("plot_wires.jl")


@recipe(NetPlot, net) do scene
    Attributes(
        linewidth=3,
        csegs=20,
        axon_margin = 0.05,
        bridge_radius = 0.025,
        neuron_margin = 0.1,
        default_axon_color=(active=:red,passive=:black),
        default_input_color=(active=:red,passive=:black),
        default_dendrite_color=(active=:red,passive=:silver),
        default_spine_color=(active=:red,passive=:silver),
        input_colors = DefaultDict{Symbol,Union{Symbol,RGBf,RGBAf}}(:black),
        axon_colors = DefaultDict{Symbol,Union{Symbol,RGBf,RGBAf}}(:black),
        dendrite_colors = DefaultDict{Tuple{Symbol,Symbol},Union{Symbol,RGBf,RGBAf}}(:silver),
        spine_colors = DefaultDict{Tuple{Symbol,Symbol,Symbol},Union{Symbol,RGBf,RGBAf}}(:silver),
    )
end

function Makie.plot!(netplot::NetPlot)
    net = netplot[1][]
    
    tpn = TPN(net)
    state = Observable(deepcopy(tpn.x₀))

    (; axon_margin, bridge_radius,csegs,linewidth,neuron_margin) = netplot.attributes
    bridge_radius = bridge_radius[]
    axon_margin = axon_margin[]
    csegs=csegs[]
    linewidth=linewidth[]
    neuron_margin=neuron_margin[]

    neurons = OrderedDict{Symbol,Combined{neuronplot, Tuple{valtype(net.neurons), Point{2, Float32}}}}()
    widths = OrderedDict{Symbol,Observable{Vector{Float32}}}()
    for (name,neuron) in pairs(net.neurons)
        n = neuronplot!(netplot, neuron, Point2f(0.0,0.0); bridge_radius,csegs,inputs_kwargs=Dict(:linewidth=>linewidth))
        neurons[name]=n
        widths[name]=n[:x_extent]
        netplot[Symbol(name)] = n
    end

    onany(values(widths)...) do _widths...
        offset = neuron_margin
        for (name,(x_min,x_max)) in zip(keys(widths),_widths)
            offset -= x_min
            neurons[name][:offset][] = Point2f(offset, 0.0)
            offset += x_max + neuron_margin
            notify(neurons[name][:offset])
        end
    end

    notify(first(values(widths)))

    # layout axons
    horizontal_axons=NamedTuple{(:name,:x_range, :tgts),Tuple{Symbol,IntervalSets.ClosedInterval{Float32}, Vector{Point2f}}}[]
    for (name,neuron) in pairs(neurons)
        targets = get(net.axons, name, [])
        x_min = Inf
        x_max = -Inf
        x_root = neuron[:soma][][1]
        
        tgts = Point2f[]
        
        # add target output if any
        if name in net.outputs
            x_min = 0
        end

        # add target neurons
        for target in targets
            x_tgt = neurons[target.target[1][].name[]][:inputs_ports][][target.target[2]][1]
            x_min = min(x_min, x_tgt, x_root)
            x_max = max(x_max, x_tgt, x_root)
            push!(tgts, neurons[target.target[1][].name[]][:inputs_ports][][target.target[2]])
        end
        sort!(tgts; by=first)
        push!(horizontal_axons, (name=name, x_range=x_min..x_max, tgts=tgts))
    end

    # simple heuristic sorting
    sort!(horizontal_axons; 
        lt=(x1,x2)->(x1.x_range ∈ x2.x_range) || 
            (maximum(x1.x_range) - minimum(x1.x_range)) ≤ (maximum(x2.x_range) - minimum(x2.x_range))
    )

    row_blockers = Vector{IntervalSets.ClosedInterval{Float32}}[]
    
    wires = OrderedDict{Symbol,WireNet}()
    row_y = 0
    axon_y = row_y
    
    input_ports = OrderedDict{Symbol,Point2f}()
    output_ports = OrderedDict{Symbol,Point2f}()
    
    # collect axons from input to neurons
    for name in net.inputs
        if name in keys(net.axons)
            targets = net.axons[name]
            axon_y -= axon_margin
            
            start = Point2f(0,axon_y)
            wires[name] = WireNet([
                [start,neurons[tgt.target[1][].name[]][:inputs_ports][][tgt.target[2]]] for tgt in targets
            ])
            input_ports[name] = start
        end
    end
    
    # collect axons between neurons
    row_y = axon_y
    for axon in horizontal_axons 
        idx=findfirst(row->all(other->isempty(axon.x_range ∩ other), row), row_blockers)
        if isnothing(idx)
            push!(row_blockers, IntervalSets.ClosedInterval{Float32}[])
            idx = length(row_blockers)
        end
        push!(row_blockers[idx], axon.x_range)
        axon_y = row_y-(idx+1)*axon_margin
        soma::Observable{Point2f} = neurons[axon.name][:soma]
        knee::Observable{Point2f} = @lift($(soma)+Point2f(0,axon_y))
    
        wrs = Vector{Q where Q}[[soma,knee]]
        # if there is exactly one target, connect directly
        if length(axon.tgts) == 1
            push!(wrs[1], axon.tgts[1])
        end

        if length(axon.tgts)>1 
            if axon.tgts[1][1] < soma[][1] < axon.tgts[end][1] 
                # soma is in between two targets -> draw a fork
                push!(wrs, [axon.tgts[1], knee, axon.tgts[end]])
                append!(wrs, [[Point2f(tgt[1], axon_y),tgt] for tgt in axon.tgts[2:end-1]])
            elseif soma[][1] < axon.tgts[1][1] 
                # soma is left of all targets -> draw a right arm
                push!(wrs, [knee, axon.tgts[end]])
                append!(wrs, [[Point2f(tgt[1], axon_y),tgt] for tgt in axon.tgts[1:end-1]])
            elseif soma[][1] > axon.tgts[end][1] 
                # soma is right of all targets -> draw a left arm
                push!(wrs, [axon.tgts[1],knee])
                append!(wrs, [[Point2f(tgt[1], axon_y),tgt] for tgt in axon.tgts[2:end]])
            end
        end

        wires[axon.name] = WireNet(wrs)

        if axon.name in net.outputs
            output_ports[axon.name] = Point2f(minimum(axon.x_range), axon_y)
        end
    end

    
    # draw wires
    w=wireplot!(netplot, wires; bridge_radius, csegs, linewidth)
    
    # set colors of all the inputs
    onany(netplot[:axon_colors]) do cols
        _to_update = Set{Symbol}()
        for key in keys(net.axons)
            value = cols[key]
            w[:net_colors][][key] = value
            for target in net.axons[key]
                name = Symbol("$(target.target[1][].name[])")
                push!(_to_update, name)
                netplot[name][][:inputs_color][][target.target[2]] = value
            end
        end
        for name in _to_update
            notify(netplot[name][][:inputs_color])
        end
        notify(w[:net_colors])
    end

    # set colors of all the dendrites
    on(netplot[:dendrite_colors]) do cols
        for n_key in keys(net.neurons)
            tmp = netplot[Symbol("$(n_key)")][][:dendrites_color]
            for d_key in keys(tmp[])
                value = cols[(n_key,d_key)]
                tmp[][d_key] = value
            end
            notify(tmp)
        end
    end

    # set colors of all the spines
    on(netplot[:spine_colors]) do cols
        for n_key in keys(net.neurons)
            tmp = netplot[Symbol("$(n_key)")][][:spines_color]
            for (d_key,s_key) in keys(tmp[])
                value = cols[(n_key,d_key,s_key)]
                tmp[][(d_key,s_key)] = value
            end
            notify(tmp)
        end
    end

    notify(netplot[:spine_colors])
    notify(netplot[:axon_colors])
    notify(netplot[:dendrite_colors])

    netplot[:wires] = w
    netplot[:input_ports] = input_ports
    netplot[:output_ports] = output_ports
    netplot[:tpn] = tpn
    netplot[:state] = state

    set_colors_to_state!(netplot)

    netplot
end


function set_colors_to_state!(netplot)
    state = netplot[:state]

    (;default_axon_color,default_dendrite_color,default_spine_color) = netplot.attributes
    
    on(state) do state
        net = netplot[:net][]
        net_tpn =  netplot[:tpn][]
        
        # Process axons
        for name in keys(netplot[:axon_colors][])
            tgt_name = name in net.inputs ? "input_$(name)_on" : "$(name)_$(net.neurons[name].dendrite.name)_on"
            tgt_idx = net_tpn.P_index[Symbol(tgt_name)]
            netplot[:axon_colors][][name] = 
                state.m[tgt_idx] > 0 ? default_axon_color.active[] : default_axon_color.passive[]
        end

        # Process spines
        for (neuron_name, dendrite_name, input_name) in keys(netplot[:spine_colors][])
            tgt_idx = net_tpn.P_index[Symbol("$(neuron_name)_$(dendrite_name)_$(input_name)_on")]
            netplot[:spine_colors][][(neuron_name, dendrite_name, input_name)] = 
                state.m[tgt_idx] > 0 ? default_spine_color.active[] : default_spine_color.passive[]
        end

        # Process dendrites
        for (neuron_name, dendrite_name) in keys(netplot[:dendrite_colors][])
            tgt_idx = net_tpn.P_index[Symbol("$(neuron_name)_$(dendrite_name)_on")]
            netplot[:dendrite_colors][][(neuron_name, dendrite_name)] = 
                state.m[tgt_idx] > 0 ? default_dendrite_color.active[] : default_dendrite_color.passive[]
        end

        notify(netplot[:spine_colors])
        notify(netplot[:axon_colors])
        notify(netplot[:dendrite_colors])
    end

    notify(state)
end
