import ..TimePetriNets
using SparseArrays

function TimePetriNets.TPN(n::Neuron{H}, name=:neuron; M=Int) where {H}
    # collect all places and transitions to glue dendrite segments + synapses together
    P = Symbol[]
    T = Symbol[]
    m₀ = M[]
    eft = H[]
    lft = H[]
    trans_probs = Float64[]

    R = Tuple{Int,Int,M}[]
    ΔF = Tuple{Int,Int,M}[]

    function add_place!(name, initial)
        push!(P, Symbol(name))
        push!(m₀, M(initial))
        length(P)
    end

    function add_transition!(name, e, l=e, p=1.0)
        push!(T, Symbol(name))
        push!(eft, H(e))
        push!(lft, H(l))
        push!(trans_probs, p)
        length(T)
    end

    function copy_transitions!(other, names...)
        idxs = indexin(Symbol.(names), other.T)
        return add_transition!.(
            names,
            other.eft[idxs],
            other.lft[idxs],
            other.trans_probs[idxs]
        )
    end

    function copy_places!(other, names...)
        idxs = indexin(Symbol.(names), other.P)
        return add_place!.(
            names,
            other.x₀.m[idxs]
        )
    end

    sub_tpns = TPN{M,H}[]

    # go through dendritic tree
    queue = Tuple{SomaOrDendrite,Union{Nothing,Int},Union{Nothing,Int}}[(n.dendrite, nothing, nothing)]
    while !isempty(queue)
        (seg, parent_counter, parent_trans) = pop!(queue)
        has_reset = !isempty(seg.inh_inputs)

        # the dendrite is mainly a monoflop
        dendrite_tpn = TimePetriNets.Monoflop(Symbol("$(name)_$(seg.name)"), 
            seg==n.dendrite ? seg.parameters.spike_duration+seg.parameters.refractory_period : seg.parameters.plateau_duration; add_reset=has_reset)
        push!(sub_tpns, dendrite_tpn)

        # need to access start, stop transition as well as the on place --> copy them
        t_names = ("$(name)_$(seg.name)_start", "$(name)_$(seg.name)_stop")
        seg_start, seg_stop = copy_transitions!(dendrite_tpn, t_names...)
        seg_reset = has_reset ? copy_transitions!(dendrite_tpn, "$(name)_$(seg.name)_reset")[1] : nothing
        (seg_on,) = copy_places!(dendrite_tpn, "$(name)_$(seg.name)_on")

        # depending on excitatory #inputs, add synaptic counter
        syn_counter_exc = length(seg.exc_inputs) == 1 ? nothing :
                          add_place!("$(name)_$(seg.name)_syn_counter_exc", 0)

        # depending on inhibitory #inputs, add synaptic counter and inverted counter
        syn_counter_inh = length(seg.inh_inputs) == 1 ? nothing :
                          add_place!("$(name)_$(seg.name)_syn_counter_inh", 0)
        syn_counter_inh_off = add_place!("$(name)_$(seg.name)_syn_counter_inh_off", length(seg.inh_inputs))

        # depending on #children, add dendritic counter
        den_counter = length(seg.children) == 1 ? nothing :
                      add_place!("$(name)_$(seg.name)_den_counter", 0)

        # if there is a parent counter, update it, too
        if !isnothing(parent_counter)
            push!(ΔF, (parent_counter, seg_start, 1), (parent_counter, seg_stop, -1))
        end

        # if there is a parent transaction, make it sensitive to this plateau
        if !isnothing(parent_trans)
            push!(R, (seg_on, parent_trans, 1))
        end

        ## Add all inhibitory synapses
        for inp in seg.inh_inputs
            # the inhibitory synapse is mainly a monoflop
            synapse_tpn = TimePetriNets.Monoflop(Symbol("$(name)_$(seg.name)_$(inp)"), seg.parameters.ipsp_duration; add_reset=false)
            push!(sub_tpns, synapse_tpn)
            # we need to access the start and stop transition and on place --> copy them 
            (t_start, t_stop) = copy_transitions!(synapse_tpn, "$(name)_$(seg.name)_$(inp)_start", "$(name)_$(seg.name)_$(inp)_stop")
            (p_on,) = copy_places!(synapse_tpn, "$(name)_$(seg.name)_$(inp)_on")
        
            # Add trigger place
            t_trigger = add_place!("$(name)_$(seg.name)_$(inp)_trigger", 0)
            push!(ΔF, (t_trigger, t_start, -1))
        
            # if there is a synapse counter, update it, too
            if !isnothing(syn_counter_inh)
                push!(ΔF, (syn_counter_inh, t_start, 1), (syn_counter_inh, t_stop, -1))
            else
                syn_counter_inh = p_on
            end
            # off-counter works in opposite direction
            push!(ΔF, (syn_counter_inh_off, t_start, -1), (syn_counter_inh_off, t_stop, 1))
        end

        if has_reset
            if !isnothing(parent_counter)
                push!(ΔF, (parent_counter, seg_reset, -1))
            end
            
            # enable reset when at least one inhibitory synapse is on
            # enable segment if all inhibitory synapses are off
            push!(R, 
                (syn_counter_inh, seg_reset, 1),
                (syn_counter_inh_off, seg_start, length(seg.inh_inputs))
            )
        end

        ## Add all excitatory synapses
        for inp in seg.exc_inputs
            # the excitatory synapse is mainly a monoflop
            synapse_tpn = TimePetriNets.Monoflop(Symbol("$(name)_$(seg.name)_$(inp)"), seg.parameters.epsp_duration; add_reset=false)
            push!(sub_tpns, synapse_tpn)
            # we need to access the start and stop transition and on place --> copy them 
            (t_start, t_stop) = copy_transitions!(synapse_tpn, "$(name)_$(seg.name)_$(inp)_start", "$(name)_$(seg.name)_$(inp)_stop")
            (p_on,) = copy_places!(synapse_tpn, "$(name)_$(seg.name)_$(inp)_on")

            # Add trigger place
            t_trigger = add_place!("$(name)_$(seg.name)_$(inp)_trigger", 0)
            push!(ΔF, (t_trigger, t_start, -1))

            # if there is a synapse counter, update it, too
            if !isnothing(syn_counter_exc)
                push!(ΔF, (syn_counter_exc, t_start, 1), (syn_counter_exc, t_stop, -1))
            else
                syn_counter_exc = p_on
            end
        end

        # Make start transition sensitive to excitatory synaptic input
        push!(R, (syn_counter_exc, seg_start, seg.parameters.synaptic_threshold))

        if !isnothing(den_counter)
            # Make start transition sensitive to dendritic input
            push!(R, (den_counter, seg_start, min(length(seg.children), seg.parameters.dendritic_threshold)))

            # Go through child branches with shared dendrite counter
            for child in seg.children
                push!(queue, (child, den_counter, nothing))
            end
        else
            # Add single child branch without extra counter -> directly connect to seg_start
            push!(queue, (seg.children[1], nothing, seg_start))
        end
    end

    m = length(P)
    n = length(T)

    R = isempty(R) ? spzeros(M,m,n) : sparse(collect.(zip(R...))..., m, n)
    ΔF = isempty(ΔF) ? spzeros(M,m,n) : sparse(collect.(zip(ΔF...))..., m, n)

    glue = TimePetriNets.TPN(P, T, R, ΔF, eft, lft, trans_probs, m₀)
    |(glue, sub_tpns...)
end


function TimePetriNets.TPN(net::NeuralNetwork, args...;kwargs...)
    # construct all neurons' tpns
    neuron_tpns = [TimePetriNets.TPN(neuron,name,args...;kwargs...) for (name,neuron) in pairs(net.neurons)]
    (H,M) = eltype(neuron_tpns).parameters

    input_trans = Monoflop.(Symbol.(Ref("input_").*String.(net.inputs)), 0)

    # combine synapse and neuron tpns into one tpn
    net_tpn = |(input_trans..., neuron_tpns...)

    # set lft of all inputs to infinity
    for input in net.inputs
        idx = net_tpn.T_index[Symbol("input_$(input)_start")]        
        net_tpn.eft[idx] = zero(H)
        net_tpn.lft[idx] = typemax(H)
    end


    function add_axon!(src_idx, axons)
        for axon in axons
            if isa(axon.target,Symbol)
                continue
            end

            (neuron,input) = axon.target
            neuron_name = neuron[].name[]

            for dendrite in neuron[].ports[input]
                tgt_idx = net_tpn.P_index[Symbol("$(neuron_name)_$(dendrite.name)_$(input)_trigger")]
                # add arc from transition to trigger place
                net_tpn.ΔF[tgt_idx, src_idx] += 1
                
                t_idx = get(net_tpn.T_index, Symbol("$(neuron_name)_$(dendrite.name)_$(input)_start"), nothing)
                if !isnothing(t_idx)
                    net_tpn.eft[t_idx] = axon.parameters.axon_delay
                    net_tpn.lft[t_idx] = axon.parameters.axon_delay
                end
                t_idx = get(net_tpn.T_index, Symbol("$(neuron_name)_$(dendrite.name)_$(input)_fail"), nothing)
                if !isnothing(t_idx)
                    net_tpn.eft[t_idx] = axon.parameters.axon_delay
                    net_tpn.lft[t_idx] = axon.parameters.axon_delay
                end
            end
        end
        nothing
    end

    # add synapses between the neurons
    for (source_name, axons) in pairs(net.axons)
        # find the spike-start transition of the source neuron or input
        source_name = if source_name ∈ net.inputs
            Symbol("input_$(source_name)")
        else
            Symbol("$(source_name)_$(net.neurons[source_name].dendrite.name)")
        end
        src_name = Symbol("$(source_name)_start")
        src_idx = net_tpn.T_index[src_name]
        add_axon!(src_idx, axons)
    end
    # update TPN
    TimePetriNets.update!(net_tpn)
    return net_tpn
end