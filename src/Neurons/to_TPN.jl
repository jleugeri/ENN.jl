import ..TimePetriNets

using SparseArrays

function TimePetriNets.TPN(n::Neuron, name=:neuron, τ_spike::H=1, τ_psp::H=τ_spike, τ_plateau::H=100; M=Int) where H
    # collect all places and transitions
    P = Symbol[]
    T = Symbol[]
    m₀= M[]
    eft=H[]
    lft=H[]

    R =  Tuple{Int,Int,M}[]
    ΔF = Tuple{Int,Int,M}[]

    function add_place(name,initial)
        push!(P, Symbol(name))
        push!(m₀, M(initial))
        length(P)
    end

    function add_transition(name,e,l)
        push!(T, Symbol(name))
        push!(eft, H(e))
        push!(lft, H(l))
        length(T)
    end

    # go through dendritic tree
    queue = Tuple{SomaOrDendrite,Union{Nothing,Int},Union{Nothing,Int}}[(n.dendrite,nothing,nothing)]
    while !isempty(queue)
        (seg,parent_counter,parent_trans) = pop!(queue)
        
        ## Add on-place and off-place
        seg_on = add_place("seg_$(name)_$(seg.name)_on", 0)
        seg_off = add_place("seg_$(name)_$(seg.name)_off", 1)

        # depending on #inputs, add synaptic counter
        syn_counter = length(seg.inputs) <= 1 ? nothing : 
            add_place("seg_$(name)_$(seg.name)_syn", 0)

        # depending on #children, add dendritic counter
        den_counter = length(seg.children) <= 1 ? nothing : 
            add_place("seg_$(name)_$(seg.name)_den", 0)

        ## Add on and off transition        
        (start_name,stop_name,reset_name,time) = if isnothing(parent_counter) && isnothing(parent_trans)
            "spike_$(name)_start", "spike_$(name)_stop", "spike_$(name)_reset", τ_spike
        else
            "seg_$(name)_$(seg.name)_start", "seg_$(name)_$(seg.name)_stop", "seg_$(name)_$(seg.name)_reset", τ_plateau
        end
        seg_start = add_transition(start_name, 0, 0)
        seg_stop = add_transition(stop_name, time, time)
        seg_reset = add_transition(reset_name, 0, 0)
        
        # if there is a parent counter, update it, too
        if !isnothing(parent_counter)
            push!(ΔF,(parent_counter,   seg_start,  1))
            push!(ΔF,(parent_counter,   seg_stop,  -1))
            push!(ΔF,(parent_counter,   seg_reset,  -1))
        end

        # if there is a parent transaction, make it sensitive to this plateau
        if !isnothing(parent_trans)
            push!(R,(seg_on,   parent_trans,  1))
        end

        # Set token changes for transitions
        push!(ΔF,(seg_off,  seg_start,  -1))
        push!(ΔF,(seg_on,   seg_start,   1))
        push!(ΔF,(seg_off,  seg_stop,    1))
        push!(ΔF,(seg_on,   seg_stop,   -1))
        push!(ΔF,(seg_off,  seg_reset,   1))
        push!(ΔF,(seg_on,   seg_reset,  -1))

        ## Add inhibitory synapse
        # Add trigger, on and off place
        inp = :inh
        inh_trigger = add_place("syn_$(name)_$(seg.name)_$(inp)_trigger", 0)
        inh_on = add_place("syn_$(name)_$(seg.name)_$(inp)_on", 0)
        inh_off = add_place("syn_$(name)_$(seg.name)_$(inp)_off", 1)
        
        # Add on and off transition
        syn_start = add_transition("seg_$(name)_$(seg.name)_$(inp)_start", 0, 0)
        syn_stop = add_transition("seg_$(name)_$(seg.name)_$(inp)_stop", τ_psp, τ_psp)
        
        # Set token changes for transitions
        push!(ΔF,(inh_trigger,  inh_start,  -1))
        push!(ΔF,(inh_off,      inh_start,  -1))
        push!(ΔF,(inh_on,       inh_start,   1))
        push!(ΔF,(inh_off,      inh_stop,    1))
        push!(ΔF,(inh_on,       inh_stop,   -1))
        
        # Make start transition sensitive to inhibitory synaptic input
        push!(R,(inh_off,seg_start,1))
        
        # Make reset transition sensitive to inhibitory synaptic input
        push!(R,(inh_on,seg_reset,1))


        
        ## Add all excitatory synapses
        for inp in seg.inputs
            # Add trigger, on and off place
            syn_trigger = add_place("syn_$(name)_$(seg.name)_$(inp)_trigger", 0)
            syn_on = add_place("syn_$(name)_$(seg.name)_$(inp)_on", 0)
            syn_off = add_place("syn_$(name)_$(seg.name)_$(inp)_off", 1)
            
            # Add on and off transition
            syn_start = add_transition("seg_$(name)_$(seg.name)_$(inp)_start", 0, 0)
            syn_stop = add_transition("seg_$(name)_$(seg.name)_$(inp)_stop", τ_psp, τ_psp)

            # Set token changes for transitions
            push!(ΔF,(syn_trigger,  syn_start,  -1))
            push!(ΔF,(syn_off,      syn_start,  -1))
            push!(ΔF,(syn_on,       syn_start,   1))
            push!(ΔF,(syn_off,      syn_stop,    1))
            push!(ΔF,(syn_on,       syn_stop,   -1))
            
            # if there is a synapse counter, update it, too
            if !isnothing(syn_counter)
                push!(ΔF,(syn_counter,  syn_start,   1))
                push!(ΔF,(syn_counter,  syn_stop,   -1))
            else
                syn_counter=syn_on
            end
        end

        # Make start transition sensitive to excitatory synaptic input
        if !isnothing(syn_counter)
            push!(R,(syn_counter,seg_start,seg.synaptic_threshold))
        end

        if !isnothing(den_counter)
            # Make start transition sensitive to dendritic input
            push!(R,(den_counter,seg_start,seg.dendritic_threshold))
            
            # Go through child branches with shared dendrite counter
            for child in seg.children
                push!(queue, (child, den_counter, nothing))
            end
        elseif length(seg.children)==1
            # Add single child branch without extra counter -> directly connect to seg_start
            push!(queue,(seg.children[1],nothing,seg_start))
        end
    end

    M = length(P)
    N = length(T)

    R = sparse(collect.(zip(R...))..., M, N)
    ΔF = sparse(collect.(zip(ΔF...))..., M, N)
    
    TimePetriNets.TPN(P,T,R,ΔF,eft,lft,m₀)
end


function TimePetriNets.TPN(net::NeuralNetwork, args...;kwargs...)
    # construct all neurons' tpns
    neuron_tpns = [TimePetriNets.TPN(neuron,name,args...;kwargs...) for (name,neuron) in pairs(net.neurons)]
    # combine synapse and neuron tpns into one tpn
    net_tpn = |(neuron_tpns...)

    # add synapses between the neurons
    for (source_name,synapses) in net.synapses
        # find the spike-start transition of the source neuron
        src_idx = findfirst(==(Symbol("spike_$(source_name)_start")), net_tpn.T)
        for synapse in synapses
            (dendrite,input) = synapse.target
            # traverse to the soma
            neuron = dendrite
            while !isa(neuron[],Neuron)
                neuron = neuron[].parent
            end
            neuron_name=findfirst(==(neuron[]),net.neurons)

            if synapse.type == :inh
                input = :inh
            end

            tgt_idx = findfirst(==(Symbol("syn_$(neuron_name)_$(dendrite[].name)_$(input)_trigger")), net_tpn.P)
            # add arc from transition to trigger place
            net_tpn.ΔF[tgt_idx, src_idx] += 1
        end
    end

    # update TPN
    TimePetriNets.update!(net_tpn)
    return net_tpn
end