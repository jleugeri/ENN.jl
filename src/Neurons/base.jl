export SomaOrDendrite, DendriteSegment, Neuron, Synapse, NeuralNetwork
using DataStructures
abstract type SomaOrDendrite end

struct DendriteSegment <: SomaOrDendrite
    children::Vector{DendriteSegment}
    parent::Ref{SomaOrDendrite}
    inputs::Vector{Symbol}
    name::Symbol
    dendritic_threshold::Int
    synaptic_threshold::Int
    function DendriteSegment(children::Vector{DendriteSegment}, inputs::Vector{Symbol}, name::Symbol, dendritic_threshold::Int=isempty(children) ? 0 : 1, synaptic_threshold::Int=length(inputs))
        if dendritic_threshold > length(children)
            @warn "Threshold '$(dendritic_threshold)' is larger than number of children '$(length(children))'."
        end
        
        this = new(children,Ref{SomaOrDendrite}(), inputs, name, dendritic_threshold, synaptic_threshold)
        for child in children
            child.parent[]=this
        end
        this
    end
end

function DendriteSegment(children, inputs::Vector{Symbol}, name::Symbol, dendritic_threshold::Int=isempty(children) ? 0 : 1, synaptic_threshold::Int=length(inputs))
    new_children = Vector{DendriteSegment}()
    for child in children
        push!(new_children, DendriteSegment(child...))
    end
    invoke(DendriteSegment, Tuple{Vector{DendriteSegment},Vector{Symbol},Symbol,Int,Int}, new_children, inputs, name, dendritic_threshold, synaptic_threshold)
end

struct Neuron <: SomaOrDendrite
    all_segments::Vector{DendriteSegment}
    dendrite::DendriteSegment
    name::Ref{Symbol}
    function Neuron(args...; kwargs...)
        dendrite = DendriteSegment(args...; kwargs...)

        # store all segments in BFS-order
        head=0
        all_segments = [dendrite]
        while head < length(all_segments)
            head+=1
            d=all_segments[head]
            append!(all_segments, d.children)
        end

        this = new(all_segments,dendrite,Ref(:neuron))
        dendrite.parent[]=this
        this
    end
end


struct Synapse
    source::Union{Ref{Neuron},Nothing}
    target::Tuple{Ref{Neuron},Symbol}
    target_dendrites::Vector{Ref{DendriteSegment}}
    type::Symbol
end

struct NeuralNetwork
    neurons::OrderedDict{Symbol,Neuron}
    inputs::OrderedDict{Symbol,Vector{Synapse}}
    outputs::OrderedDict{Symbol,Symbol}
    synapses::OrderedDict{Symbol,Vector{Synapse}}
end

function NeuralNetwork(neurons, inputs::OrderedDict{Symbol,Vector{X}}, outputs::OrderedDict{Symbol,Symbol}, synapses::OrderedDict{Symbol,Vector{X}}) where X <:Union{Tuple,NTuple,NamedTuple}
    for (name,neuron) in pairs(neurons)
        neuron.name[] = name
    end

    _inputs = OrderedDict{Symbol,Vector{Synapse}}()
    for (source_name, terminals) in pairs(inputs)
        _terminals = getkey(_inputs, source_name, Synapse[])
        for (target_neuron_name,target_port_name, typ, attrs...) in terminals
            tgt_n = neurons[target_neuron_name]
            idxs = findall(seg->target_port_name ∈ seg.inputs,tgt_n.all_segments)
            syn = Synapse(nothing, (Ref(tgt_n), target_port_name), [Ref(tgt_n.all_segments[idx]) for idx in idxs], typ, attrs...)
            push!(_terminals, syn)
        end
        _inputs[source_name] = _terminals
    end

    _synapses = OrderedDict{Symbol,Vector{Synapse}}()
    for (source_name, terminals) in pairs(synapses)
        _terminals = getkey(_synapses, source_name, Synapse[])
        src = neurons[source_name]
        for (target_neuron_name,target_port_name, typ, attrs...) in terminals
            tgt_n = neurons[target_neuron_name]
            idxs = findall(seg->target_port_name ∈ seg.inputs,tgt_n.all_segments)
            syn = Synapse(Ref(src), (Ref(tgt_n), target_port_name), [Ref(tgt_n.all_segments[idx]) for idx in idxs], typ, attrs...)
            push!(_terminals, syn)
        end
        _synapses[source_name] = _terminals
    end
    NeuralNetwork(neurons, _inputs, outputs, _synapses)
end
