export SomaOrDendrite, DendriteSegment, Neuron, Synapse, NeuralNetwork

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

        this = new(all_segments,dendrite)
        dendrite.parent[]=this
        this
    end
end


struct Synapse
    source::Ref{Neuron}
    target::Tuple{Ref{DendriteSegment},Symbol}
    type::Symbol
end

struct NeuralNetwork
    neurons::Dict{Symbol,Neuron}
    synapses::Dict{Symbol,Vector{Synapse}}
end

function NeuralNetwork(neurons, synapses::Dict{Symbol,Vector{X}}) where X <:Union{Tuple,NTuple,NamedTuple}
    _synapses = Dict{Symbol,Vector{Synapse}}()
    for (source_name, terminals) in pairs(synapses)
        _terminals = getkey(_synapses, source_name, Synapse[])
        src = neurons[source_name]
        for (target_neuron_name,target_segment_name,target_spine_name, attrs...) in terminals
            tgt_n = neurons[target_neuron_name]
            idx = findfirst(seg->seg.name==target_segment_name,tgt_n.all_segments)
            tgt = tgt_n.all_segments[idx]
            typ = target_spine_name==:inh ? :inh : :excitatory
            syn = Synapse(Ref(src), (Ref(tgt), target_spine_name), typ, attrs...)
            push!(_terminals, syn)
        end
        _synapses[source_name] = _terminals
    end
    NeuralNetwork(neurons,_synapses)
end
