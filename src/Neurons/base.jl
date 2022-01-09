export SomaOrDendrite, DendriteSegment, Neuron

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
