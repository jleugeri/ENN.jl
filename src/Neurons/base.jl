export SomaOrDendrite, DendriteSegment, Neuron, Synapse, NeuralNetwork, SynapseType, excitatory, inhibitory
using DataStructures
abstract type SomaOrDendrite{T} end

mutable struct DefaultParams{T}
    children::Vector{Ref{DefaultParams{T}}}
    epsp_duration::T
    ipsp_duration::T
    plateau_duration::T
    dendritic_threshold::UInt
    synaptic_threshold::UInt
    axon_delay::T
    refractory_period::T
    p_trans::Union{Float64,Vector{Float64}}
    _manually_set::DefaultDict{Symbol,Bool}
    function DefaultParams{T}(; kwargs...) where T
        obj = new{T}(Ref{DefaultParams{T}}[], T(1),T(1),T(100), 1, 1, T(0), T(1), 1.0, DefaultDict{Symbol,Bool}(false))

        for (key,value) in kwargs
            setproperty!(obj, key, value)
            obj._manually_set[key] = true
        end
        return obj
    end
end

function inherit!(p::DefaultParams, other::DefaultParams; connect=true)
    for (key,value) in other._manually_set
        if !p._manually_set[key]
            setproperty!(p, key, getproperty(other, key))
        end
    end
    for child in p.children
        inherit!(child[], p; connect=false)    
    end
    
    if connect
        push!(other.children, Ref(p))
    end
    nothing
end

##
struct DendriteSegment{T} <: SomaOrDendrite{T}
    children::Vector{DendriteSegment{T}}
    parent::Ref{SomaOrDendrite{T}}
    neuron::Ref{SomaOrDendrite{T}}
    exc_inputs::Vector{Symbol}
    name::Symbol
    parameters::DefaultParams{T}
    function DendriteSegment{T}(children::Vector{DendriteSegment{T}}, exc_inputs::Vector{Symbol}, name::Symbol; kwargs...) where T

        this = new{T}(children, Ref{SomaOrDendrite{T}}(), Ref{SomaOrDendrite{T}}(), exc_inputs, name, DefaultParams{T}(;kwargs...))
        for child in children
            child.parent[]=this
            inherit!(child.parameters, this.parameters)
        end
        this
    end
end

function DendriteSegment{T}(children, exc_inputs::Vector{Symbol}, name::Symbol, kwargs=Dict{Symbol,Any}()) where T
    new_children = DendriteSegment{T}[DendriteSegment{T}(child...) for child in children]
    invoke(DendriteSegment{T}, Tuple{Vector{DendriteSegment{T}},Vector{Symbol},Symbol}, new_children, exc_inputs, name; kwargs...)
end
##
struct Neuron{T} <: SomaOrDendrite{T}
    all_segments::Vector{DendriteSegment{T}}
    dendrite::DendriteSegment{T}
    ports::Dict{Symbol,Vector{DendriteSegment{T}}}
    name::Ref{Symbol}
    parameters::DefaultParams{T}

    function Neuron{T}(args...; kwargs...) where T
        dendrite = DendriteSegment{T}(args...)
        params = DefaultParams{T}(;kwargs...)

        # store all segments in BFS-order
        head=0
        all_segments = [dendrite]
        while head < length(all_segments)
            head+=1
            d=all_segments[head]
            append!(all_segments, d.children)
        end

        ports = Dict{Symbol,Vector{DendriteSegment{T}}}()
        for segment in all_segments
            for input in segment.exc_inputs
                pp=get(ports, input, DendriteSegment{T}[])
                push!(pp, segment)
                ports[input]=pp
            end
            ports[Symbol("inh_$(segment.name)")] = [segment]
        end

        this = new{T}(all_segments,dendrite,ports, Ref(:neuron), params)
        
        function recursively_set_neuron!(d::DendriteSegment{T}, n::Neuron{T})
            d.neuron[] = n
            for child in d.children
                recursively_set_neuron!(child, n)
            end
        end

        dendrite.parent[]=this
        recursively_set_neuron!(dendrite, this)
        inherit!(dendrite.parameters, this.parameters)

        this
    end
end

##
struct Axon{T}
    source::Union{Ref{Neuron{T}},Symbol}
    target::Union{Tuple{Ref{Neuron{T}},Symbol},Symbol}
    parameters::DefaultParams{T}
    function Axon{T}(source, target; kwargs...) where T
        params = DefaultParams{T}(;kwargs...)
        if isa(source, Neuron)
            inherit!(params, source.parameters)
        end
        new{T}(source, target, params)
    end
end

struct NeuralNetwork{T}
    neurons::OrderedDict{Symbol,Neuron{T}}
    inputs::Vector{Symbol}
    outputs::Vector{Symbol}
    axons::OrderedDict{Symbol,Vector{Axon{T}}}
    parameters::DefaultParams{T}
    function NeuralNetwork{T}(neurons, inputs, outputs, axons; kwargs...) where T
        params = DefaultParams{T}(;kwargs...)
        
        for neuron in values(neurons)
            inherit!(neuron.parameters,params)
        end
        for axon_bundle in values(axons)
            for axon in axon_bundle
                inherit!(axon.parameters,params)
            end
        end

        new{T}(neurons, inputs, outputs, axons, params)
    end
end

function NeuralNetwork(neurons::OrderedDict{Symbol,Neuron{T}}, inputs::Vector{Symbol}, outputs::Vector{Symbol}, axons::OrderedDict{Symbol,Vector{X}}; kwargs...) where {T,X <:Union{Tuple,NTuple,NamedTuple}}
    for (name,neuron) in pairs(neurons)
        neuron.name[] = name
    end

    _axons = OrderedDict{Symbol,Vector{Axon}}()
    for (source_name, terminals) in pairs(axons)
        _terminals = getkey(_axons, source_name, Axon[])
        src = source_name ∈ inputs ? source_name : Ref(neurons[source_name])
        for (target_neuron_name, target_port_name, axon_kwargs...) in terminals
            tgt = target_neuron_name ∈ outputs ? target_neuron_name : (Ref(neurons[target_neuron_name]),target_port_name)
            syn = Axon{T}(src, tgt; axon_kwargs...)
            push!(_terminals, syn)
        end
        _axons[source_name] = _terminals
    end
    NeuralNetwork{T}(neurons, inputs, outputs, _axons; kwargs...)
end
