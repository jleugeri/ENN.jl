using CairoMakie
using Base.Iterators, LinearAlgebra
using ENN.Neurons, ENN.TimePetriNets
CairoMakie.activate!()

include("plot_neuron.jl")
include("plot_wires.jl")
##


neuron1 = Neuron((((((),[:E1,:E2,:E3,:E4],:E),((),[:E1,:E2],:F),((),[:G1,:G2],:G),((),[:H1,:H2],:H)),[:E1,:E2],:B),((((),Symbol[:E1],:I),),Symbol[],:C),(((((((),[:K1,:K2],:K),((),[:L1,:L2],:L),((),[:A1,:A2],:M))),Symbol[],:J),),Symbol[:H1],:D)),[:A1,:A2],:A)
neuron2 = Neuron((((((),[:E1,:E2,:E3,:E4],:E),((),[:E1,:E2],:F),((),[:G1,:G2],:G),((),[:H1,:H2],:H)),[:E1,:E2],:B),((((),Symbol[:E1],:I),),Symbol[],:C),(((((((),[:K1,:K2],:K),((),[:L1,:L2],:L),((),[:A1,:A2],:M))),Symbol[],:J),),Symbol[:H1],:D)),[:A1,:A2],:A)


net = NeuralNetwork(
    Dict(:n1=>neuron1,:n2=>neuron2),    
    Dict(:n1=>[(:n2, :E1, :excitatory)],:n2=>[(:n1, :E1, :excitatory)])
)
##

f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())

neurons = Dict{Symbol,Combined{neuronplot, Tuple{Neuron, Point{2, Float32}}}}()
widths = Dict{Symbol,Observable{Vector{Float32}}}()
for (name,neuron) in pairs(net.neurons)
    n = neuronplot!(ax, neuron, Point2f(0.0,0.0))
    neurons[name]=n
    widths[name]=n[:x_extent]
end

neuron_margin = 0.1

onany(values(widths)...) do _widths...
    offset = neuron_margin
    for (name,(x_min,x_max)) in zip(keys(widths),_widths)
        offset -= x_min
        neurons[name][:offset][] = Point2f(offset, 0.0)
        offset += x_max + neuron_margin
        notify(neurons[name][:offset])
    end

    xlims!(ax, [0,offset])
end

notify(first(values(widths)))
##

horizontal_axons = NamedTuple{(:neuron,:x_min,:x_max),Tuple{Symbol,Float32,Float32}}[]
for (name,neuron) in pairs(neurons)
    targets = get(net.synapses, name, [])
    x_min = Inf
    x_max = -Inf
    ports = [neurons[target.target[1][].name[]][:axon_in_ports][][target.target[2]] for target in targets]
    println(ports)
end

# simple heuristic sorting


neurons[1][:soma]

#res= netplot!(ax, neuron, Point2f(0,0); portside=:both)
display(f)