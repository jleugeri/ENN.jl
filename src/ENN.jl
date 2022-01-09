module ENN

export TimedAutomata, TimePetriNets, Neurons

include(joinpath("TimedAutomata/","TimedAutomata.jl"))
include(joinpath("TimePetriNets/","TimePetriNets.jl"))
include(joinpath("Neurons/","Neurons.jl"))

end