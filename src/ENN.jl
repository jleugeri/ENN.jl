module ENN

export TimedAutomata, TimePetriNets

include(joinpath("TimedAutomata/","TimedAutomata.jl"))
include(joinpath("TimePetriNets/","TimePetriNets.jl"))

end