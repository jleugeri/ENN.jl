module Neurons
using Requires

include("base.jl")
include("to_TPN.jl")
include("to_TA.jl")

function __init__()
    # Optional dependency for makie-based plotting
    @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" begin
        include("plot_network.jl")
    end
end

end