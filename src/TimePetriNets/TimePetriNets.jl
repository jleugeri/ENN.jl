module TimePetriNets

using Requires

include("base.jl")
include("Monoflop.jl")

function __init__()
    # Optional dependency for makie-based plotting
    @require GraphMakie="1ecd5474-83a3-4783-bb4f-06765db800d2" begin
        @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" begin
            include("plotting_makie.jl")
        end
    end
end

end