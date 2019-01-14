module TimedAutomata
using LightGraphs, GraphPlot, Base.Iterators
export TimeSlice,TimePoint,TimeInterval,TimeLeftHalfLine,TimeRightHalfLine, ClockRegion, collect_constants, @Clk_str,
    MA,TMA,MATransition,TMATransition, to_graph, untime, language, ⊆, ==, PartialOrder, TotalOrder, OrderRelation, gplot,
    RegularExpression, RegAtom, RegUnion, RegConjunction, RegStar, RegEmpty, ∪, *, repr, substitute
#, TMARun, RunState

##########################
#    Clock Conditions    #
##########################
include("clock_constraints.jl")


##########################
#     Timed Automata     #
##########################

struct TMATransition{TMAState}
    s::TMAState
    a::Symbol
    λ::Vector{Int}
    δ::ClockRegion
end

struct TMA{TMAState}
    Σ::Vector{TMAState}
    Σ′::Vector{TMAState}
    num_clocks::Int
    E::Dict{TMAState, Vector{TMATransition{TMAState}}}
    I::Dict{TMAState, <:ClockRegion}
    s0::TMAState
end

#= 
mutable struct RunState{TMAState}
    state::TMAState
    clocks::Vector{Float64}
end

Base.:(==)(s1::RunState{TMAState}, s2::RunState{TMAState}) where {TMAState} = all(s1.state .== s2.state) && all(s1.clocks .≈ s2.clocks)
 =#

##########################
#   Nontimed Automata    #
##########################
struct MATransition{MAState}
    s::MAState
    a::Symbol
end

struct MA{MAState}
    Σ::Vector{MAState}
    Σ′::Vector{MAState}
    E::Dict{MAState, Vector{MATransition{MAState}}}
    s0::MAState
end

struct UMAState{TMAState, N, T}
    tma_state::TMAState
    clock_region::ClockRegion{N,T}
end
Base.:(==)(s1::UMAState{TMAState,N,T},s2::UMAState{TMAState,N,T}) where{TMAState,N,T} = (s1.tma_state == s2.tma_state) && (s1.clock_region == s2.clock_region)
Base.hash(s::UMAState{TMAState,N,T}) where {TMAState,N,T} = hash((hash(UMAState),hash(TMAState),hash(N),hash(T),hash(s.tma_state),hash(s.clock_region)))

include("analysis.jl")

function to_graph(ma::T) where T<:Union{MA,TMA}
    edges = Dict{Pair{Int,Int},Any}()
    states_u = collect(ma.Σ)
    l = LightGraphs.SimpleDiGraph(length(states_u))

    for (s,ee) ∈ pairs(ma.E)
        for e ∈ ee
            s_idx = findfirst(x->x==s,states_u)
            t_idx = findfirst(x->x==e.s,states_u)
            add_edge!(l,s_idx , t_idx)
            edges[s_idx=>t_idx] = e
        end
    end

    return (graph=l, edges=edges, vertices=states_u)
end

function GraphPlot.gplot(ma::MA{UMAState{TMAState,N,T}}) where {TMAState,N,T}
    g = to_graph(ma)

    tma_states = [v.tma_state for v ∈ g.vertices]
    clock_regions = [v.clock_region for v ∈ g.vertices]
    tma_states_u = unique(tma_states)
    clock_regions_u = unique(clock_regions)
    s_id = Vector{Int}(indexin(tma_states,tma_states_u))
    r_id = Vector{Int}(indexin(clock_regions,clock_regions_u))
    vertex_labels = ["S$(s)-R$(r)" for (s,r) ∈ zip(s_id, r_id)]
    edge_labels = ["$(g.edges[Pair(e)].a)" for e ∈ edges(g.graph)]
    
    gplot(
        g.graph, Vector{Float64}(r_id), Vector{Float64}(s_id);
        nodelabel=vertex_labels,
        edgelabel=edge_labels
    )
end

function GraphPlot.gplot(m::Union{TMA,MA})
    g = to_graph(m)

    states_u = unique(g.vertices)
    s_id = indexin(g.vertices,states_u)
    
    vertex_labels = ["S$(s)" for s ∈ s_id]
    edge_labels = ["$(g.edges[Pair(e)].a)" for e ∈ edges(g.graph)]
    
    gplot(
        g.graph;
        nodelabel=vertex_labels,
        edgelabel=edge_labels
    )
end



end