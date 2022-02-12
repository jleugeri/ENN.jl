#using GLMakie
using CairoMakie
using Base.Iterators, LinearAlgebra
using ENN.Neurons, ENN.TimePetriNets
CairoMakie.activate!()
##

struct WireNet{P}
    specified_connections::Vector{Vector{Observable{P}}}
    connections::Observable{Vector{Vector{P}}}
    function WireNet(connections::Vector{Vector{Q}}) where {P,Q<:Union{P,Observable{P}}}
        obs = filter(x->isa(x,Observable), collect(flatten(connections)))
        
        new_connections = Observable([
            [zero(P) for i in 1:(length(connection)*2-1)] for connection in connections
        ])

        function update!(args...)
            for (connection,spec_connection) in zip(new_connections[],connections)
                end_horizontal = if length(spec_connection) >= 3
                    n1,n2,n3 = to_value.(spec_connection[end-2:end])
                    !(n1[1] ≤ n2[1] ≤ n3[1] || n1[1] ≥ n2[1] ≥ n3[1])
                else
                    false
                end
            
                for i in reverse(eachindex(spec_connection[1:end-1]))
                    n1,n2 = to_value.(spec_connection[i:i+1])
                    intermediate = if n1[1] == n2[1]
                        end_horizontal = true
                        (n1 + n2)/2
                    elseif n1[2] == n2[2]
                        end_horizontal = false
                        (n1 + n2)/2
                    elseif end_horizontal
                        end_horizontal = false
                        P(n1[1],n2[2])
                    else
                        end_horizontal = true
                        P(n2[1],n1[2])
                    end
                    connection[2i-1] = n1
                    connection[2i] = intermediate
                end
                connection[end] = to_value(spec_connection[end])
            end
            notify(new_connections)
        end

        onany(update!, obs...)
        update!()

        new{P}(connections, new_connections)
    end
end

##
##
function intersectSegments(A::P, B::P, C::P, D::P)::Union{Nothing, P} where P
    ΔAB = B-A
    ΔAC = C-A
    ΔAD = D-A

    # fast check of bounding box
    lu1 = min.(A,B)
    ro1 = max.(A,B)
    lu2 = min.(C,D)
    ro2 = max.(C,D)
    if any(ro2 .< lu1) || any(ro1 .< lu2)
        # bounding boxes don't intersect -> no chance for overlap
        return nothing
    end

    ΔAB′ = ΔAB/(ΔAB⋅ΔAB)

    # get first projection
    ρ1 = (ΔAB′⋅ΔAC)
    δ1 = ΔAC-ρ1 * ΔAB
    # get second projection
    ρ2 = (ΔAB′⋅ΔAD)
    δ2 = ΔAD-ρ2 * ΔAB

    # check if both δ point in the same direction 
    if δ1⋅δ2 > 0
        # they point in the same direction -> no intersection
        return nothing
    else
        # they point in opposite directions -> compute intersection
        d1 = √(δ1⋅δ1)
        d2 = √(δ2⋅δ2)
        
        d12 = d1+d2
        if d12 ≈ 0
            # CD  must lie on the line AB -> Error
            AssertionError("Line segments $((A,B)) and $((C,D)) are colinear!")
        end

        ρ3 = d1 / d12 * (ρ2-ρ1) + ρ1
        if 0.0 ≤ ρ3 ≤ 1.0
            return A + ρ3 * ΔAB
        else
            return nothing
        end
    end
end

function splitSegment(p_start,p_end,p_split,r_split)
    Δ1 = p_split-p_start
    Δ2 = p_split-p_end
    n1 = norm(Δ1)
    n2 = norm(Δ2)

    ret1 = n1 < r_split ? nothing : (p_start, p_split-r_split*Δ1/n1)
    ret2 = n2 < r_split ? nothing : (p_split-r_split*Δ2/n2, p_end)
    return ret1, ret2
end

##
    
(anglepoint(center::P, angle, radius) where {P}) = center + P(cos(angle),sin(angle))*radius

@recipe(WirePlot, nets) do scene
    Attributes(
        bridge_radius=0.1,
        linewidth=3,
        csegs=20,
    )
end

function Makie.plot!(wireplot::WirePlot)
    nets = wireplot[1][]
    net_names = keys(nets)
    all_wires = Dict((k=>Observable(Point2f[]) for k in net_names)...)
    function update_wires!(nets...)
        all_segments=[[Tuple{Point2f,Point2f}[] for connection in net] for net in nets]
        for (i,(key,net)) in enumerate(zip(net_names, nets))
            line = all_wires[i][]
            empty!(line)
            for (j,connection) in enumerate(net)
                for k in 1:length(connection)-1
                    segments = Union{Nothing, Tuple{Point2f,Point2f}}[(connection[k],connection[k+1])]
                    while !isempty(segments)
                        segment = pop!(segments)
                        split=false
                        for other_net_connections in all_segments[1:i-1]
                            for other_net_connection in other_net_connections
                                for other_segment in other_net_connection
                                    pt = intersectSegments(segment[1],segment[2],other_segment[1],other_segment[2])
                                    if !isnothing(pt) 
                                        (seg1,seg2) = splitSegment(segment[1], segment[2], pt, bridge_radius)
                                        if !isnothing(seg2)
                                            push!(segments, seg2)
                                        end
                                        if !isnothing(seg1)
                                            push!(segments, seg1)
                                        end
                                        split=true
                                        break
                                    end
                                end
                                if split
                                    break
                                end
                            end
                            if split
                                break
                            end
                        end
    
                        if !split
                            push!(all_segments[i][j], segment)
    
                            if isempty(line) || isnan(line[end])
                                append!(line,segment)
                            elseif segment[1] ≈ line[end]
                                push!(line, segment[end])
                            else
                                Δ = segment[1]-line[end]
                                slope = Δ[1] ≈ 0.0 ? sign(Δ[2])*π/2 : atan(Δ[2]/Δ[1])
                                append!(line,anglepoint.(Ref(0.5*line[end] + 0.5*segment[1]), LinRange(slope+π, slope, csegs), 0.5*norm(Δ)))
                                push!(line, segment[end])
                            end
                        end
                    end
                end
                push!(line, Point2f(NaN,NaN))
            end
            notify(all_wires[i])
        end
    end
    onany(update_wires!, (net.connections for net in nets)...) 
    notify(net1.connections)

    col = get(wireplot.attributes,:color,Observable(:black))
    wireplot[:net_color] = Dict(name=>to_value(col) for name in keys(all_wires))

    merged_wires = lift(values(all_wires)...) do wires...
        [wires...;]
    end
    merged_colors = lift(wireplot[:net_color], values(all_wires)...) do cols, wires...
        colors = []
        for (key,wire) in zip(keys(all_wires),wires)
            append!(colors, fill(cols[key], length(wire)))
        end
        colors
    end

    l=lines!(wireplot, merged_wires; wireplot.attributes..., color=merged_colors)

    wireplot
end

##
#=
f = Figure()
ax = Axis(f[1, 1])

A = Point2f(0,0)
B = Point2f(1,0)
C = Point2f(0.5,1)
D = Point2f(0.5,-0.5)
E = Point2f(0.75,0.75)
net1 = WireNet([[A,B,C],[A,C]])
net2 = WireNet([[D,E]])
nets = [net1,net2]

res=wireplot!(ax, nets)
display(f)
=#