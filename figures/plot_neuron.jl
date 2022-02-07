using GLMakie, Base.Iterators, LinearAlgebra
using ENN.Neurons, ENN.TimePetriNets
##
struct WireNet{P}
    wires::Vector{Vector{P}}
    wire_segments::Vector{Tuple{P,P}}
end

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

#=
for w in w1.wires
    start1 = first(w)
    for end1 in drop(w,1)
        for v in w2.wires
            start2 = first(v)
            for end2 in drop(v,1)

                a1 = l1.e.y - l1.s.y
                b1 = l1.s.x - l1.e.x
                c1 = a1 * l1.s.x + b1 * l1.s.y
            
                a2 = l2.e.y - l2.s.y
                b2 = l2.s.x - l2.e.x
                c2 = a2 * l2.s.x + b2 * l2.s.y
            
                Δ = a1 * b2 - a2 * b1
                # If lines are parallel, intersection point will contain infinite values
                return Point((b2 * c1 - b1 * c2) / Δ, (a1 * c2 - a2 * c1) / Δ)
                start2=end2
            end
        end
        start1=end1
    end
end

f = Figure()
ax = Axis(f[1, 1])

for wn in [w1,w2]
    for w in wn.wires
        lines!(ax, [w[1], w[2]])
        scatter!(ax, w[3])
    end
end

f
=#
##
w1 = WireNet{Point2f}([[Point2f(0,0),Point2f(1,1)]],[])
w2 = WireNet{Point2f}([[Point2f(0,1),Point2f(1,0),Point2f(0.5,1.0)]],[])

f = Figure()
ax = Axis(f[1, 1])

wireNets = [w1,w2]

function plot_wireNets(wireNets, bridge_radius=0.1)
    wireSegments = [[[Tuple{Point2f,Point2f}[] for line in wire[1:end-1]] for wire in wireNet.wires] for wireNet in wireNets]
    for (i,wireNet) in enumerate(wireNets)
        for (j,wire) in enumerate(wireNet.wires)
            for k in 1:length(wire)-1
                segments = Union{Nothing, Tuple{Point2f,Point2f}}[(wire[k],wire[k+1])]
                while !isempty(segments)
                    segment = pop!(segments)
                    split=false
                    for wireNet2 in wireSegments[1:i-1]
                        for wire2 in wireNet2
                            for line2 in wire2
                                for segment2 in line2
                                    pt = intersectSegments(segment[1],segment[2],segment2[1],segment2[2])
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
                        if split
                            break
                        end
                    end

                    if !split
                        push!(wireSegments[i][j][k], segment)
                    end
                end
            end
        end

        l=linesegments!(ax, collect(flatten(flatten.(wireSegments[i]))))
        
        for wire in wireSegments[i]
            for line in wire
                for k in eachindex(line[1:end-1])
                    p_start = line[k][2]
                    p_end   = line[k+1][1]

                    Δ = p_end-p_start
                    slope = Δ[1] ≈ 0.0 ? sign(Δ[2])*π/2 : atan(Δ[2]/Δ[1])
                    arc!(ax, 0.5*p_start + 0.5*p_end, 0.5*norm(Δ), slope+0, slope+π, color=l[:color])
                end
            end
        end
    end
end

struct NeuronModule{P}
    ports_left::Dict{Symbol,P}
    ports_right::Dict{Symbol,P}
    ports_top::Dict{Symbol,P}
    ports_bottom::Dict{Symbol,P}

end

##

struct TreeLayout{T,F}
    rows::Vector{Vector{Vector{NamedTuple{(:node,:position),Tuple{T,Observable{Point{2,F}}}}}}}
    x_min::Vector{F}
    x_max::Vector{F}
    y_min::Vector{F}
    y_max::Vector{F}
    x_margin::Observable{F}
    y_margin::Observable{F}
end
TreeLayout(root::T; pos::Observable{Point{2,F}}=Observable(zero(Point2f)), x_margin=0.1,y_margin=0.1, height=(_)->1.0, width=(_)->1.0) where {T,F} = 
    TreeLayout{T,F}([[[(node=root,position=pos)]]],F[pos[][1]-width(root)/2],F[pos[][1]+width(root)/2],F[pos[][2]],F[pos[][2]+height(root)],x_margin,y_margin)


##
function add_right_branch!(tree::TreeLayout{T,F}, branch::TreeLayout{T,F}; keep_centered=true) where {T,F}
    # add dummy root to right branch
    pushfirst!(branch.x_min, Inf)
    pushfirst!(branch.x_max, -Inf)
    pushfirst!(branch.y_min, zero(F))
    pushfirst!(branch.y_max, zero(F))

    # get all arrays to the same size
    new_length = maximum(length.((
                                    tree.x_min, tree.x_max, tree.y_min, tree.y_max, 
                                    branch.x_min, branch.x_max, branch.y_min, branch.y_max
                                )))
    append!(tree.x_min, fill(Inf, new_length-length(tree.x_min)))
    append!(tree.x_max, fill(-Inf, new_length-length(tree.x_max)))
    append!(tree.y_min, fill(0, new_length-length(tree.y_min)))
    append!(tree.y_max, fill(0, new_length-length(tree.y_max)))
    append!(branch.x_min, fill(Inf, new_length-length(branch.x_min)))
    append!(branch.x_max, fill(-Inf, new_length-length(branch.x_max)))
    append!(branch.y_min, fill(0, new_length-length(branch.y_min)))
    append!(branch.y_max, fill(0, new_length-length(branch.y_max)))

    append!(tree.rows, [[NamedTuple{(:node,:position),Tuple{T,Observable{Point{2,F}}}}[]] for i in 1:(new_length-length(tree.rows))])

    # compute x-offset between the trees
    tree.x_margin[] = max(tree.x_margin[], branch.x_margin[])
    tree.y_margin[] = max(tree.y_margin[], branch.y_margin[])
    x_offset = max(maximum(tree.x_max .+ tree.x_margin[] .- branch.x_min),zero(F)) 
       

    # compute new "hull" of the combined tree
    tree.x_min .= min.(tree.x_min, branch.x_min.+x_offset)
    tree.x_max .= max.(tree.x_max, branch.x_max.+x_offset)

    # raise row
    min_height = max.(tree.y_max .- tree.y_min, branch.y_max .- branch.y_min)
    
    # shift branch nodes to right positions and integrate into tree
    tmp_y = tree.y_max[1] + tree.y_margin[]
    for (i,(tree_row,branch_row,row_height)) in enumerate(zip(tree.rows[2:end], branch.rows, min_height[2:end]))
        for clique in branch_row
            for node in clique
                node.position[] = Point{2,F}(node.position[][1]+x_offset, tmp_y)
            end
        end

        # on the lowest level, cliques are merged
        if i==1
            append!(tree_row[1], branch_row[1])
        else
            append!(tree_row, branch_row)
        end
        tree.y_min[i+1] = tmp_y
        tree.y_max[i+1] = tmp_y+row_height
        
        tmp_y = tree.y_max[i+1] + tree.y_margin[]
    end

    # compute new common root position
    new_root_x = (tree.x_min[2] + tree.x_max[2])/2
    root_width = tree.x_max[1] - tree.x_min[1]
    tree.x_min[1] = new_root_x - root_width/2
    tree.x_max[1] = new_root_x + root_width/2
    tree.rows[1][1][1].position[] = Point{2,F}(new_root_x,zero(F))

    #maybe re-center the entire tree
    if keep_centered
        tree.x_min .-= new_root_x
        tree.x_max .-= new_root_x
        for row in tree.rows
            for clique in row
                for node in clique
                   node.position[] -= Point{2,F}(new_root_x,zero(F))
                end
            end
        end
    end
    nothing
end
##

function layout(dendrite::DendriteSegment; params...)
    tree = TreeLayout(dendrite; params...)

    for child in dendrite.children
        child_branch = layout(child; params...)
        add_right_branch!(tree,child_branch)
    end
    return tree
end


##

neuron = Neuron((((((),[:A1,:A2],:A),),[:B1,:B2],:B),),[:C1,:C2],:C)
tree = layout(neuron.dendrite)
add_right_branch!(tree,TreeLayout(neuron.dendrite))
tree
##
neuron = Neuron((((((),[:E1,:E2,:E3,:E4],:E),((),[:E1,:E2],:F),((),[:G1,:G2],:G),((),[:H1,:H2],:H)),[:E1,:E2],:B),((((),Symbol[],:I),),Symbol[],:C),(((((((),[:K1,:K2],:K),((),[:L1,:L2],:L),((),[:A1,:A2],:M))),Symbol[],:J),),Symbol[],:D)),[:A1,:A2],:A)
tree = layout(neuron.dendrite)
##

f = Figure()
ax = Axis(f[1, 1])

for row in tree.rows
    for clique in row
        scatter!(ax, [node.position[] for node in clique])
    end
end
#hlines!(ax, tree.y_min, linestyle=:dash, color=:red)
#hlines!(ax, tree.y_max, linestyle=:dash, color=:blue)
#vlines!(ax, tree.x_min, linestyle=:dash, color=:red)
#vlines!(ax, tree.x_max, linestyle=:dash, color=:blue)

display(f)
##

f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())

(anglepoint(center::P, angle, radius) where {P}) = center + P(cos(angle),sin(angle))*radius


function excSynapse(foot::Point{2,F}, dist, radius; csegs=20) where F
    arc = anglepoint.(Ref(foot+Point{2,F}(-dist,0)), LinRange(pi/4, pi*7/4, csegs), radius)
    res = Point{2,F}[]
    push!(res, foot+Point{2,F}(0,sqrt(2)*radius/2), arc..., foot-Point{2,F}(0,sqrt(2)*radius/2))
    return res
end
    
function drawDendrite(bottom::Point{2,F}, width, height, cr; csegs=20, soma=false) where F
    cr = min(width/2-eps(F), height/2-eps(F), cr)
    
    # inner corners
    ictl = bottom .+ Point2(cr-width/2, -cr+height)
    ictr = bottom .+ Point2(-cr+width/2, -cr+height)
    icbl = bottom .+ Point2(cr-width/2, cr)
    icbr = bottom .+ Point2(-cr+width/2, cr)

    
    cstr = anglepoint.(Ref(ictr), LinRange(0, pi/2, csegs), cr)
    cstl = anglepoint.(Ref(ictl), LinRange(pi/2, pi, csegs), cr)
    
    res = Point{2,F}[]
    append!(res, cstr, cstl)
    if soma 
        csb = anglepoint.(Ref(bottom), LinRange(3pi/4, 9pi/4, 3*csegs), sqrt(2)*width/2)
        append!(res, csb)
        return res
    else
        csbl = anglepoint.(Ref(icbl), LinRange(pi, 3pi/2, csegs), cr)
        csbr = anglepoint.(Ref(icbr), LinRange(3pi/2, 2pi, csegs), cr)
        append!(res, csbl, csbr)
        return res
    end  
end


function draw_neuron!(ax, neuron; 
        route_to_bottom=true, 
        F=Float32, 
        dend_width=F(0.1), 
        syn_radius=dend_width/4, 
        dend_margin=F(0.1), 
        row_margin=F(0.1), 
        syn_margin=dend_width, 
        connector_kwargs=Dict(:linewidth=>5, :color=>:black), 
        axons_in_kwargs=Dict(:linewidth=>3, :color=>:black), 
        dendrite_kwargs=Dict(:strokewidth=>2, :color=>:silver, :strokecolor=>:black), 
        synapse_kwargs=Dict(:strokewidth=>2, :strokecolor=>:black), 
        csegs=20
    )

    hidedecorations!(ax)
    ax.aspect[]=DataAspect()


    function height(node::DendriteSegment)
        if isa(node.parent[],Neuron)
            (length(unique(node.inputs))+1)*syn_margin
        else
            (length(unique(flatten(n.inputs for n in node.parent[].children)))+1)*syn_margin
        end
    end
    width(node::DendriteSegment) = dend_width

    tree = layout(neuron.dendrite; height, width, x_margin=dend_margin, y_margin=row_margin)

    connector_lines = Point{2,F}[]
    dendrites = Dict{Symbol,Vector{Point{2,F}}}()
    spines = Dict{Tuple{Symbol,Symbol},Vector{Point{2,F}}}()
    
    
    all_inputs = unique(flatten(node.node.inputs for row in tree.rows for clique in row for node in clique))
    all_input_ys = Dict(name=>F[] for name in all_inputs)
    
    axons_in = Dict(name=>Point{2,F}[] for name in all_inputs)
    input_x = if route_to_bottom
        Dict(name=>minimum(tree.x_min) - dend_margin - syn_margin/2 * i for (i,name) in enumerate(all_inputs))
    else
        Dict(name=>minimum(tree.x_min) - dend_margin for name in all_inputs)
    end

    input_intersections = Dict(name=>F[] for name in all_inputs)
    input_num = 0
    for (i,row) in enumerate(tree.rows)

        row_inputs = unique(flatten(node.node.inputs for clique in row for node in clique))
        # sort in same order as overall
        order = sortperm(Vector{Int}(indexin(row_inputs, all_inputs)))
        permute!(row_inputs, order)
        input_y = Dict(v=>k*syn_margin for (k,v) in enumerate(row_inputs))
        input_joint = Dict(
            inp_name=>Point{2,F}(input_x[inp_name], input_y[inp_name] + tree.y_min[i] + syn_margin/2) for inp_name in row_inputs
        )


        row_input_num = input_num
        row_input_ys = Dict(name=>F[] for name in row_inputs)
        for clique in row
            if isempty(clique)
                continue
            end
        
            # target point
            pt5 = (first(clique).position[] + last(clique).position[])/2 - Point{2,F}(zero(F),tree.y_margin[])
            # joint
            pt4 = pt5 + Point{2,F}(zero(F),tree.y_margin[]/3)
            
            for (j,node) in enumerate(clique)
                # start points
                pt1 = node.position[]
                # knee points
                pt2 = pt1 - Point{2,F}(zero(F),tree.y_margin[]/3)
                # foot points
                pt3 = (pt2 + pt5)/2

                if i>1
                    # add leg
                    push!(connector_lines, pt1, pt2, pt3, pt4, Point{2,F}(NaN,NaN))
                end

                # add rounded rect + synapses
                yy = [inp=>input_y[inp] for inp  in node.node.inputs]
                sort!(yy, rev=true; by=last)
                
                for (inp_name,y) in yy
                    foot = pt1+Point{2,F}(-dend_width/2, y)

                    # check if this horizontal connector cuts through another vertical connector
                    idx=searchsortedfirst(all_inputs, inp_name)
                    row_input_num = max(idx,row_input_num)
                    if idx <= input_num
                        append!(input_intersections[inp_name],[y for l in (idx+1):input_num for y in all_input_ys[all_inputs[l]]])
                    end

                    push!(row_input_ys[inp_name], input_joint[inp_name][2])
                    
                    # add the horizontal connector
                    pushfirst!(axons_in[inp_name], input_joint[inp_name], foot+Point{2,F}(-dend_margin/3-syn_margin/2,syn_margin/2), foot+Point{2,F}(-dend_margin/3-syn_radius*sqrt(2)/2,syn_radius*sqrt(2)/2), Point{2,F}(NaN,NaN))
                    
                    spines[(node.node.name,inp_name)] = excSynapse(foot, dend_margin/3, syn_radius; csegs)
                end

                dendrites[node.node.name] = drawDendrite(pt1, dend_width, tree.y_max[i]-tree.y_min[i], F(Inf); csegs, soma=i==1)
            end

            if i>1
                # add vertical line
                push!(connector_lines, pt4, pt5, Point{2,F}(NaN,NaN))
            end
        end
        mergewith!(union, all_input_ys, row_input_ys)
        input_num = row_input_num
    end

    # prepare outputs
    all_dendrites = Dict{Symbol,Any}()
    all_axons = Dict{Symbol,Any}()
    all_spines = Dict{Tuple{Symbol,Symbol},Any}()
    all_ports = Dict{Symbol,Point{2,F}}()

    # complete input axons
    if route_to_bottom
        for inp_name in all_inputs
            arcs = [anglepoint.(Ref(Point{2,F}(input_x[inp_name], inter)), LinRange(3pi/2,pi/2,csegs), syn_margin/4) for inter in sort(unique(input_intersections[inp_name]))]
            
            port = Point{2,F}(input_x[inp_name], 0)
            pushfirst!(axons_in[inp_name], port, flatten(arcs)...)
            all_ports[inp_name] = port
        end
    end
    
    #draw clique connectors
    lines!(ax, connector_lines; connector_kwargs...)


    # draw synapses
    for (key,spine) in pairs(spines)
        p=poly!(ax, spine; synapse_kwargs...)
        all_spines[key] = p
    end

    # draw branches
    for (key,dendrite) in dendrites
        p=poly!(ax, dendrite; dendrite_kwargs...)
        all_dendrites[key] = p
    end

    #draw incoming axons
    for (key,axon_in) in pairs(axons_in)
        l = lines!(ax, axon_in; axons_in_kwargs...)
        all_axons[key] = l
    end
    return (;all_dendrites, all_axons, all_spines, all_ports)
end
res= draw_neuron!(ax, neuron)
display(f)
##

# parse neuron into row-digestible form
rows = [(offset=0.0, ports=neuron.dendrite.inputs, segments=DendriteSegment[neuron.dendrite])]
ports_by_row = []







seg_dist = 1.0
syn_dist = 0.1
i = 1
while i ≤ length(rows)    
    next_segments = DendriteSegment[]

    for segment in rows[i].segments
        append!(next_segments, segment.children)
    end

    if !isempty(next_segments)
        ports = unique(vcat((seg.inputs for seg in next_segments)...))
        max_height = length(ports)*syn_dist
        push!(rows, (offset = rows[i].offset + max_height + seg_dist, ports=ports, segments=next_segments))
    end
    i += 1
end

# compute ports row-by-row
for row in rows
    hlines!(ax, [row.offset])
    hlines!(ax, row.offset .+ (seg_dist/2) .+ (0:(length(row.ports)-1)) .* syn_dist, linestyle=:dash )
end

display(f)
##
function plot(neuron)


    #=
    \foreach \dendrites [count=\row, remember=\row as \parentrow] in #3 {    
        \foreach \parent/\connections [count=\branch] in \dendrites {    
            \xdef\numbranches{\branch};
        }

        \coordinate[above left=1cm and (\numbranches-1)*3.75mm of \prevL] (start_row_\row);

        \xdef\prevD{start_row_\row}


        \foreach \parent/\connections [count=\branch] in \dendrites {    
            \coordinate[right=7.5mm of \prevD] (start_branch_\branch);
            \coordinate[below=3mm of start_branch_\branch] (tmp);
            \xdef\prevD{start_branch_\branch}
            \xdef\prev{tmp}

            \foreach \val [count=\num] in \connections
            {
                \coordinate[above=3mm of \prev](base_\neuron_\row_\branch_\num);
                \ifthenelse{\equal{\val}{0}}{}{
                    \ifthenelse{\val>0}{
                        \draw (base_\neuron_\row_\branch_\num)++(-1mm,0.7071mm) -- ++(-0.5mm,0) arc(45:315:1mm) -- ++(0.5mm,0) -- cycle; 
                    }{
                        \fill (base_\neuron_\row_\branch_\num)++(-1mm,0.7071mm) -- ++(-0.5mm,0) arc(45:315:1mm) -- ++(0.5mm,0) -- cycle; 
                    }
                    \draw (base_\neuron_\row_\branch_\num)++(-2.207mm,0)+(135:1mm) -- ++(-1.5mm,1.5mm) 
                        coordinate(syn_\neuron_\row_\branch_\num);
                };
                \xdef\prev{base_\neuron_\row_\branch_\num}
            };
            \coordinate (tip_\neuron_\row_\branch) at (\prev);


            \ifthenelse{\row>1}{
                \draw[black, rounded corners=1mm] (start_branch_\branch)+(-1mm,-3mm) rectangle ($(tip_\neuron_\row_\branch)+(1mm,3mm)$);
            }{
                \draw[black] (start_branch_\branch)+(-1mm,-3mm) arc(117.56:423.43:2mm) [rounded corners=1mm] -- ($(tip_\neuron_\row_\branch)+(1mm,3mm)$) -- ++(-2mm,0) [rounded corners=0] -- cycle;
                \draw (start_branch_\branch)++(0mm,-6.75mm) -- ++(0,-1.5mm) coordinate(axon_\neuron);
            }
            
            \ifthenelse{\row>1}{
                \path ($(start_branch_\branch)!0.5!(tip_\neuron_\parentrow_\parent)$) coordinate(tmp1);
                \path (tmp1 -| tip_\neuron_\parentrow_\parent) coordinate(tmp2);
                \draw[black] (start_branch_\branch)+(0,-3mm) |- (tmp2) -- ($(tip_\neuron_\parentrow_\parent)+(0,3mm)$);
            }{}
        }

        \path (\prev -| start_row_\row) coordinate(stop_row_\row);
        \xdef\prevL{stop_row_\row}
    };
    \end{scope}
    =#
end



plot_wireNets([w1,w2])

display(f)
