mutable struct TreeLayout{T,F}
    rows::Vector{Vector{Vector{NamedTuple{(:node,:position),Tuple{T,Ref{Point{2,F}}}}}}}
    x_min::Vector{F}
    x_max::Vector{F}
    y_min::Vector{F}
    y_max::Vector{F}
    x_margin::F
    y_margin::F
end
TreeLayout(root::T; pos::Point{2,F}=zero(Point2f), x_margin=0.1,y_margin=0.1, height=(_)->1.0, width=(_)->1.0) where {T,F} = 
    TreeLayout{T,F}([[[(node=root,position=Ref(pos))]]],F[pos[1]-width(root)/2],F[pos[1]+width(root)/2],F[pos[2]],F[pos[2]+height(root)],x_margin,y_margin)

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

    append!(tree.rows, [[NamedTuple{(:node,:position),Tuple{T,Point{2,F}}}[]] for i in 1:(new_length-length(tree.rows))])

    # compute x-offset between the trees
    tree.x_margin = max(tree.x_margin, branch.x_margin)
    tree.y_margin = max(tree.y_margin, branch.y_margin)
    x_offset = max(maximum(tree.x_max .+ tree.x_margin .- branch.x_min),zero(F)) 
       

    # compute new "hull" of the combined tree
    tree.x_min .= min.(tree.x_min, branch.x_min.+x_offset)
    tree.x_max .= max.(tree.x_max, branch.x_max.+x_offset)

    # raise row
    min_height = max.(tree.y_max .- tree.y_min, branch.y_max .- branch.y_min)
    
    # shift branch nodes to right positions and integrate into tree
    tmp_y = tree.y_max[1] + tree.y_margin
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
        
        tmp_y = tree.y_max[i+1] + tree.y_margin
    end

    # compute new common root position
    new_root_x = (tree.x_min[2] + tree.x_max[2])/2
    root_width = tree.x_max[1] - tree.x_min[1]
    tree.x_min[1] = new_root_x - root_width/2
    tree.x_max[1] = new_root_x + root_width/2
    tree.rows[1][1][1].position[] = Point{2,F}(new_root_x,tree.rows[1][1][1].position[][2])

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

function layout(dendrite::DendriteSegment; pos=zero(Point2f), params...)
    tree = TreeLayout(dendrite;pos, params...)

    for child in dendrite.children
        child_branch = layout(child; params...)
        add_right_branch!(tree,child_branch)
    end
    return tree
end

(anglepoint(center::P, angle, radius) where {P}) = center + P(cos(angle),sin(angle))*radius

function excSynapse(foot::Point{2,F}, dist, radius; csegs=20, dir=:left) where F
    res = Point{2,F}[]
    if dir == :left
        arc = anglepoint.(Ref(foot+Point{2,F}(-dist,0)), LinRange(pi/4, pi*7/4, csegs), radius)
        push!(res, foot+Point{2,F}(0,sqrt(2)*radius/2), arc..., foot-Point{2,F}(0,sqrt(2)*radius/2))
    else
        arc = anglepoint.(Ref(foot+Point{2,F}( dist,0)), LinRange(3pi/4, -3pi/4, csegs), radius)
        push!(res, foot+Point{2,F}(0,sqrt(2)*radius/2), arc..., foot-Point{2,F}(0,sqrt(2)*radius/2))
    end
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

@recipe(NeuronPlot, neuron, offset) do scene
    Attributes(
        portside=:both,
        route_to_bottom=true, 
        dend_width=0.1f0, 
        syn_radius=0.025f0, 
        dend_margin=0.1f0, 
        row_margin=0.1f0, 
        syn_margin=0.1f0, 
        syn_stem=0.03f0,
        connector_kwargs=Dict(:linewidth=>5, :color=>:black), 
        inputs_kwargs=Dict(:linewidth=>3, :color=>:black), 
        dendrite_kwargs=Dict(:strokewidth=>2, :color=>:silver, :strokecolor=>:black), 
        spine_kwargs=Dict(:strokewidth=>2, :strokecolor=>:black), 
        csegs=20
    )
end

function Makie.plot!(neuronplot::NeuronPlot)
    F = Float32
    neuron = neuronplot[1]
    offset = neuronplot[2]

    (;
        portside, route_to_bottom, dend_width, syn_radius, dend_margin, 
        row_margin, syn_margin, syn_stem, connector_kwargs, inputs_kwargs, 
        dendrite_kwargs, spine_kwargs, csegs
    ) = neuronplot.attributes

    
    # prepare outputs
    all_ports = Dict{Symbol,Point{2,F}}()

    #hidedecorations!(ax)
    #neuronplot.aspect[]=DataAspect()

    connector_lines = Observable(Point{2,F}[])
    dendrites = Observable(Dict{Symbol,Vector{Point{2,F}}}())
    spines = Observable(Dict{Tuple{Symbol,Symbol},Vector{Point{2,F}}}())
    inputs = Observable(Dict{Symbol,Vector{Point{2,F}}}())
    y_extent = Observable(zeros(F,2))
    x_extent = Observable(zeros(F,2))

    onany(neuron,
        portside, route_to_bottom, dend_width, syn_radius, dend_margin, 
        row_margin, syn_margin, syn_stem, csegs
    ) do neuron,
        portside, route_to_bottom, dend_width, syn_radius, dend_margin, 
        row_margin, syn_margin, syn_stem, csegs

        empty!(connector_lines[])
        empty!(dendrites[])
        empty!(spines[])
        empty!(inputs[])

        # do initial layout to figure out horizontal placement
        tree = layout(neuron.dendrite; pos=Point{2,F}(zero(F),√(2)*dend_width/2), height=_->one(F), width=_->dend_width, x_margin=dend_margin, y_margin=row_margin)

        all_inputs = unique(flatten(node.node.exc_inputs for row in tree.rows for clique in row for node in clique))
        all_input_ys = Dict(name=>F[] for name in all_inputs)
        merge!(inputs[], Dict(name=>Point{2,F}[] for name in all_inputs))
        
        # figure out on which side to draw the symbol wires
        input_side = if portside == :left
            Dict(name => :left for name in all_inputs)
        elseif portside == :right
            Dict(name => :right for name in all_inputs)
        elseif portside == :both
            res = Dict{Symbol,Symbol}()
            # compute center of mass of synapse x-positions for each input symbol
            center_of_mass = Dict{Symbol,F}()
            for row in tree.rows
                for clique in row
                    for node in clique
                        for inp in node.node.exc_inputs
                            center_of_mass[inp] = get(center_of_mass, inp, zero(F))+node.position[][1]
                        end
                    end
                end
            end

            # assign side based on center of mass
            for name in all_inputs
                res[name] =  center_of_mass[name] <= 0 ? :left : :right
            end
            res
        else
            ArgumentError("'portside' must be either ':left', ':right' or ':both', not '$(portside)'.")
        end  

        isleft(name) = input_side[name]==:left
        all_inputs_left = filter(isleft, all_inputs)
        all_inputs_right = filter(!isleft, all_inputs)

        # figure out which inputs can coexist on the left and right side
        row_inputs_left_padded = Vector{Union{Nothing,Symbol}}[]
        row_inputs_right_padded = Vector{Union{Nothing,Symbol}}[]
        for row in tree.rows
            from_left = Dict{Symbol,F}()
            from_right = Dict{Symbol,F}()
            # find rightmost and leftmost extent of left and right inputs, respectively
            for clique in row
                for node in clique
                    node_x = node.position[][1]
                    for input in node.node.exc_inputs
                        if isleft(input)
                            from_left[input] = max(node_x,get(from_left,input,typemin(F)))
                        else
                            from_right[input] = min(node_x,get(from_right,input,typemax(F)))
                        end
                    end
                end
            end

            # sort left and right inputs to conform with the overall order
            row_inputs = unique(flatten(node.node.exc_inputs for clique in row for node in clique))
            row_inputs_left = filter(isleft, row_inputs)
            row_inputs_right = filter(!isleft, row_inputs)
            order_left = sortperm(collect(indexin(row_inputs_left, all_inputs_left)))
            order_right = sortperm(collect(indexin(row_inputs_right, all_inputs_right)))
            permute!(row_inputs_left, order_left)
            permute!(row_inputs_right, order_right)

            # greedy algorithm: alternatingly place left and right if possible
            left=true
            order = NTuple{2,Union{Nothing,Symbol}}[]
            while length(row_inputs_left) + length(row_inputs_right) > 0
                next_left = isempty(row_inputs_left) ? nothing  : first(row_inputs_left)
                next_right = isempty(row_inputs_right) ? nothing  : first(row_inputs_right)

                if !isnothing(next_left) && !isnothing(next_right) && from_left[next_left] + syn_radius + syn_stem ≤ from_right[next_right]
                    # both inputs fit
                    push!(order,(next_left,next_right))
                    popfirst!(row_inputs_left)
                    popfirst!(row_inputs_right)
                elseif left && !isnothing(next_left)
                    push!(order,(next_left,nothing))
                    popfirst!(row_inputs_left)
                else
                    push!(order,(nothing,next_right))
                    popfirst!(row_inputs_right)
                end
                left != left
            end
            rl,rr = zip(order...)
            push!(row_inputs_left_padded,collect(rl))
            push!(row_inputs_right_padded,collect(rr))
        end

        # store for each node how tall the corresponding row is
        node_height = Dict(
            node.node.name => (length(row_inputs_left_padded[i])+1)*syn_margin
            for (i,row) in enumerate(tree.rows) for clique in row for node in clique
        )
        
        # re-layout with correct assigment of synapses to left and right side
        tree = layout(neuron.dendrite; pos=Point{2,F}(zero(F),√(2)*dend_width/2), height=node->node_height[node.name], width=x->dend_width, x_margin=dend_margin, y_margin=row_margin)

        # if we route the signals all the way to the bottom, they should be staggered
        if route_to_bottom
            input_x = Dict(
                (name=>minimum(tree.x_min) - dend_margin - syn_margin/2 * i for (i,name) in enumerate(all_inputs_left))...,
                (name=>maximum(tree.x_max) + dend_margin + syn_margin/2 * i for (i,name) in enumerate(all_inputs_right))...
            )
        else
            input_x = Dict(
                (name=>minimum(tree.x_min) - dend_margin for name in all_inputs_left)...,
                (name=>maximum(tree.x_max) + dend_margin for name in all_inputs_right)...
            )
        end

        x_extent[] .= extrema([collect(values(input_x)); tree.x_min; tree.x_max; zero(F)])
        y_extent[] .= (zero(F), tree.y_max[end])
        notify(x_extent)
        notify(y_extent)

        # compute, for each input, which other inputs (on the same side) it's vertical connector intersects with
        input_intersections = Dict(name=>F[] for name in all_inputs)
        
        # iterate all rows
        for (i,row) in enumerate(tree.rows)
            # compute y-position for all inputs
            input_y = Dict(
                (v=>k*syn_margin + tree.y_min[i] + syn_margin/2 for (k,v) in enumerate(row_inputs_left_padded[i]) if !isnothing(v))...,
                (v=>k*syn_margin + tree.y_min[i] + syn_margin/2 for (k,v) in enumerate(row_inputs_right_padded[i]) if !isnothing(v))...
            )
            mergewith!(union, all_input_ys, input_y)

            # iterate all cliques in the row
            for clique in row
                if isempty(clique)
                    continue
                end
            
                # target point
                pt5 = (first(clique).position[] + last(clique).position[])/2 - Point{2,F}(zero(F),tree.y_margin[])
                # joint
                pt4 = pt5 + Point{2,F}(zero(F),tree.y_margin[]/3)
                
                # iterate all nodes in the clique
                for (j,node) in enumerate(clique)
                    # start points
                    pt1 = node.position[]
                    node_x = pt1[1]

                    # knee points
                    pt2 = pt1 - Point{2,F}(zero(F),tree.y_margin[]/3)
                    # foot points
                    pt3 = (pt2 + pt5)/2

                    # add leg (if not root of the tree)
                    if i>1
                        push!(connector_lines[], pt1, pt2, pt3, pt4, Point{2,F}(NaN,NaN))
                    end

                    # sort synapses by y-position
                    order = sortperm([input_y[inp] for inp in node.node.exc_inputs])
                    node_inputs = node.node.exc_inputs[order]
                    
                    # draw synapses
                    for inp_name in node_inputs
                        foot = Point{2,F}(node_x + (isleft(inp_name) ? -dend_width/2 : dend_width/2), input_y[inp_name]- syn_margin/2)
                        
                        # add the horizontal connector
                        pushfirst!(inputs[][inp_name], 
                            Point{2,F}(input_x[inp_name],input_y[inp_name]), 
                            foot+Point{2,F}((isleft(inp_name) ? -1 : 1)*(syn_stem+syn_margin/2), syn_margin/2), 
                            foot+Point{2,F}((isleft(inp_name) ? -1 : 1)*(syn_stem+syn_radius*sqrt(2)/2), syn_radius*sqrt(2)/2), 
                            Point{2,F}(NaN,NaN)
                        )
                        
                        spines[][(node.node.name,inp_name)] = excSynapse(foot, syn_stem, syn_radius; csegs, dir=input_side[inp_name])
                    end

                    dendrites[][node.node.name] = drawDendrite(pt1, dend_width, tree.y_max[i]-tree.y_min[i], F(Inf); csegs, soma=i==1)
                end

                if i>1
                    # add vertical line
                    push!(connector_lines[], pt4, pt5, Point{2,F}(NaN,NaN))
                end
            end
        end

        # complete input axons
        if route_to_bottom
            potential_intersections = F[]
            for inp_name in reverse(all_inputs_left)
                input_intersections = sort(unique(filter(<(maximum(all_input_ys[inp_name])),potential_intersections)))
                arcs = [anglepoint.(Ref(Point{2,F}(input_x[inp_name], inter)), LinRange(3pi/2,pi/2,csegs), syn_margin/4) for inter in input_intersections]
                port = Point{2,F}(input_x[inp_name], 0)
                pushfirst!(inputs[][inp_name], port, flatten(arcs)...)
                all_ports[inp_name] = port

                append!(potential_intersections, all_input_ys[inp_name])
            end
            potential_intersections = F[]
            for inp_name in reverse(all_inputs_right)
                input_intersections = sort(unique(filter(<(maximum(all_input_ys[inp_name])),potential_intersections)))
                arcs = [anglepoint.(Ref(Point{2,F}(input_x[inp_name], inter)), LinRange(3pi/2,pi/2,csegs), syn_margin/4) for inter in input_intersections]
                port = Point{2,F}(input_x[inp_name], 0)
                pushfirst!(inputs[][inp_name], port, flatten(arcs)...)
                all_ports[inp_name] = port

                append!(potential_intersections, all_input_ys[inp_name])
            end
        end
    end

    
    # force update
    neuron[] = neuron[]

    neuronplot[:x_extent] = x_extent
    neuronplot[:y_extent] = y_extent

    #draw clique connectors
    lines!(neuronplot, @lift($(connector_lines).+Ref($(offset))); connector_kwargs[]...)

    # draw spines
    col = pop!(spine_kwargs[], :color, :gray)
    neuronplot[:spines_color] = @lift(Dict(name=>col for name in keys($(spines))))
    neuronplot[:spines]=poly!(neuronplot, @lift([val.+ Ref($(offset)) for val in values($(spines))]); spine_kwargs[]...,
        color=@lift([$(neuronplot[:spines_color])[key] for key in keys($(spines))]), spine_kwargs[]...)
    neuronplot[:spine_ids]=@lift(keys($(spines)))

    # draw branches
    col = pop!(dendrite_kwargs[], :color, :gray)
    neuronplot[:dendrites_color] = @lift(Dict(name=>col for name in keys($(dendrites))))
    neuronplot[:dendrites]=poly!(neuronplot, @lift([val.+ Ref($(offset)) for val in values($(dendrites))]); 
        color=@lift([$(neuronplot[:dendrites_color])[key] for key in keys($(dendrites))]), dendrite_kwargs[]...)
    neuronplot[:dendrite_ids]=@lift(keys($(dendrites)))
    
    #draw incoming axons
    col = pop!(inputs_kwargs[], :color, :gray)
    neuronplot[:inputs_color] = @lift(Dict(name=>col for name in keys($(inputs))))
    neuronplot[:axons]=lines!(neuronplot, @lift(Point{2,F}[([val.+ Ref($(offset)); [Point{2,F}(NaN,NaN)]] for val in values($(inputs)))...;]); 
        color=@lift([(fill($(neuronplot[:inputs_color])[key], length(val)+1) for (key,val) in pairs($(inputs)))...;]), inputs_kwargs[]...)
    neuronplot[:inputs_ids]=@lift(collect(zip(keys($(inputs)),length.(values($(inputs))))))
    neuronplot[:inputs_ports] = @lift(Dict(key=>first(value)+$(offset) for (key, value) in pairs($(inputs))))

    neuronplot[:soma] = offset

    return neuronplot
end