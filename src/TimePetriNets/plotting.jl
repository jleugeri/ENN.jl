using GraphMakie, GLMakie, Graphs, Printf

"""
GraphMakie.graphplot(pn::TPN, x=pn.x₀, args...; nlabels_align=(:center,:bottom), nlabels_distance=25, elabels_distance=5, kwargs...)

Creates an interactive plot of the TPN `pn` in the initial state `x`. Additional `args...` and `kwargs...` are directly passed on to `graphplot(::SimpleDiGraph)`.
For a description of keyword arguments, see the `GraphMakie.jl` documentation.
"""
function GraphMakie.graphplot(pn::TPN, x=pn.x₀, args...; nlabels_align=(:center,:bottom), nlabels_distance=25, elabels_distance=5, kwargs...)   
g = SimpleDiGraph(length(pn.P)+length(pn.T))

nlabels = ([
      string.(pn.P)
    ; string.(pn.T)])
node_marker = [fill(:circle,length(pn.P));fill(:rect,length(pn.T))]

elabels=Dict{Graphs.SimpleGraphs.SimpleEdge,String}()
rows = rowvals(pn.ΔF)
vals = nonzeros(pn.ΔF)
for t_id in eachindex(pn.T)
    t_id′ = t_id + length(pn.P)
    for i in nzrange(pn.ΔF, t_id)
        p_id = rows[i]
        w = vals[i]
        e = if w > 0
            Graphs.SimpleGraphs.SimpleEdge(t_id′,p_id)
        elseif w < 0
            Graphs.SimpleGraphs.SimpleEdge(p_id,t_id′)
        else
            error("Found edge with 0 token change for id $(p_id)")
        end

        add_edge!(g,e)
        elabels[e] = w==1 ? "" : repr(Int(w))
    end
end

# draw plot
f,ax,p= graphplot(g, args...; nlabels, elabels=getindex.(Ref(elabels),edges(g)), nlabels_align, nlabels_distance, elabels_distance, node_marker, node_size=0.5, node_attr=(markerspace = SceneSpace, strokecolor=:black, strokewidth=1, glowwidth=5, glowcolor=fill(:transparent, length(pn.T)+length(pn.P))), node_color=fill(:gray,length(pn.T)+length(pn.P)), kwargs...)
ax.autolimitaspect = 1

# add important attributes to plot
tokens = Observable(Point2f[])
p.attributes.tokens=tokens
p.attributes.p_marking=Observable(copy(x.m))
p.attributes.t_marking=Observable(copy(x.h))
p.attributes.max_delay=lift((m,h)->max_delay(pn, TPNState(m,h)), p.attributes.p_marking,p.attributes.t_marking)
p.attributes.selected=1

# draw menu
f[1,2] = sidebar = GridLayout(width=Fixed(200), tellheight=false)
sidebar[1,1] =  Label(f, "Selected: ")
l_idx = sidebar[1,2] =  Label(f, lift(
    idx-> (idx>length(pn.P) ? (string(pn.T[idx-length(pn.P)])*"\n(transition)") : (string(pn.P[idx]))*"\n(place)"),
    p.attributes.selected)
)

place_menu = GridLayout()
transition_menu = GridLayout()
delay_menu = GridLayout()
save_menu = GridLayout()
sidebar[2,1:2] = place_menu
sidebar[3,1:2] = transition_menu
sidebar[4,1:2] = delay_menu
sidebar[5,1:2] = save_menu

# save menu
save_menu[1,1:2] = Label(f, "Save frame:")
t_fn = save_menu[2,1] = Textbox(f, placeholder = "filename")
b_save = save_menu[2,2] = Button(f, label="Save")

# place menu
l_nt_t= place_menu[1,1] = Label(f, "#tokens: ")
l_nt = place_menu[1,2] =  Label(f, lift((idx,m)-> (idx>length(pn.P) ? "N.A." : string(Int(m[idx]))),p.attributes.selected,p.attributes.p_marking))
b_at = place_menu[2,1] =  Button(f, label="+1")
b_st = place_menu[2,2] =  Button(f, label="-1")

# transition menu
l_en_t= transition_menu[1,1] =  Label(f, "Status: ")
l_en = transition_menu[1,2] =  Label(f, lift(
    (idx,m,h)-> if idx ≤ length(pn.P)
        "N.A."
    elseif isready(pn, TPNState(m,h), idx-length(pn.P))
        "ready"
    elseif isenabled(pn, TPNState(m,h), idx-length(pn.P))
        "enabled"
    else
        "disabled"
    end,
    p.attributes.selected,
    p.attributes.p_marking,
    p.attributes.t_marking
))
l_h_t = transition_menu[2,1] = Label(f, "Clock: ")
l_h = transition_menu[2,2] = Label(f, lift((idx,h)->if idx ≤ length(pn.P)
        "N.A."
    elseif isnothing(h[idx-length(pn.P)])
        "♮"
    else
        @sprintf "%.2f" h[idx-length(pn.P)]
    end,
    p.attributes.selected,
    p.attributes.t_marking
))
b_f  = transition_menu[3,1] = Button(f,label="Fire!")
b_fa  = transition_menu[3,2] = Button(f,label="Fire all necessary!")

# delay menu
delay = Observable(0.0)
l_dm_t = delay_menu[1,1] = Label(f, "Max. delay:")
l_dm   = delay_menu[1,2] = Label(f, lift(dm->(@sprintf "%.2f" dm),p.attributes.max_delay))
s_d    = delay_menu[2,1:2] = Slider(f, snap=true, startvalue=1.0, range=LinRange(0.0,1.0,11))
t_d    = delay_menu[3,1] = Textbox(f, displayed_string=lift(d->(@sprintf "%.2f" d),delay), validator=(s->(res=tryparse(Float64,s);!isnothing(res) && 0≤res≤p.attributes.max_delay[])), reset_on_defocus=false)
b_d    = delay_menu[3,2] = Button(f, label="Elapse!")
t_d.stored_string[] = @sprintf "%.2f" delay[]

on(t_d.stored_string) do v
    f = parse(Float64,v)
    if p.attributes.max_delay[] ≈ 0.0
        set_close_to!(s_d, 0.0)
    else
        set_close_to!(s_d, f/p.attributes.max_delay[])
    end
end

on(s_d.value) do v
    delay[] = v*p.attributes.max_delay[]
end


# draw tokens
scatter!(ax, p.attributes.tokens; color=:black, markersize=0.1, markerspace=SceneSpace, zindex=1)

# change tokens with nodes, markings
onany(p.attributes.node_pos, p.attributes.p_marking) do pos,mm
    _tokens = Point2f[]
    for (i,t) in enumerate(mm)
        if t == 0
        elseif t == 1
            push!(_tokens, pos[i])
        else
            r = 0.125
            ϕ = LinRange(0,2π,t+1)[1:end-1]
            append!(_tokens, Ref(pos[i]) .+ r.*Point2f.(sin.(ϕ), cos.(ϕ)))
        end
    end
    empty!(p.attributes.tokens[])
    append!(p.attributes.tokens[],_tokens)
    p.attributes.tokens[] = p.attributes.tokens[]
end

# change display based on marking
onany(p.attributes.p_marking, p.attributes.t_marking) do m,t
    x = TPNState(p.attributes.p_marking[],p.attributes.t_marking[])
    # color enabled transitions
    colors = fill(:gray, length(pn.P)+length(pn.T))
    for t_id in eachindex(pn.T)
        if isready(pn, x, t_id)
            colors[length(pn.P)+t_id] = :red
        elseif isenabled(pn, x, t_id)
            colors[length(pn.P)+t_id] = :blue
        end
    end
    p.attributes.node_color[] .= colors
    p.attributes.node_color[] = p.attributes.node_color[]
end

# change display & menu when selected node changes
on(p.attributes.selected) do idx
    p.attributes.node_attr.glowcolor[] .= :transparent
    p.attributes.node_attr.glowcolor[][idx]=:red
    p.attributes.node_attr.glowcolor[] = p.attributes.node_attr.glowcolor[]
end

# fire transition on button click
on(b_f.clicks) do n
    idx = p.attributes.selected[]-length(pn.P)
    x = TPNState(p.attributes.p_marking[],p.attributes.t_marking[])
    if idx <= 0
        println("No transition selected!")
    elseif !isready(pn, x, idx)
        println("Transition not ready!")
    else
        println("Fired $(pn.T[idx])")
        fire!(pn, x, idx)
        p.attributes.p_marking[] .= x.m
        p.attributes.t_marking[] .= x.h
        p.attributes.p_marking[] = p.attributes.p_marking[] 
        p.attributes.t_marking[] = p.attributes.t_marking[] 
    end
end

# fire all transitions on button click
on(b_fa.clicks) do n
    x = TPNState(p.attributes.p_marking[],p.attributes.t_marking[])
    transitions = fire_necessary!(pn, x)

    p.attributes.p_marking[] .= x.m
    p.attributes.t_marking[] .= x.h
    p.attributes.p_marking[] = p.attributes.p_marking[] 
    p.attributes.t_marking[] = p.attributes.t_marking[] 
    
    if isempty(transitions)
        println("No ready transitions to fire.")
    else
        println("Fired transitions "*join(transitions, ",", " and ")*".")
    end
end

# add or remove tokens on button click
on(b_at.clicks) do n
    idx = p.attributes.selected[]
    if idx > length(pn.P)
        println("No place selected!")
    else
        p.attributes.p_marking[][idx] += 1
        x = TPNState(p.attributes.p_marking[],p.attributes.t_marking[])
        for t_id ∈ eachindex(pn.T)
            update_clock!(pn,x,t_id)
        end
        p.attributes.t_marking[] .= x.h
        p.attributes.p_marking[] = p.attributes.p_marking[]
        p.attributes.t_marking[] = p.attributes.t_marking[]
    end
end
on(b_st.clicks) do n
    idx = p.attributes.selected[]
    if idx > length(pn.P)
        println("No place selected!")
    elseif p.attributes.p_marking[][idx] <= 0
        println("Number of tokens must be ≥ 0")
    else
        p.attributes.p_marking[][idx] -= 1
        x = TPNState(p.attributes.p_marking[],p.attributes.t_marking[])
        for t_id ∈ eachindex(pn.T)
            update_clock!(pn,x,t_id)
        end
        p.attributes.t_marking[] .= x.h
        p.attributes.p_marking[] = p.attributes.p_marking[]
        p.attributes.t_marking[] = p.attributes.t_marking[]
    end
end

# delay time on button press
on(b_d.clicks) do n
    if delay[] > p.attributes.max_delay[]
        println("Delay $(delay[]) exceeds max. delay of $(p.attributes.max_delay).")
    else
        x = TPNState(p.attributes.p_marking[],p.attributes.t_marking[])
        elapse!(pn, x, delay[])
        println("Elapsed τ=$(delay[]).")
        p.attributes.p_marking[] .= x.m
        p.attributes.t_marking[] .= x.h
        p.attributes.p_marking[] = p.attributes.p_marking[] 
        p.attributes.t_marking[] = p.attributes.t_marking[] 
    end
end

# save on button press
on(b_save.clicks) do n
    save(t_fn.stored_string[], ax.scene)
end

# add click interaction to select node
function select_node(idx, event, axis)
    p.attributes.selected=idx
end
register_interaction!(ax, :nodeclick, NodeClickHandler(select_node))

# add drag interaction
deregister_interaction!(ax, :rectanglezoom)
register_interaction!(ax, :ndrag, NodeDrag(p))

hidedecorations!(ax)

# force execution of some updates
p.attributes.selected[] = p.attributes.selected[]
p.attributes.p_marking[] = p.attributes.p_marking[]
p.attributes.t_marking[] = p.attributes.t_marking[]

return f,ax,p
end

