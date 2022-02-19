export Monoflop

function Monoflop(name::Symbol, duration::H; start_delay=zero(H), reset_delay=zero(H), reset_duration=zero(H), M=Int, add_reset=false, reset_name=:reset) where H
    P = Symbol[]
    T = Symbol[]
    m₀= M[]
    eft=H[]
    lft=H[]
    trans_probs=Float64[]

    R =  Tuple{Int,Int,M}[]
    ΔF = Tuple{Int,Int,M}[]

    function add_place!(name,initial)
        push!(P, Symbol(name))
        push!(m₀, M(initial))
        length(P)
    end

    function add_transition!(name,e,l=e,p=1.0)
        push!(T, Symbol(name))
        push!(eft, H(e))
        push!(lft, H(l))
        push!(trans_probs, p)
        length(T)
    end

    # add places
    p_on,p_off = add_place!.((("$(name)_$(state)") for state in (on,:off)),(0,1))
    if add_reset
        p_reset_on,p_reset_off = 
            add_place!.((("$(name)_$(reset_name)_$(state)") for state in (:on,:off)),(0,1))
    end

    # add transitions and arcs
    t_start,t_stop = add_transition!.((("$(name)_$(trans)") for trans in (:start,:stop)),(start_delay,duration))
    push!(ΔF,
        (p_off, t_start, -1),
        (p_on,  t_start,  1),
        (p_off, t_stop,   1),
        (p_on,  t_stop,  -1)
    )

    # add reset-related transitions and arcs
    t_reset_start = nothing
    t_reset_stop = nothing
    t_reset = nothing
    if add_reset
        # add transitions
        t_reset_start,t_reset_stop,t_reset = 
            add_transition!.((("$(name)_reset_$(trans)") for trans in (:start,:stop,:reset)),(reset_delay, reset_duration, zero(H)),(reset_delay, reset_duration, zero(H)),(1.0,0.0,1.0))
        
        # add read arcs
        push!(R,(p_reset_on, t_reset, 1),(p_reset_off, t_start, 1))
        # add state-changing arcs
        push!(ΔF,
            (p_on, t_reset, -1),(p_off, t_reset, 1),
            (p_reset_on, t_reset_start, 1), (p_reset_off, t_reset_start, -1), 
            (p_reset_on, t_reset_stop, -1), (p_reset_off, t_reset_stop, 1), 
        )
    end

    ΔF = sparse(collect.(zip(ΔF...))..., length(P), length(T))
    R  = isempty(R) ? zero(ΔF) : sparse(collect.(zip(R...))..., length(P), length(T))
    
    TPN(P,T,R,ΔF,eft,lft,trans_probs,m₀)
end

