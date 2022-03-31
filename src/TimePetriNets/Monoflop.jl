export Monoflop

function Monoflop(name::Symbol, duration::H; start_delay=zero(H), M=Int, add_reset=false) where H
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

    # add transitions and arcs
    t_start,t_stop = add_transition!.((("$(name)_$(trans)") for trans in (:start,:stop)),(start_delay,duration))
    push!(ΔF,
        (p_off, t_start, -1),
        (p_on,  t_start,  1),
        (p_off, t_stop,   1),
        (p_on,  t_stop,  -1)
    )

    # add reset-related transitions and arcs
    if add_reset
        # add transitions
        t_reset = add_transition!("$(name)_reset",zero(H),zero(H),1.0)
        # add state-changing arcs
        push!(ΔF,
            (p_on, t_reset, -1),
            (p_off, t_reset, 1)
        )
    end

    ΔF = sparse(collect.(zip(ΔF...))..., length(P), length(T))
    R  = isempty(R) ? zero(ΔF) : sparse(collect.(zip(R...))..., length(P), length(T))
    
    TPN(P,T,R,ΔF,eft,lft,trans_probs,m₀)
end

