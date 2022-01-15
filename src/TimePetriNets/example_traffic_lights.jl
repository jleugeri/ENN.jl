pn = TimePetriNet(
    [:red_1,:yellow_1,:green_1, :red_2, :yellow_2, :green_2, :free], 
    [:ry1,:yg1,:gr1, :ry2,:yg2,:gr2], 
    sparse([
        1;0;0 ; 0;0;0 ; 1 ;; # traffic light 1 red->yellow
        0;1;0 ; 0;0;0 ; 0 ;; # traffic light 1 yellow->green
        0;0;1 ; 0;0;0 ; 0 ;; # traffic light 1 green->red
        0;0;0 ; 1;0;0 ; 1 ;; # traffic light 2 red->yellow
        0;0;0 ; 0;1;0 ; 0 ;; # traffic light 2 yellow->green
        0;0;0 ; 0;0;1 ; 0    # traffic light 2 green->red
    ]), 
    sparse([
        -1; 1; 0 ;  0; 0; 0 ; -1 ;; # traffic light 1 red->yellow
         0;-1; 1 ;  0; 0; 0 ;  0 ;; # traffic light 1 yellow->green
         1; 0;-1 ;  0; 0; 0 ;  1 ;; # traffic light 1 green->red
         0; 0; 0 ; -1; 1; 0 ; -1 ;; # traffic light 2 red->yellow
         0; 0; 0 ;  0;-1; 1 ;  0 ;; # traffic light 2 yellow->green
         0; 0; 0 ;  1; 0;-1 ;  1    # traffic light 2 green->red
    ]), 
    ones(6),
    ones(6),
    [5;0;0;0;1;0;0]
)

f,ax,p = graphplot(pn); f

#=
trange = 1:60*length(states)
record(f, "traffic_light.mp4", trange; framerate=30) do tt
    mm = states[1+div(tt-1,60)]
    animate_step!(p, pn, mm)
    update!(f.scene)
end
=#