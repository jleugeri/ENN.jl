using SparseArrays, ENN

pn=TPN(
    [:is_trig1, :is_on1, :is_off1, :is_trig2, :is_on2, :is_off2, :counter],
    [:in1, :in2, :on1, :off1, :on2, :off2],
    sparse([
        0; 0; 0;  0; 0; 0;  0 ;; # fire first input
        0; 0; 0;  0; 0; 0;  0 ;; # fire second input
        1; 0; 1;  0; 0; 0;  0 ;; # turn 1 on
        0; 1; 0;  0; 0; 0;  1 ;; # turn 1 off
        0; 0; 0;  1; 0; 1;  0 ;; # turn 2 on
        0; 0; 0;  0; 1; 0;  1    # turn 2 off
    ]),
    sparse([
        1; 0; 0;  0; 0; 0;  0 ;; # fire first input
        0; 0; 0;  1; 0; 0;  0 ;; # fire second input
       -1; 1;-1;  0; 0; 0;  1 ;; # turn 1 on and increase counter
        0;-1; 1;  0; 0; 0; -1 ;; # turn 1 off and decrease counter
        0; 0; 0; -1; 1;-1;  1 ;; # turn 2 on and increase counter
        0; 0; 0;  0;-1; 1; -1    # turn 2 off and decrease counter
    ]),
    [0.0, 0.0, 0.0, 1.0, 0.0, 1.0],
    [Inf, Inf, 0.0, 1.0, 0.0, 1.0],
    [0, 0, 1, 0, 0, 1, 0]
)

f,ax,p = graphplot(pn); f
