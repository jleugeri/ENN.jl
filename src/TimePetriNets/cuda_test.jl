using CUDA, SparseArrays

P   = Symbol[Symbol("P$(i)") for i in 1:num_places]
T   = Symbol[Symbol("T$(i)") for i in 1:num_transitions]
C   = sprand(Bool, num_places, num_transitions, 0.3)
ΔF  = sprand(Int, num_places, num_transitions, 0.3)
eft = zeros(Float64, num_transitions)
lft = ones(Float64, num_transitions)
m₀  = zeros(Int, num_places)
h₀  = zeros(Float64, num_transitions)
