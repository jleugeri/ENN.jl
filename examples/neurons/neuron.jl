using TimedAutomata, Base.Iterators
using GraphMakie
##

##
neuron_sequential = Neuron((((((),:A),),:B),),:C);
neuron_sequential_ta = TA(neuron_sequential)
neuron_sequential_p = prune(neuron_sequential_ta)

#zone = TAConstraint(neuron_sequential_p.clocks, zeros(TABound{Int},length(neuron_sequential_p.clocks)+1,length(neuron_sequential_p.clocks)+1))
zone = @constraint clk_1-clk_2==0 && clk_2-clk_3==0 && clk_3-clk_1==0 Int
TimedAutomata.set_clock_list!(zone, neuron_sequential_ta.clocks)
neuron_sequential_tts = TTS(neuron_sequential_p; initial_zone=zone)
@assert all(neuron_sequential_tts.states[end][1].components[1])
sequential_language_full=language(neuron_sequential_tts; final_states=neuron_sequential_tts.states[end:end], ignore_messages=(_)->false)
sequential_language_no_silent=language(neuron_sequential_tts; final_states=neuron_sequential_tts.states[end:end])
sequential_language_input=language(neuron_sequential_tts; final_states=neuron_sequential_tts.states[end:end], forbid_messages=(m)->m.symbol∈(:↻,), ignore_messages=(m)->m.direction!=MSG_IN)
sequential_language_no_reset_no_spike=language(neuron_sequential_tts; final_states=neuron_sequential_tts.states[end:end], forbid_messages=(m)->m.symbol∈(:spike,:↻))

##
neuron_parallel = Neuron((((),:A),((),:B)),:C);
neuron_parallel_ta = TA(neuron_parallel)
neuron_parallel_p = prune(neuron_parallel_ta)
neuron_parallel_tts = TTS(neuron_parallel_p; initial_zone=zone)
id = findfirst(((l,_),)->all(l.components[1]), neuron_parallel_tts.states)
parallel_language_full=language(neuron_parallel_tts; final_states=neuron_parallel_tts.states[id:id], ignore_messages=(_)->false)
parallel_language_no_silent=language(neuron_parallel_tts; final_states=neuron_parallel_tts.states[id:id])
parallel_language_input=language(neuron_parallel_tts; final_states=neuron_parallel_tts.states[id:id], forbid_messages=(m)->m.symbol∈(:↻,), ignore_messages=(m)->m.direction!=MSG_IN)
parallel_language_no_reset_no_spike=language(neuron_parallel_tts; final_states=neuron_parallel_tts.states[id:id], forbid_messages=(m)->m.symbol∈(:spike,:↻))

## Plotting

f,ax,_ = graphplot(neuron_sequential_ta);
hidedecorations!(ax)
display(f)
print("Press key to continue: ")
read(stdin, 1)
save(joinpath("examples","neurons", "figures","sequential_unpruned.png"), f)

f,ax,_ = graphplot(neuron_sequential_p);
hidedecorations!(ax)
display(f)
print("Press key to continue: ")
read(stdin, 1)
save(joinpath("examples","neurons", "figures","sequential_pruned.png"), f)

f,ax,_ = graphplot(neuron_parallel_p);
hidedecorations!(ax)
display(f)
print("Press key to continue: ")
read(stdin, 1)
save(joinpath("examples","neurons", "figures","parallel_pruned.png"), f)
