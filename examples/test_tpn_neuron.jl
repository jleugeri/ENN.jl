using ENN.Neurons, ENN.TimePetriNets, GLMakie, GraphMakie

neuron_sequential = Neuron((((((),[:A1,:A2],:A),),[:B1,:B2],:B),),[:C1,:C2],:C);
tpn=TPN(neuron_sequential)

f,ax,p = graphplot(tpn); f