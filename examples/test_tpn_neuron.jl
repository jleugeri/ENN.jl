using ENN.Neurons, ENN.TimePetriNets, GLMakie, GraphMakie

n1 = Neuron((((((),[:A1,:A2],:A),),[:B1,:B2],:B),),[:C1,:C2],:C)
n2 = Neuron((((),[:B1,:B2],:B),),[:trig],:A)
#tpn=TPN(neuron_sequential)

net = NeuralNetwork(
    Dict(
        :n1 => n1,
        :n2 => n2
    ),
    Dict(
        :n1 => [(:n1, :A, :A1), (:n2, :A, :trig)],
        :n2 => [(:n1, :B, :inh), (:n1, :B, :B2)],
    )
)
net_tpn = TPN(net)


f,ax,p = graphplot(net_tpn); f
