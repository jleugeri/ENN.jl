using ENN.Neurons, GLMakie

neuron1 = Neuron((((((),[:E1,:E2,:E3,:E4],:E),((),[:E1,:E2],:F),((),[:G1,:G2],:G),((),[:H1,:H2],:H)),[:E1,:E2],:B),((((),Symbol[:E1],:I),),Symbol[],:C),(((((((),[:K1,:K2],:K),((),[:L1,:L2],:L),((),[:A1,:A2],:M))),Symbol[],:J),),Symbol[:H1],:D)),[:A1,:A2],:A)
neuron2 = Neuron((((((),[:E1,:E2,:E3,:E4],:E),((),[:E1,:E2],:F),((),[:G1,:G2],:G),((),[:H1,:H2],:H)),[:E1,:E2],:B),((((),Symbol[:E1],:I),),Symbol[],:C),(((((((),[:K1,:K2],:K),((),[:L1,:L2],:L),((),[:A1,:A2],:M))),Symbol[],:J),),Symbol[:H1],:D)),[:A1,:A2],:A)

net = NeuralNetwork(
    Dict(:n1=>neuron1,:n2=>neuron2),
    Dict(:in1=>[(:n2, :E2, :excitatory)]),
    Dict(:n2=>:out1),
    Dict(:n1=>[(:n2, :E1, :excitatory)],:n2=>[(:n1, :E1, :excitatory),(:n1, :H1, :excitatory)])
)

f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())
res=netplot!(ax, net, linewidth=2)
display(f)
