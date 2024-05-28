using TensorNetworkAD2
using CairoMakie
using TensorNetworkAD2.Enzyme
using TensorNetworkAD2.OMEinsum
using TensorNetworkAD2.LinearAlgebra

Dcut = 20
n = 20

lnZ = []
K = 0.5
t = 1.0/K
#T = Ising( K )
y = trg(K, Dcut, n);
#@show lnZ
println(K, " ", y/2^n)

