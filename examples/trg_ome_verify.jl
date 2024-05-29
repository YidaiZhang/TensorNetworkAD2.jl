using TensorNetworkAD2
using CairoMakie
using TensorNetworkAD2.Enzyme
using TensorNetworkAD2.OMEinsum
using TensorNetworkAD2.LinearAlgebra


lnZ = []
K = 0.5
t = 1.0/K
M = [[sqrt(cosh(K)) sqrt(sinh(K))];
         [sqrt(cosh(K)) -sqrt(sinh(K))];
         ]
T = ein"(ai, aj), (ak, al) -> ijkl"(M, M, M, M)
lp = exp(0.5)+exp(-0.5)
lm = exp(0.5)-exp(-0.5)
lp^2/2
lm^2/2
lp*lm/2


Dcut = 20
n = 20

#T = Ising( K )
y = trg(K, Dcut, n);
#@show lnZ
println(K, " ", y/2^n)

D = 2
T_rand = rand(2,2,2,2)
T11 = reshape(T_rand, (D^2, D^2))