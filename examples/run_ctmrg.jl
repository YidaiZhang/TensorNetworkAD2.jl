using TensorNetworkAD2
using TensorNetworkAD2.OMEinsum


h = hamiltonian()

ipeps = SquareIPEPS(rand(2,2,2,2,2));

A = TensorNetworkAD2.indexperm_symmetrize(ipeps);

@time TensorNetworkAD2.energy(h,ipeps, χ=20, tol=1e-6,maxit=100)




