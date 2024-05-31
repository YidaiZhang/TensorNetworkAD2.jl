using TensorNetworkAD2
using TensorNetworkAD2.OMEinsum

h = hamiltonian(Heisenberg())
ipeps = SquareIPEPS(rand(2,2,2,2,2));
ipeps = TensorNetworkAD.indexperm_symmetrize(ipeps);
TensorNetworkAD.energy(h,ipeps, χ=20, tol=1e-6,maxit=100)