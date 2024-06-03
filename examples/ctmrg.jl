using TensorNetworkAD2
using TensorNetworkAD2.OMEinsum

h = hamiltonian(Heisenberg())

ipeps = SquareIPEPS(rand(2,2,2,2,2));

A = TensorNetworkAD2.indexperm_symmetrize(ipeps);
A
A = ein"abcdx,ijkly -> aibjckdlxy"(ipeps.bulk, conj(ipeps.bulk))
E = OMEinsum.ein"ijkl -> ijk"(B)
C = OMEinsum.ein"ijkl -> ij"(B)



C_new = ein"ij,ijk,ijk,klij -> ijkl"(C,E,E,B)
E_new = ein"ij,ijkl -> ik"(C,E,B)



