using TensorNetworkAD2
using TensorNetworkAD2.OMEinsum
using Zygote
using LinearAlgebra
using Optim

h = hamiltonian()

ipeps = SquareIPEPS(rand(2,2,2,2,2));

A = TensorNetworkAD2.indexperm_symmetrize(ipeps);
println("A = ", A)

e = TensorNetworkAD2.energy(h,ipeps, χ=20, tol=1e-6,maxit=100)

function en(ipeps)
    return TensorNetworkAD2.energy(h, ipeps, χ=20, tol=1e-6, maxit=100)
end

grad = Zygote.gradient(en, ipeps)

while true
    grad = Zygote.gradient(en, ipeps)
    ipeps.bulk = ipeps - 0.01*grad
    println("Energy = ", en(ipeps))
end

res = optimiseipeps(ipeps, h; χ=20, tol=1e-6, maxit=100,
        optimargs = (Optim.Options(f_tol=1e-6, show_trace=true),));
