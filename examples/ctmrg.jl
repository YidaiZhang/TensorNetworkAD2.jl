using TensorNetworkAD2
using TensorNetworkAD2.OMEinsum

h = hamiltonian(Heisenberg())
ipeps = SquareIPEPS(rand(2,2,2,2,2));
ipeps = TensorNetworkAD2.indexperm_symmetrize(ipeps);
TensorNetworkAD2.energy(h,ipeps, χ=20, tol=1e-6,maxit=100)

function energy(h::AbstractArray{T,4}, ipeps::IPEPS; χ::Int, tol::Real, maxit::Int) where T
    ipeps = indexperm_symmetrize(ipeps)  # NOTE: this is not good
    D = getd(ipeps)^2
    s = gets(ipeps)
    ap = ein"abcdx,ijkly -> aibjckdlxy"(ipeps.bulk, conj(ipeps.bulk))
    ap = reshape(ap, D, D, D, D, s, s)
    a = ein"ijklaa -> ijkl"(ap)

    rt = SquareCTMRGRuntime(a, Val(:raw), χ)
    rt  = ctmrg(rt; tol=tol, maxit=maxit)
    e = expectationvalue(h, ap, rt)
    return e
end
energy(h,ipeps, χ=20, tol=1e-6,maxit=100)