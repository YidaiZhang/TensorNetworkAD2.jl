using TensorNetworkAD2
using TensorNetworkAD2.ITensors


β = 0.5
d = 2
s = Index(d)
sₕ = addtags(s, "horiz")
sᵥ = addtags(s, "vert")
T = ising_mpo(sₕ, sᵥ, β)

χmax = 20
nsteps = 20
κ, T = trg_i(T; χmax=χmax, nsteps=nsteps, svd_alg="divide_and_conquer")

κ_exact = exp(-β * ising_free_energy(β))
@show κ, κ_exact
@show abs(κ - κ_exact)

lnZ = log(κ)