module TensorNetworkAD2

using LinearAlgebra
using ITensors
using OMEinsum
using TensorOperations
using QuadGK
using Random
using IterTools: iterated,imap
using Base.Iterators: take, drop
using Optim, LineSearches
using Zygote


export TRG
export trg
export trg_i
export ising_mpo, ising_free_energy
export hamiltonian, Heisenberg
export IPEPS, SquareIPEPS
export energy, expectationvalue, optimiseipeps
export StopFunction
export fixedpoint



include("autodiff.jl")
include("trg_to.jl")
include("trg_ome.jl")
include("trg_itensors.jl")
include("2d_ising.jl")
include("2d_heisenbergmodel.jl")
include("ipeps.jl")
include("ctmrg.jl")
include("energy.jl")
include("fixpoint.jl")

end
