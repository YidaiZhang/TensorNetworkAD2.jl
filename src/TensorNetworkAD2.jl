module TensorNetworkAD2

using LinearAlgebra
using ITensors
using OMEinsum
using TensorOperations
using QuadGK
using Random
using IterTools: iterated,imap
using Base.Iterators: take, drop



export TRG
export trg
export trg_i
export ising_mpo, ising_free_energy
export hamiltonian, Heisenberg
export IPEPS, SquareIPEPS
export energy, expectationvalue
export StopFunction
export fixedpoint



include("trg_to.jl")
include("trg_ome.jl")
include("trg_itensors.jl")
include("2d_ising.jl")
include("2d_heisenbergmodel.jl")
include("ctmrg.jl")


end
