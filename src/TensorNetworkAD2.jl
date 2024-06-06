module TensorNetworkAD2

using LinearAlgebra
using ITensors
using OMEinsum
using TensorOperations
using QuadGK
using Random
using IterTools: iterated,imap
using Base.Iterators: take, drop
using Zygote
using Optim

import LinearAlgebra



export TRG
export trg
export trg_i
export ising_mpo, ising_free_energy
export hamiltonian
export IPEPS, SquareIPEPS
export energy, expectationvalue
export StopFunction
export fixedpoint
export num_grad
export optimiseipeps
export svd_back, trg_svd



include("trg_to.jl")
include("trg_ome.jl")
include("trg_itensors.jl")
include("2d_ising.jl")
include("ctmrg.jl")
include("autodiff.jl")


end
