module TensorNetworkAD2

using LinearAlgebra
using ITensors
using OMEinsum
using TensorOperations
using QuadGK
using Random
using IterTools: iterated
using Base.Iterators: take, drop
using Optim, LineSearches



export TRG
export trg
export trg_i
export ising_mpo, ising_free_energy
export hamiltonian, Heisenberg
export IPEPS, SquareIPEPS
export diaglocalhamiltonian, expectationvalue, optimiseipeps
export energy, indexperm_symmetrize
export state 


include("trg_to.jl")
include("trg_ome.jl")
include("trg_itensors.jl")
include("2d_ising.jl")
include("heisenbergmodel.jl")
include("ipeps.jl")
include("ctmrg.jl")
include("variation.jl")

end
