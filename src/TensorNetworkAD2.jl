module TensorNetworkAD2

using LinearAlgebra
using ITensors
using OMEinsum
using TensorOperations
using QuadGK

export TRG
export trg
export trg_i
export ising_mpo, ising_free_energy
export hamiltonian, Heisenberg

include("trg_to.jl")
include("trg_ome.jl")
include("trg_itensors.jl")
include("2d_ising.jl")
include("heisenbergmodel.jl")

end
