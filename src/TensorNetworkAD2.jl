module TensorNetworkAD2

using LinearAlgebra
using TensorOperations
using Enzyme
using OMEinsum

export TRG
export trg
export trg_i
export ising_mpo, ising_free_energy

include("trg_to.jl")
include("trg_ome.jl")
include("trg_itensors.jl")
include("2d_ising.jl")

end
