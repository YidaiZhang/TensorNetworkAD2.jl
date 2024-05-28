module TensorNetworkAD2

using LinearAlgebra
using TensorOperations
using Enzyme
using OMEinsum

export TRG
export trg


include("trg_to.jl")
include("trg_ome.jl")


end
