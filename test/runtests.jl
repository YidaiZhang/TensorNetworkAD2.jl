using TensorNetworkAD2
using Test
using Random

@testset "TensorNetworkAD2.jl" begin
    @testset "autodiff" begin
        println("autodiff tests running...")
        include("autodiff.jl")
    end

@testset "ctmrg" begin
    println("ctmrg tests running...")
    include("ctmrg.jl")
end

@testset "fixedpoint" begin
    println("fixedpoint tests running...")
    include("fixpoint.jl")
end

@testset "svd" begin
    println("svd tests running...")
    include("svd.jl")
end

@testset "ipeps" begin
    println("ipeps tests running...")
    include("ipeps.jl")
end

end