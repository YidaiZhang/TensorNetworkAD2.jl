using TensorNetworkAD2
using Test, Random
using Zygote
using LinearAlgebra: norm

@testset "autodiff" begin
    a = randn(10,10)
    @test Zygote.gradient(norm, a)[1] ≈ num_grad(norm, a)

end