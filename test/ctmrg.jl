using TensorNetworkAD2
using TensorNetworkAD2:  getχ, SquareCTMRGRuntime, CTMRGRuntime, getd, gets, ctmrg
using Test, Random
using Zygote

@testset "IPEPS" begin
    @test SquareLattice() isa AbstractLattice
    sq = SquareIPEPS(randn(3,3,3,3,2))
    @test sq isa SquareIPEPS
    @test getd(sq) == 3
    @test gets(sq) == 2
    @test_throws DimensionMismatch SquareIPEPS(randn(3,3,4,3,2))
end

@testset "runtime" begin
    # NOTE: previous initializet and initializec should have been tested.
    a = randn(ComplexF64, 2,2,2,2)
    env1 = SquareCTMRGRuntime(a, Val(:random), 10)
    @test env1 isa CTMRGRuntime
    @test getχ(env1) == 10
    env2 = SquareCTMRGRuntime(a, Val(:raw), 10)
    @test env1 isa CTMRGRuntime
    @test getχ(env2) == 10
end

@testset "ctmrg unit test" begin
    rt = SquareCTMRGRuntime(randn(2,2,2,2), Val(:random), 10)
    rt = ctmrg(rt; tol=1e-6, maxit=10)
    @test rt isa CTMRGRuntime
    @test getχ(rt) == 10
end

