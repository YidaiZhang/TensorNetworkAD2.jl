using Test
using TensorNetworkAD2
using TensorNetworkAD2:  energy, expectationvalue, optimiseipeps,
                       hamiltonian, indexperm_symmetrize, num_grad
using OMEinsum, Zygote, Random
using LinearAlgebra: svd, norm
using Optim, LineSearches



@testset "gradient" begin
    Random.seed!(0)
    h = hamiltonian()
    ipeps = SquareIPEPS(randn(2,2,2,2,2))
    a = indexperm_symmetrize(ipeps)
    gradzygote = first(Zygote.gradient(a) do x
        energy(h,x; χ=4, tol=0, maxit=100)
    end).bulk
    gradnum = num_grad(a.bulk, δ=1e-3) do x
        energy(h, SquareIPEPS(x); χ=4, tol=0, maxit=100)
    end

    @test isapprox(gradzygote, gradnum, atol=1e-3)
end

@testset "complex" begin
    Random.seed!(2)
    h = hamiltonian()
    ipeps = SquareIPEPS(randn(2,2,2,2,2))
    a = indexperm_symmetrize(ipeps)
    ca = SquareIPEPS(a.bulk .+ 0im)
    @test energy(h,a; χ=4, tol=1e-12, maxit=100) ≈ energy(h,ca; χ=4, tol=1e-12, maxit=100)
    ϕ = exp(1im * rand()* 2π)
    ca = SquareIPEPS(a.bulk .* ϕ)
    @test energy(h,ca; χ=4, tol=1e-12, maxit=100) ≈ energy(h,a; χ=4, tol=1e-12, maxit=100)

    gradzygote = first(Zygote.gradient(a) do x
        real(energy(h,x; χ=4, tol=1e-12,maxit=100))
    end)

    ca = SquareIPEPS(a.bulk .+ 0im)
    @test gradzygote.bulk ≈ first(Zygote.gradient(ca) do x
        real(energy(h,x; χ=4, tol=1e-12,maxit=100))
    end).bulk

    Random.seed!(2)
    # real
    h = hamiltonian()
    ipeps = SquareIPEPS(randn(2,2,2,2,2))
    a = indexperm_symmetrize(ipeps)
    res1 = optimiseipeps(a, h; χ=20, tol=1e-12, maxit=100,
        optimargs = (Optim.Options(f_tol=1e-6, store_trace = true, show_trace=false),));

    # complex
    ipeps = SquareIPEPS(randn(2,2,2,2,2) .+ randn(2,2,2,2,2) .* 1im)
    a = indexperm_symmetrize(ipeps)
    res2 = optimiseipeps(a, h; χ=20, tol=1e-12, maxit=100,
        optimargs = (Optim.Options(f_tol=1e-6,store_trace = true,  show_trace=false, allow_f_increases=true),),
        optimmethod = Optim.LBFGS(
            m = 10,
            alphaguess = LineSearches.InitialStatic(alpha=1, scaled=true),
            linesearch = LineSearches.Static())
        );
    @test isapprox(minimum(res1), minimum(res2), atol = 1e-3)
end