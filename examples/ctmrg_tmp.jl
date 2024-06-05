using TensorNetworkAD2
using TensorNetworkAD2.LinearAlgebra
using TensorNetworkAD2.OMEinsum
using TensorNetworkAD2.Zygote

using Test

D, χ = 2, 30
max_iters = 50
tol = 1e-10

A = Array{ComplexF64, 5}(randn(ComplexF64, 2, D, D, D, D))
A += ein"iabcd -> iadcb"(A) # left-right
A += ein"iabcd -> icbad"(A) # up-down
A += ein"iabcd -> ibadc"(A) # diagonal
A += ein"iabcd -> idcba"(A) # rotation
A /= norm(A)

@test norm(A - ein"iabcd -> iadcb"(A)) < 1e-10
@test norm(A - ein"iabcd -> icbad"(A)) < 1e-10

heis = hamiltonian(Heisenberg())

#=
function envionment(A)
    bulk = reshape(ein"pabcd, pijkl -> aibjckdl"(A, conj(A)), D^2, D^2, D^2, D^2)
    corner = Array{ComplexF64, 2}(randn(ComplexF64, χ, χ))
    corner += corner'
    corner /= norm(corner)
    #corner += transpose(corner)
    #corner /= norm(corner)
    edge = Array{ComplexF64, 3}(randn(ComplexF64, χ, χ, D^2))
    edge += ein"abi -> bai"(conj(edge))
    edge /= norm(edge)

    old_S = Array{ComplexF64, 1}(randn(ComplexF64, χ *D^2))

    for n in 1:max_iters
        print("Iteration $n ")
        new_corner = ein"ab,aci,dbj,ijkl -> cldk"(corner, edge, edge, bulk)
        new_corner_matrix = reshape(new_corner, χ * D^2, χ * D^2)
        new_corner_matrix += new_corner_matrix'
        U, S, _ = svd(new_corner_matrix)

        if (norm(old_S - S) < 1e-10)
            println("Converged at iteration $n")
            break
        end
        old_S = S

        isometry = reshape(U[:, 1:χ], χ, D^2, χ)
        tmp_corner = corner
        corner = ein"aibj, aic, bjd -> cd"(new_corner, conj(isometry), isometry)

        new_edge = ein"abi,ijkl -> ajblk"(edge,bulk)

        edge = ein"aibjk,aic,bjd -> cdk"(new_edge, isometry, conj(isometry))

        corner += corner'
        edge += ein"abi -> bai"(conj(edge))

        corner /= norm(corner)

        edge /= norm(edge)

        #println(ein"((ab, bci), (dc, del)), ((fe, fgk), (hg, hbj)), ijkl->"(corner, edge, corner, edge, corner, edge, corner, edge,  bulk)[])
    end

    return corner, edge
end
=#

#=
function PEPS_norm(PEPS, corner, edge)
    bulk = reshape(ein"pabcd, pijkl -> aibjckdl"(PEPS, conj(PEPS)), D^2, D^2, D^2, D^2)

    println("norm: ", ein"(ab, acr), (cdv, ed), efx, gf, ghw, hit, ji, jbs, rstu, vuwx->"(corner, edge, edge, corner, edge, corner, edge, edge, corner, edge, bulk, bulk)[])

    return ein"((ab, aci), (dc, del)), ((fe, fgk), (hg, hbj)), ijkl->"(corner, edge, corner, edge, corner, edge, corner, edge, bulk)[]
end
=#

function PEPS_energy_expectation(PEPS, h, corner, edge)
    braket = reshape(ein"pabcd, qijkl -> pqaibjckdl"(PEPS, conj(PEPS)), 2, 2, D^2, D^2, D^2, D^2)
    #braket /= norm(braket)
    half = ein"((ab, aci), dbj, (de, fek)), pqijkl->pqclf"(corner, edge, edge, corner, edge, braket)

    peps_norm = ein"ppabc, qqcba->"(half, half)[]

    return real(ein"pqabc, rscba, pqrs->"(half, half, h)[]/peps_norm)
end

env = envionment(A)

println(PEPS_norm(A, env...))

println(PEPS_energy_expectation(A, heis, env...))

using Optim
using LineSearches

function fixed_point(PEPS)
    env = envionment(PEPS)
    h = hamiltonian(Heisenberg())
    @show energy = PEPS_energy_expectation(PEPS, h, env...)
    return energy
end

corner = Array{ComplexF64, 2}(randn(ComplexF64, χ, χ))
corner += corner'
corner /= norm(corner)
edge = Array{ComplexF64, 3}(randn(ComplexF64, χ, χ, D^2))
edge += ein"abi -> bai"(conj(edge))
edge /= norm(edge)

function iterate(PEPS)
    new_corner = ein"ab,aci,dbj,ijkl -> cldk"(corner, edge, edge, PEPS)
    new_corner_matrix = reshape(new_corner, χ * D^2, χ * D^2)
    new_corner_matrix += new_corner_matrix'
    U, S, _ = svd(new_corner_matrix)

    isometry = reshape(U[:, 1:χ], χ, D^2, χ)
    corner = ein"aibj, aic, bjd -> cd"(new_corner, conj(isometry), isometry)

    new_edge = ein"abi,ijkl -> ajblk"(edge,bulk)

    edge = ein"aibjk,aic,bjd -> cdk"(new_edge, isometry, conj(isometry))

    corner += corner'
    edge += ein"abi -> bai"(conj(edge))

    corner /= norm(corner)

    edge /= norm(edge)
    
    return real(sum(corner))
end

corner = Array{ComplexF64, 2}(rand(ComplexF64,χ, χ))
edge = Array{ComplexF64, 3}(rand(ComplexF64,χ, χ, D^2))

function tmp_iterate(PEPS, c, e)
    bulk = reshape(ein"pabcd, pijkl -> aibjckdl"(PEPS, conj(PEPS)), D^2, D^2, D^2, D^2)
    corner = c
    edge = e

    for n in 1:30
        println("Iteration $n")
        new_corner = ein"ab,aci,dbj,ijkl -> cldk"(corner, edge, edge, bulk)
        new_corner_matrix = reshape(new_corner, χ * D^2, χ * D^2)
        new_corner_matrix += new_corner_matrix'
        U, S, _ = svd(new_corner_matrix)

        isometry = reshape(U[:, 1:χ], χ, D^2, χ)

        corner = ein"aibj, aic, bjd -> cd"(new_corner, conj(isometry), isometry)
        corner += corner'
        corner /= norm(corner)

        new_edge = ein"abi,ijkl -> ajblk"(edge,bulk)

        edge = ein"aibjk,aic,bjd -> cdk"(new_edge, isometry, conj(isometry))
        edge += ein"abi -> bai"(conj(edge))
        edge /= norm(edge)
    end

    energy = PEPS_energy_expectation(PEPS, hamiltonian(Heisenberg()), corner, edge)
    println("Energy: ", energy)
    
    return energy
end

while true 
    g = Zygote.gradient(x ->  begin
        energy= tmp_iterate(x, corner, edge)
    end , A)[1]
    A -= 1e-2 * g
end


#=
@show Optim.optimize(x -> begin
        energy, c, e = tmp_iterate(x, corner, edge)
    end, Δ -> Zygote.gradient(x ->  begin
        energy, c, e = tmp_iterate(x, corner, edge)
    end, A)[1], A, LBFGS(), tol = 1e-10)

optimizer = LBFGS(; m = 20,
        linesearch = LineSearches.HagerZhang(),
        P = nothing,
        precondprep = (P, x) -> nothing)
=#