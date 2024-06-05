using TensorNetworkAD2
using TensorNetworkAD2.LinearAlgebra
using TensorNetworkAD2.OMEinsum
using TensorNetworkAD2.Zygote

using Optim
using Test

D, χ = 2, 30
max_iters = 100
tol = 1e-12

A = Array{ComplexF64, 5}(randn(D, D, D, D, 2))
A += permutedims(A, (1,4,3,2,5))    # left-right
A += permutedims(A, (3,2,1,4,5))    # up-down
A += permutedims(A, (2,1,4,3,5))    # diagonal
A += permutedims(A, (4,3,2,1,5))    # rotation
#=
A += ein"iabcd -> iadcb"(A) # left-right
A += ein"iabcd -> icbad"(A) # up-down
A += ein"iabcd -> ibadc"(A) # diagonal
A += ein"iabcd -> idcba"(A) # rotation
A /= norm(A)
=#

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
    braket = reshape(ein"abcdp, ijklq -> aibjckdlpq"(PEPS, conj(PEPS)), D^2, D^2, D^2, D^2, 2, 2)
    braket /= norm(braket)
    half = ein"((ab, aci), dbj, (de, fek)), ijklpq->clfpq"(corner, edge, edge, corner, edge, braket)

    peps_norm = ein"abcpp, abcqq->"(half, half)[]

    return real(ein"abcpq, (abcrs, pqrs)->"(half, half, h)[]/peps_norm)
end

#=
function fixed_point(PEPS)
    env = envionment(PEPS)
    h = hamiltonian(Heisenberg())
    @show energy = PEPS_energy_expectation(PEPS, h, env...)
    return energy
end
=#

corner = Array{ComplexF64, 2}(randn(χ, χ))
corner += corner'
corner /= norm(corner)

edge = Array{ComplexF64, 3}(randn(χ, χ, D^2))
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

function tmp_iterate(PEPS, c, e)
    bulk = reshape(ein"abcdp, ijklp -> aibjckdl"(PEPS, conj(PEPS)), D^2, D^2, D^2, D^2)
    corner = c
    edge = e

    S_old = Vector{ComplexF64}(rand(χ * D^2))

    for n in 1:max_iters
        println("Iteration $n")
        new_corner = ein"((ab,aci),dbj),ijkl -> cldk"(corner, edge, edge, bulk)
        new_corner_matrix = reshape(new_corner, χ * D^2, χ * D^2)
        new_corner_matrix += new_corner_matrix'
        U, S, _ = svd(new_corner_matrix)

        if (norm(S_old ./ S_old[1] - S./S[1]) < tol)
            println("Converged at iteration $n")
            break
        end
        S_old = S

        isometry = reshape(U[:, 1:χ], χ, D^2, χ)

        corner = ein"(aibj, aic), bjd -> cd"(new_corner, conj(isometry), isometry)
        corner += corner'
        corner /= norm(corner)

        new_edge = ein"abi,ijkl -> ajblk"(edge,bulk)

        edge = ein"(aibjk,aic),bjd -> cdk"(new_edge, conj(isometry), isometry)
        edge += ein"abi -> bai"(conj(edge))
        edge /= norm(edge)
    end

    energy = PEPS_energy_expectation(PEPS, hamiltonian(Heisenberg()), corner, edge)
    println("Energy: ", energy)
    
    return energy
end

while true 
    g = Zygote.gradient(x ->  begin
        energy = tmp_iterate(x, corner, edge)
    end , A)[1]
    A -= 1e-2 * g
end

let energy = x -> tmp_iterate(x, corner, edge)
    @show res = Optim.optimize(energy,
        Δ -> Zygote.gradient(energy,Δ)[1], A, LBFGS(m = 20), inplace = false)
end
