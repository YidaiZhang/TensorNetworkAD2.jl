abstract type HamiltonianModel end

struct Heisenberg{T<:Real} <: HamiltonianModel
    Jz::T
    Jx::T
    Jy::T
end
Heisenberg() = Heisenberg(1.0,1.0,1.0)

const σx = Float64[0 1; 1 0]
const σy = ComplexF64[0 -1im; 1im 0]
const σz = Float64[1 0; 0 -1]
const id2 = Float64[1 0; 0 1]

function hamiltonian(model::Heisenberg)
    h = model.Jz * ein"ij,kl -> ijkl"(σz,σz) -
        model.Jx * ein"ij,kl -> ijkl"(σx, σx) -
        model.Jy * ein"ij,kl -> ijkl"(σy, σy)
    h = ein"ijcd,kc,ld -> ijkl"(h,σx,σx')
    real(h ./ 2)
end
