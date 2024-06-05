abstract type AbstractLattice end

struct SquareLattice <: AbstractLattice end

struct IPEPS{LT<:AbstractLattice, T, N, AT<:AbstractArray{T, N}}
    bulk::AT
end
IPEPS{LT}(bulk::AT) where {LT,T,N,AT<:AbstractArray{T,N}} = IPEPS{LT,T,N,AT}(bulk)

const SquareIPEPS{T} = IPEPS{SquareLattice, T, 5}
function SquareIPEPS(bulk::AT) where {T,AT<:AbstractArray{T, 5}}
    size(bulk,1) == size(bulk,2) == size(bulk,3) == size(bulk,4) || throw(DimensionMismatch(
        "size of tensor error, should be `(d, d, d, d, s)`, got $(size(bulk))."))
    IPEPS{SquareLattice,T,5,AT}(bulk)
end
getd(ipeps::SquareIPEPS) = size(ipeps.bulk, 1)
gets(ipeps::SquareIPEPS) = size(ipeps.bulk, 5)

function indexperm_symmetrize(ipeps::SquareIPEPS)
    x = ipeps.bulk
    x += permutedims(x, (1,4,3,2,5)) # left-right
    x += permutedims(x, (3,2,1,4,5)) # up-down
    x += permutedims(x, (2,1,4,3,5)) # diagonal
    x += permutedims(x, (4,3,2,1,5)) # rotation
    return SquareIPEPS(x / norm(x))
end

struct CTMRGRuntime{LT,T,N,AT<:AbstractArray{T,N},CT,ET}
    bulk::AT
    corner::CT
    edge::ET
    function CTMRGRuntime{LT}(bulk::AT,
        corner::AbstractArray{T}, edge::AbstractArray{T}) where {LT<:AbstractLattice,T,N,AT<:AbstractArray{T,N}}
        new{LT,T,N,AT,typeof(corner), typeof(edge)}(bulk,corner,edge)
    end
end

const SquareCTMRGRuntime{T,AT} = CTMRGRuntime{SquareLattice,T,4,AT}
SquareCTMRGRuntime(bulk::AT,corner,edge) where {T,AT<:AbstractArray{T, 4}} = CTMRGRuntime{SquareLattice}(bulk,corner,edge)

const σx = Float64[0 1; 1 0]
const σy = ComplexF64[0 -1im; 1im 0]
const σz = Float64[1 0; 0 -1]
const id2 = Float64[1 0; 0 1]

function hamiltonian()
    h =  ein"(ij,kl) -> ijkl"(σz,σz) -
         ein"(ij,kl) -> ijkl"(σx, σx) -
         ein"(ij,kl) -> ijkl"(σy, σy)
    h = ein"ijcd,kc,ld -> ijkl"(h,σx,σx')
    h = real(h ./ 2)
end

function SquareCTMRGRuntime(bulk::AbstractArray{T,4}, env::Val, χ::Int) where T
    return SquareCTMRGRuntime(bulk, _initializect_square(bulk, env, χ)...)
end

function _initializect_square(bulk::AbstractArray{T,4}, env::Val{:random}, χ::Int) where T
    corner = randn(T, χ, χ)
    edge = randn(T, χ, size(bulk,1), χ)
    corner += adjoint(corner)
    edge += permutedims(conj(edge), (3,2,1))
    corner, edge
end

function _initializect_square(bulk::AbstractArray{T,4}, env::Val{:raw}, χ::Int) where T
    cinit = ein"ijkl -> ij"(bulk)
    tinit = ein"ijkl -> ijk"(bulk)
    corner = [i in CartesianIndices(cinit) ? cinit[i] : zero(T) for i in CartesianIndices(zeros(T, χ, χ))]
    edge = [i in CartesianIndices(tinit) ? tinit[i] : zero(T) for i in CartesianIndices(zeros(T, χ, size(bulk,1), χ))]
    return corner, edge
end

getχ(rt::CTMRGRuntime) = size(rt.corner, 1)
getD(rt::CTMRGRuntime) = size(rt.bulk, 1)

function ctmrg(rt::CTMRGRuntime; tol::Real, maxit::Integer)
    # initialize
    oldvals = fill(Inf, getχ(rt)*getD(rt))

    stopfun = StopFunction(oldvals, -1, tol, maxit)
    rt, vals = fixedpoint(res->ctmrgstep(res...), (rt, oldvals), stopfun)
    return rt
end


function ctmrgstep(rt::SquareCTMRGRuntime, vals)
    # grow
    bulk, corner, edge = rt.bulk, rt.corner, rt.edge
    D, χ = getD(rt), getχ(rt)
    cp = ein"((ad,iba),dcl),jkcb -> ijlk"(corner, edge, edge, bulk)
    tp = ein"(iam,jkla) -> ijklm"(edge,bulk)

    # renormalize
    cpmat = reshape(cp, χ*D, χ*D)
    cpmat += adjoint(cpmat)
    u, s, v = svd(cpmat)
    z = reshape(u[:, 1:χ], χ, D, χ)

    corner = ein"((abcd,abi),cdj) -> ij"(cp, conj(z), z)
    edge = ein"((abjcd,abi),dck) -> ijk"(tp, conj(z), z)

    vals = s ./ s[1]

    # indexperm_symmetrize
    corner += corner'
    edge += ein"ijk -> kji"(conj(edge))

    # normalize
    corner /= norm(corner)
    edge /= norm(edge)

    return SquareCTMRGRuntime(bulk, corner, edge), vals
end

function energy(h::AbstractArray{T,4}, ipeps::IPEPS; χ::Int, tol::Real, maxit::Int) where T
    ipeps = indexperm_symmetrize(ipeps)  
    D = getd(ipeps)^2
    s = gets(ipeps)
    ap = ein"(abcdx,ijkly) -> aibjckdlxy"(ipeps.bulk, conj(ipeps.bulk))
    ap = reshape(ap, D, D, D, D, s, s)
    a = ein"ijklaa -> ijkl"(ap)

    rt = SquareCTMRGRuntime(a, Val(:raw), χ)
    rt  = ctmrg(rt; tol=tol, maxit=maxit)
    e = expectationvalue(h, ap, rt)
    return e
end

function expectationvalue(h, ap, rt::SquareCTMRGRuntime)
    corner, edge = rt.corner, rt.edge
    ap /= norm(ap)
    l = ein"((ab,ica),(bde,eg)),cjfdlm,gfk -> ijklm"(corner,edge,edge,corner,ap,edge)
    e = ein"(abcij,abckl),ijkl -> "(l,l,h)[]
    n = ein"(ijkaa,ijkbb) -> "(l,l)[]
    return e/n
end

function fixedpoint(f, guess, stopfun)
    for state in iterated(f, guess)
        stopfun(state) && return state
    end
end

mutable struct StopFunction{T,S}
    oldvals::T
    counter::Int
    tol::S
    maxit::Int
end


function (st::StopFunction)(state)
    st.counter += 1
    st.counter > st.maxit && return true

    vals = state[2]
    diff = norm(vals - st.oldvals)
    diff <= st.tol && return true
    st.oldvals = vals

    return false
end

function optimiseipeps(ipeps::IPEPS{LT}, h; χ::Int, tol::Real, maxit::Int,
    optimargs = (),
    optimmethod = LBFGS(m = 20)) where LT
bulk = ipeps.bulk
let energy = x -> real(energy(h, IPEPS{LT}(x); χ=χ, tol=tol, maxit=maxit))
    res = optimize(energy,
        Δ -> Zygote.gradient(energy,Δ)[1], bulk, optimmethod, inplace = false, optimargs...)
end
end