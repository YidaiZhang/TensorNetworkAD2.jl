@Zygote.nograd StopFunction
@Zygote.nograd _initializect_square

@Zygote.adjoint function IPEPS{LT,T,N,AT}(bulk) where {LT,T,N,AT}
    IPEPS{LT,T,N,AT}(bulk), dy -> (dy.bulk,)
end

@Zygote.adjoint function CTMRGRuntime{LT}(bulk::AT,
        corner::AbstractArray{T}, edge::AbstractArray{T}) where {LT<:AbstractLattice,T,N,AT<:AbstractArray{T,N}}
        return CTMRGRuntime{LT}(bulk,corner,edge), dy->(dy.bulk, dy.corner, dy.edge)
end

@Zygote.adjoint function LinearAlgebra.norm(A::AbstractArray, p::Real = 2)
    n = norm(A,p)
    back(Δ) = let n = n
                    (Δ .* A ./ (n + eps(0f0)),)
                end
    return n, back
end

num_grad(f, K::Real; δ::Real = 1e-5) = (f(K+δ/2) - f(K-δ/2))/δ

function num_grad(f, a::AbstractArray; δ::Real = 1e-5)
    map(CartesianIndices(a)) do i
        foo = x -> (ac = copy(a); ac[i] = x; f(ac))
        num_grad(foo, a[i], δ = δ)
    end
end

function trg_svd(t, dmax, tol)
    d1, d2, d3, d4 = size(t)
    tmat = reshape(t, d1*d2, d3*d4)
    u, s, v = LinearAlgebra.svd(tmat)
    dmax = min(searchsortedfirst(s, tol, rev=true), dmax, length(s))
    FS = s[1:dmax]
    sqrtFSp = sqrt.(FS)
    u = reshape(ein"ij,j -> ij"(u[:,1:dmax],  sqrtFSp), (d1, d2, dmax))
    v = reshape(ein"ij,i -> ij"(copy(v')[1:dmax,:], sqrtFSp), (dmax, d3, d4))

    return u, v
end

struct ZeroAdder end
Base.:+(a, zero::ZeroAdder) = a
Base.:+(zero::ZeroAdder, a) = a
Base.:-(a, zero::ZeroAdder) = a
Base.:-(zero::ZeroAdder, a) = -a
Base.:-(zero::ZeroAdder) = zero

mpow2(a::AbstractArray) = a .^ 2

Zygote.@adjoint function LinearAlgebra.svd(A)
    res = LinearAlgebra.svd(A)
    res, function (dy)
        dU, dS, dVt = dy
        return (svd_back(res.U, res.S, res.V, dU, dS, dVt === nothing ? nothing : dVt'),)
    end
end


function svd_back(U::AbstractArray, S::AbstractArray{T}, V, dU, dS, dV; η::Real=1e-40) where T
    all(x -> x isa Nothing, (dU, dS, dV)) && return nothing
    η = T(η)
    NS = length(S)
    S2 = mpow2(S)
    Sinv = @. S/(S2+η)
    F = S2' .- S2
    F ./= (mpow2(F) .+ η)

    res = ZeroAdder()
    if !(dU isa Nothing)
        UdU = U'*dU
        J = F.*(UdU)
        res += (J+J')*LinearAlgebra.Diagonal(S) + LinearAlgebra.Diagonal(1im*imag(LinearAlgebra.diag(UdU)) .* Sinv)
    end
    if !(dV isa Nothing)
        VdV = V'*dV
        K = F.*(VdV)
        res += LinearAlgebra.Diagonal(S) * (K+K')
    end
    if !(dS isa Nothing)
        res += LinearAlgebra.Diagonal(dS)
    end

    res = U*res*V'

    if !(dU isa Nothing) && size(U, 1) != size(U, 2)
        res += (dU - U* (U'*dU)) * LinearAlgebra.Diagonal(Sinv) * V'
    end

    if !(dV isa Nothing) && size(V, 1) != size(V, 2)
        res = res + U * LinearAlgebra.Diagonal(Sinv) * (dV' - (dV'*V)*V')
    end
    res
end