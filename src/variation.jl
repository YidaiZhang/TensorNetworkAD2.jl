function diaglocalhamiltonian(diag::Vector)
    n = length(diag)
    h = ein"i -> ii"(diag)
    id = Matrix(I,n,n)
    reshape(h,n,n,1,1) .* reshape(id,1,1,n,n) .+ reshape(h,1,1,n,n) .* reshape(id,n,n,1,1)
end


function energy(h::AbstractArray{T,4}, ipeps::IPEPS; χ::Int, tol::Real, maxit::Int) where T
    ipeps = indexperm_symmetrize(ipeps)  # NOTE: this is not good
    D = getd(ipeps)^2
    s = gets(ipeps)
    ap = ein"abcdx,ijkly -> aibjckdlxy"(ipeps.bulk, conj(ipeps.bulk))
    ap = reshape(ap, D, D, D, D, s, s)
    a = ein"ijklaa -> ijkl"(ap)

    rt = SquareCTMRGRuntime(a, Val(:raw), χ)
    rt  = ctmrg(rt; tol=tol, maxit=maxit)
    e = expectationvalue(h, ap, rt)
    return e
end

function expectationvalue(h, ap, rt::SquareCTMRGRuntime) where T
    corner, edge = rt.corner, rt.edge
    ap /= norm(ap)
    l = ein"ab,ica,bde,cjfdlm,eg,gfk -> ijklm"(corner,edge,edge,ap,corner,edge)
    e = ein"abcij,abckl,ijkl -> "(l,l,h)[]
    n = ein"ijkaa,ijkbb -> "(l,l)[]
    return e/n
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