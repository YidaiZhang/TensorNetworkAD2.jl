function energy(h::AbstractArray{T,4}, ipeps::IPEPS; χ::Int, tol::Real, maxit::Int) where T
    ipeps = indexperm_symmetrize(ipeps)  
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
    l = ein"(ab,bde),ica,(eg,gfk),cjfdlm -> ijklm"(corner,edge,edge,corner,edge,ap)
    e = ein"(abcij,abckl),ijkl -> "(l,l,h)[]
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