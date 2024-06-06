function trg(K::Real, Dcut::Int, no_iter::Int)
    D = 2

    T = zeros(Float64, D, D, D, D)
    M = [[sqrt(cosh(K)) sqrt(sinh(K))];
         [sqrt(cosh(K)) -sqrt(sinh(K))];
         ]
    
    T = ein"(ai, aj), (ak, al) -> ijkl"(M, M, M, M)
    lnZ = 0.0

    for n in collect(1:no_iter)

        #println(n, " ", maximum(T), " ", minimum(T))
        maxval = maximum(abs.(T))
        T = T/maxval
        lnZ += 2.0^(no_iter-n+1)*log(maxval)

        D_new = min(D^2, Dcut)

        T1 = ein"urdl -> drul"(T)
        T11 = reshape(T1, (D^2, D^2))
        T2 = ein"urdl -> ldru"(T)
        T22 = reshape(T2, (D^2, D^2))
        F = svd(T11)
        F1 = reshape(F.U[:,1:D_new]*Diagonal(sqrt.(F.S[1:D_new])), (D, D, D_new))
        F3 = reshape(Diagonal(sqrt.(F.S[1:D_new]))*F.Vt[1:D_new, :], (D_new, D, D))
        F = svd(T22)
        F2 = reshape(F.U[:,1:D_new]*Diagonal(sqrt.(F.S[1:D_new])), (D, D, D_new))
        F4 = reshape(Diagonal(sqrt.(F.S[1:D_new]))*F.Vt[1:D_new, :], (D_new, D, D))
        T_new = ein"(war,abu),(lbg,dgw) -> ruld"(F1, F2, F3, F4)

        D = D_new
        T = T_new
    end

    trace = ein"ijij -> "(T)[]
    lnZ += log(trace)
end



