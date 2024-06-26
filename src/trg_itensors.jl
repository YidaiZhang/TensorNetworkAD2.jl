function trg_i(T::ITensor; χmax::Int, nsteps::Int, cutoff=0.0, svd_alg="divide_and_conquer")
  sₕ, sᵥ = filterinds(T; plev=0)
  @assert hassameinds((sₕ, sₕ', sᵥ, sᵥ'), T)

  # Keep track of the partition function per site
  κ = 1.0
  for n in 1:nsteps
    Fₕ, Fₕ′ = factorize(
      T, (sₕ', sᵥ'); ortho="none", maxdim=χmax, cutoff, tags=tags(sₕ), svd_alg
    )

    s̃ₕ = commonind(Fₕ, Fₕ′)
    Fₕ′ *= δ(dag(s̃ₕ), s̃ₕ')

    Fᵥ, Fᵥ′ = factorize(
      T, (sₕ, sᵥ'); ortho="none", maxdim=χmax, cutoff, tags=tags(sᵥ), svd_alg
    )

    s̃ᵥ = commonind(Fᵥ, Fᵥ′)
    Fᵥ′ *= δ(dag(s̃ᵥ), s̃ᵥ')

    T =
      (Fₕ * δ(dag(sₕ'), sₕ)) *
      (Fᵥ * δ(dag(sᵥ'), sᵥ)) *
      (Fₕ′ * δ(dag(sₕ), sₕ')) *
      (Fᵥ′ * δ(dag(sᵥ), sᵥ'))

    sₕ, sᵥ = s̃ₕ, s̃ᵥ

    trT = abs((T * δ(sₕ, sₕ') * δ(sᵥ, sᵥ'))[])
    T = T / trT
    κ *= trT^(1 / 2^n)
  end
  return κ, T
end