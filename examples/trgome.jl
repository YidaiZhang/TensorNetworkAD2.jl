using TensorNetworkAD2
using CairoMakie
using TensorNetworkAD2.Enzyme
using TensorNetworkAD2.OMEinsum
using TensorNetworkAD2.LinearAlgebra

Dcut = 20
n = 20
β = 0.39:0.002:0.51;

lnZ = []
for K in β
    #T = Ising( K )
    y = trg(K, Dcut, n);
    #@show lnZ
    println(K, " ", y/2^n)
    push!(lnZ,y/2^n)
end


# taking the first derivative
dF = - diff(lnZ)./diff(β)

# plotting
# Find the indices where β is between 0.40 and 0.50
indices = findall(x -> 0.40 <= x <= 0.50, β[1:end-1])

# Plot only the selected range
fig = Figure()
ax = Axis(fig[1, 1], xlabel="β", ylabel="energy density", title="reproduce")

lines!(ax, β[indices], dF[indices], label="energy density")
scatter!(ax, β[indices], dF[indices], label="-∂lnZ/∂β")

fig

# save figure
save("specific_heat.png", fig)

# taking the second derivative with β^2
dS = - β[1:end-2].^2 .* diff(dF)./diff(β[1:end-1])
# Find the indices where β is between 0.40 and 0.50
indices = findall(x -> 0.40 <= x <= 0.50, β[1:end-2])

# Plot only the selected range
fig2 = Figure()
ax2 = Axis(fig2[1, 1], xlabel="β", ylabel="specific heat", title="reproduce")

lines!(ax2, β[indices], dS[indices], label="specific heat")
scatter!(ax2, β[indices], dS[indices], label="β^2(∂^2lnZ/∂β^2)")

fig2

# save figure
save("specific_heat_to.png", fig2)

