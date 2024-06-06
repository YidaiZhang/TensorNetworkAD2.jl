using TensorNetworkAD2
using CairoMakie
using TensorNetworkAD2.OMEinsum
using TensorNetworkAD2.LinearAlgebra
using TensorNetworkAD2.Zygote


Dcut = 20
n = 20

β_array = collect(0.40:0.001:0.50)

# compute the gradient
function calc_lnZ(β)
    lnZ = 0.0
    for K in β
        y = trg(K, Dcut, n)
        lnZ += y/2^n
    end
    return lnZ
end
println("lnZ = ", calc_lnZ(0.5))

grad = Zygote.gradient(calc_lnZ, β_array)


# grad2 = Zygote.gradient(x -> Zygote.gradient(calc_lnZ, x)[1], β_array)


# analyze the gradient
neg_grad_values = -grad[1]

# verify the result
grad_verify = Zygote.gradient(calc_lnZ, 0.5)


#= second derivative
grad_array = collect(first(grad))
dF = diff(grad_array)./diff(β_range)

S =  β_range[1:end-1].^2 .* dF

grad_co = collect(grad)
grad2 = Zygote.gradient(grad, β_range)

#= plotting the result
fig = Figure()
ax = Axis(fig[1, 1], xlabel="β", ylabel="specific_heat", title="reproduce")
lines!(ax, β_range[1:end-1], S, color=:black)
scatter!(ax, β_range[1:end-1], S, label="∂lnZ/∂β", color=:red)
axislegend(ax; position = :rt, labelsize = 12)

fig





# plotting
fig = Figure()
ax = Axis(fig[1, 1], xlabel="β", ylabel="energy density", title="reproduce")

lines!(ax, β_range, neg_grad_values, color=:black)
scatter!(ax, β_range, neg_grad_values, color=:red, label="-∂lnZ/∂β")
axislegend(ax; position = :rt, labelsize = 12)
fig

# save figure
save("specific_heat_ad.png", fig)
=#










