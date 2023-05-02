using KrylovKit, Plots, MPSKitModels
using TensorKit
using TensorOperations
using MPSKit

#=
Create a heatmap with g on the x-axis and hz on the y-axis.
The color indicates the absolute value of the magnetisation
=#

# Integrators
VUMPS_alg = VUMPS(maxiter=200, verbose=false)

D = 20  # Bond dimension
Ψ = InfiniteMPS(2, D)

hz = collect(range(0, 1, length=21))
g = collect(range(0, 2, length=21))
g[1] = 0.01  # g should not be 0

Sz = @mpoham sum(σᶻ{i} for i in vertices(InfiniteChain(1)))

function magnetisation(g, hz)
    H = transverse_field_ising(; hx=g, hz=hz)
    Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
    clear_cache()
    return abs(real(sum(expectation_value(Ψ_ground, Sz))))
end

plt = heatmap(g, hz, magnetisation, c = :thermal, xlabel = "g", ylabel = "hz", colorbar_title = "magnetisation", title="magnetisation(g,hz)")
savefig(plt, "Mheatmap.png")
@show plt