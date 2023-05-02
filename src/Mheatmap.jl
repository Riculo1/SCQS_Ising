using KrylovKit, Plots, MPSKitModels
using TensorKit
using TensorOperations
using MPSKit

VUMPS_alg = VUMPS(maxiter=200, verbose=false)

bound = 1e-12
D = 20
Ψ = InfiniteMPS(2, D)

hz = collect(range(0, 1, length=21))
g = collect(range(0, 2, length=21))
g[1] = 0.01

function magnetisation(g, hz)
    H = transverse_field_ising(; hx=g, hz=hz)
    Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
    clear_cache()
    return abs(real(sum(expectation_value(Ψ_ground, σᶻ))))
end

plt = heatmap(g, hz, magnetisation, c = :thermal, xlabel = "g", ylabel = "hz", colorbar_title = "magnetisation", title="magnetisation(g,hz)")
savefig(plt, "figuren/Mheatmap.png")
@show plt