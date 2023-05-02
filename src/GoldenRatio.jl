using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

VUMPS_alg = VUMPS(maxiter=200, verbose=false)

gs = collect(range(0.4, 0.5, length=21))
hzs = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
momenta = [0]
D=40
Ψ = InfiniteMPS(2, D)

Φ = (1 + sqrt(5))/2
plt = plot([gs[1], gs[end]], [Φ, Φ], ls=:dash, label="Φ", title="Golden Ratio", xlabel="g", ylabel="E₂/E₁", legend=:outertopright, color=:red)

for hz in hzs
    ratio = []
    for (i,g) in enumerate(gs)
        H = transverse_field_ising(; hx=g, hz=hz)

        Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
        energies2, quasiparticles = excitations(H, QuasiparticleAnsatz(), momenta, Ψ_ground, envs, num=2)

        append!(ratio, real(energies2[2])/real(energies2[1]))
        @info "completed g = $g"
        clear_cache()
    end

    @info "completed hz = $hz"
    plot!(plt, gs, ratio; label="hz=$hz")

end

savefig(plt, "Golden Ratio for D=$D.png")

@show plt