using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

VUMPS_alg = VUMPS(maxiter=200, verbose=false)

gs = collect(range(0.4, 0.5, length=21))
hzs = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
momenta = [0]
D=60
Ψ = InfiniteMPS(2, D)
best = 0
best_hz = 0
best_g = 0
Φ = (1 + sqrt(5))/2
plt = plot([gs[1], gs[end]], [Φ, Φ], ls=:dash, label="Φ", title="Golden Ratio", xlabel="g", ylabel="E₂/E₁", legend=:outertopright, color=:red)

for hz in hzs
    ratio = []
    for (i,g) in enumerate(gs)
        H = transverse_field_ising(; hx=g, hz=hz)

        Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
        energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), momenta, Ψ_ground, envs, num=2)

        append!(ratio, real(energies[2])/real(energies[1]))
        if abs(real(energies[2])/real(energies[1]) - Φ) < abs(best - Φ)
            global best = real(energies[2])/real(energies[1])
            global best_hz = hz
            global best_g = g
        end
        @info "completed g = $g"
        clear_cache()
    end

    @info "completed hz = $hz"
    plot!(plt, gs, ratio; label="hz=$hz")

end
@info "The closest I got to the golden ratio was $best for hz = $best_hz at g = $best_g"

savefig(plt, "Golden Ratio for D=$D.png")

@show plt