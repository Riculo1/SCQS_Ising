using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

#=
Plot the ratio of the second and first exitation in function of g for different hz.
Check if this ratio approaches the golden ratio if g approaches 0.5.
=#

# Integrators
VUMPS_alg = VUMPS(maxiter=200, verbose=false)
QuasiparticleAnsatz_alg = QuasiparticleAnsatz()

gs = collect(range(0.4, 0.5, length=21))  # g values
hzs = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]  # hz values
momenta = [0]
D = 60  # Bond dimesnion
Ψ = InfiniteMPS(2, D)

# Parameters which got closest to the golden ratio
best = 0  # Golden ratio approximation
best_hz = 0  # hz
best_g = 0  # g

# Plot the golden ratio as a dotted red line
Φ = (1 + sqrt(5))/2
plt = plot([gs[1], gs[end]], [Φ, Φ], ls=:dash, label="Φ", title="Golden Ratio", xlabel="g", ylabel="E₂/E₁", legend=:outertopright, color=:red)

for hz in hzs
    ratio = []  # Ratio of the second and first excited state
    for (i,g) in enumerate(gs)
        H = transverse_field_ising(; hx=g, hz=hz)

        Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
        energies, quasiparticles = excitations(H, QuasiparticleAnsatz_alg, momenta, Ψ_ground, envs, num=2)

        append!(ratio, real(energies[2])/real(energies[1]))

        # Is this value closer to the golden ratio?
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