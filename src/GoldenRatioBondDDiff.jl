using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

#=
Plot the absolute error with the golden ratio in function of the bond dimension at g = 0.5 for different hz.
=#

# Integrators
VUMPS_alg = VUMPS(maxiter=200, verbose=false)
QuasiparticleAnsatz_alg = QuasiparticleAnsatz()

g = 0.5
hzs = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]  # hz values
momenta = [0]
Ds = range(2, 60)  # Bond dimensions

Φ = (1 + sqrt(5))/2  # Golden ratio
plt = plot(title="Difference with Golden Ratio", xlabel="D", ylabel="|E₂/E₁ - Φ|", legend=:outertopright, yaxis=:log)

for hz in hzs
    ratio = []
    for (i,D) in enumerate(Ds)
        Ψ = InfiniteMPS(2, D)
        H = transverse_field_ising(; hx=g, hz=hz)

        Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
        energies, quasiparticles = excitations(H, QuasiparticleAnsatz_alg, momenta, Ψ_ground, envs, num=2)

        append!(ratio, abs(real(energies[2])/real(energies[1]) - Φ))
        @info "completed D = $D"
        clear_cache()
    end

    @info "completed hz = $hz"
    plot!(plt, Ds, ratio; label="hz=$hz")

end

plot!(plt, xlabel="Bond dimension")

savefig(plt, "Golden Ratio bondD diff.png")

@show plt