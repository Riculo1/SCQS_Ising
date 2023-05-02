using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

VUMPS_alg = VUMPS(maxiter=200)

gs = [0.01, 0.1, 0.2, 0.5, 1, 2]
momenta = range(-π, π, 25)
D=20
Ψ = InfiniteMPS(2, D)

function theoreticalDispersion(k; J=1, g=0.01)
    return J/2*sqrt((cos(k)-2*g)^2 + sin(k)^2)
end

plt = plot(; title="Dispersion relations", xlabel="Momenta", ylabel="ΔE", legend=:outertopright)

for (i,g) in enumerate(gs)
    H = transverse_field_ising(; hx=g)
    
    Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
    energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), momenta, Ψ_ground, envs)

    energies = real.(energies)
    @info "completed g = $g"
    plot!(plt, momenta, energies; label="MPS g=$g", color=i, ls=:solid)
    plot!(plt, momenta, theoreticalDispersion.(momenta; g=g); label="Analytical g=$g", color=i, ls=:dash)
end

savefig(plt, "figuren/dispersion_theory_wrong.png")
@show plt