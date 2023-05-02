using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

VUMPS_alg = VUMPS(maxiter=200)

g = 0.5  # Is critical
k = 0
bondD= range(2, 32, step=2)

function theoreticalDispersion(k; J=1, g=0.01)
    return J/2*sqrt((cos(k)-2*g)^2 + sin(k)^2)
end

ΔE = []
for (i,D) in enumerate(bondD)
    Ψ = InfiniteMPS(2, D)

    H = transverse_field_ising(; hx=g)
    
    Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
    MPS_energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), [k], Ψ_ground, envs)

    MPS_energies = real(MPS_energies[1])
    analytical_energies = theoreticalDispersion(k; g=g)

    append!(ΔE, abs(MPS_energies - analytical_energies))
    @info "completed D = $D"
end

plt = plot(bondD, ΔE; title="Error in function of bond dimension", xlabel="Bond dimension", ylabel="ΔE", legend=false, yaxis=:log)
savefig(plt, "figuren/dispersionBondD.png")
@show plt