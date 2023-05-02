using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

#=
Difference between the calculated first excitation energy and the theoretical one in function of the bond dimension.
=#

# Integrators
VUMPS_alg = VUMPS(maxiter=200)
QuasiparticleAnsatz_alg = QuasiparticleAnsatz()

g = 0.5  # Is critical
k = 0  # Momentum
bondD= range(2, 32, step=2)  # Bond dimensions

function theoreticalDispersion(k; J=1, g=0.01)  # Theoretical dispersion relation for the first excited state
    return J/2*sqrt((cos(k)-2*g)^2 + sin(k)^2)
end

ΔE = []  # Difference of calculated and theoretical first excitation energies for differend bond dimensions
for (i,D) in enumerate(bondD)
    Ψ = InfiniteMPS(2, D)

    H = transverse_field_ising(; hx=g)
    
    Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
    MPS_energies, quasiparticles = excitations(H, QuasiparticleAnsatz_alg, [k], Ψ_ground, envs)

    MPS_energies = real(MPS_energies[1])
    analytical_energies = theoreticalDispersion(k; g=g)

    append!(ΔE, abs(MPS_energies - analytical_energies))
    @info "completed D = $D"
end

plt = plot(bondD, ΔE; title="Error in function of bond dimension", xlabel="Bond dimension", ylabel="ΔE", legend=false, yaxis=:log)
savefig(plt, "figuren/dispersionBondD.png")
@show plt