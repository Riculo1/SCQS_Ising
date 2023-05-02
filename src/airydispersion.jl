using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

#=
This code calculates the dispersion relation for multiple excitations.
It also compares this result with the theoretical approximation.
=#

VUMPS_alg = VUMPS(maxiter=200)  # Ground state algorithm
QuasiparticleAnsatz_alg = QuasiparticleAnsatz()  # Excitation state algorithm

airy = [2.33811, 4.08795, 5.52056, 6.7867144, 7.94413, 9.02265]  # zeros of airy function

function theoreticalDispersion(k; J=1, g=0.01)  # Theoretical dispersion relation for the first excited state
    return J/2*sqrt((cos(k)-2*g)^2 + sin(k)^2)
end

g = 0.3
hz = 0.2
num = 5  # Number of excitations
momenta = range(0, π, 11)
D = 40  # Bond dimension
Ψ = InfiniteMPS(2, D)

H = transverse_field_ising(; hx=g, hz=hz)  # Hamiltonian

Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
energies, quasiparticles = excitations(H, QuasiparticleAnsatz_alg, momenta, Ψ_ground, envs, num=num)

energies = real.(energies)

labels = reshape(["Calculated Excitation $i" for i in range(1,num)], (1, num))
plt = plot(momenta, energies, label=labels, title="Dispersion relations for g = $g, hz = $hz", xlabel="Momenta", ylabel="ΔE", legend=:outertopright)

m0 = theoreticalDispersion.(momenta./2, g=g)
A = (hz)^(2/3)*(g/(1-2g))^(1/3)
for n in range(1,num)
    plot!(plt, momenta, 2 .* m0 .+ airy[n] * A, label="Approximated Excitation $n", color=n, ls=:dash)
end

savefig(plt, "figuren/dispersion_airy.png")
@show plt