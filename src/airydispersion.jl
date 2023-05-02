using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

VUMPS_alg = VUMPS(maxiter=200)

airy = [2.33811, 4.08795, 5.52056, 6.7867144, 7.94413, 9.02265]  # zeros of airy function

function theoreticalDispersion(k; J=1, g=0.01)
    return J/2*sqrt((cos(k)-2*g)^2 + sin(k)^2)
end

g = 0.3
hz = 0.2
num = 5
momenta = range(0, π, 11)
D = 40
Ψ = InfiniteMPS(2, D)

H = transverse_field_ising(; hx=g, hz=hz)

Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), momenta, Ψ_ground, envs, num=num)

energies = real.(energies)
@info "completed g = $g"

labels = reshape(["Calculated Excitation $i" for i in range(1,num)], (1, num))
plt = plot(momenta, energies, label=labels, title="Dispersion relations for g = $g, hz = $hz", xlabel="Momenta", ylabel="ΔE", legend=:outertopright)

E0 = theoreticalDispersion.(momenta./2, g=g)
A = (hz)^(2/3)*(g/(1-2g))^(1/3)
for n in range(1,num)
    plot!(plt, momenta, 2 .* E0 .+ airy[n] * A, label="Approximated Excitation $n", color=n, ls=:dash)
end

savefig(plt, "figuren/dispersion_airy.png")
@show plt