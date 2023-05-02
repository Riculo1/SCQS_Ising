using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

clear_cache()

airy = [2.33811, 4.08795, 5.52056, 6.7867144, 7.94413, 9.02265]  # zeros of airy function

g = 0.09
hz = 0.01
momenta = [0]
D = 100
num = 5

E0 = (1-2g)/2
A = (hz)^(2/3)*(g/(1-2g))^(1/3)

Ψ = InfiniteMPS(2, D)
H = transverse_field_ising(; hx=g, hz=hz)
Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS())

energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), momenta, Ψ_ground, envs; num = num)
energies = real.(energies)

plt = scatter(zeros(length(energies)), vec(energies), markershape=:circle, label="Calculated excitations", xaxis=false, ylabel="ΔE")

theory = (airy .* A .+ 2*E0)[1:length(energies)]
scatter!(plt, 2 .* ones(length(theory)), theory, markershape=:square, label="2E0 + A*zn")

savefig(plt, "AiryTest.png")
@info energies
@info theory
@show plt
