# Comparing Airy spectrum with gained spectrum at k=0
# g should NOT be 0 or excitations won't work. (product states don't work)
# but the airy approximation is probably better when g is smaller

# Should also plot ΔE(hz) 

using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

airy = [2.33811, 4.08795, 5.52056, 6.7867144, 7.94413, 9.02265]  # zeros of airy function

g = 0.01  # g should NOT be 0 exactly, but Airy approximation is probably worse if
hz = 1.0  # TODO: vary

D = 100  # TODO: check
num = 5
momenta = [0]  # TODO: is it possible to have momentum with Airy?

Ψ = InfiniteMPS(2, D)
H = transverse_field_ising(; hx=g, hz=hz)
Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS())

energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), momenta, Ψ_ground; num = num)
energies = real.(energies)


plt = scatter(zeros(length(energies)), vec(energies), markershape=:circle, label="Calculated excitations", xaxis=false, ylabel="ΔE")

scaled_airy = (airy  ./ airy[1] .* energies[1])[1:length(energies)]
scatter!(plt, ones(length(scaled_airy)), scaled_airy, markershape=:diamond, label="Airy (scaled)")

shifted_scaled_airy = ((airy .- airy[1]) ./ (airy[2] - airy[1]) .* energies[1])[1:length(energies)+1]
scatter!(plt, 2 .* ones(length(shifted_scaled_airy)), shifted_scaled_airy, markershape=:square, label="Airy (shifted & scaled)")

savefig(plt, "AiryTest.png")
@show plt