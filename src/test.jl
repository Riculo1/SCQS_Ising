using TensorKit, TensorOperations
using MPSKit, MPSKitModels

g = 1.0
hz = 0.0

H = transverse_field_ising(; hz=hz, hx=g)
Ψ = InfiniteMPS(2, 20)
Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS())

energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), [0], Ψ_ground; num=1)

# println(sum(expectation_value(Ψ_ground, σᶻ)))
# println(sum(expectation_value(InfiniteMPS(quasiparticles[1].VLs), σᶻ)))
println(expectation_value(Ψ_ground, H))
println(expectation_value(InfiniteMPS(quasiparticles[1].VLs), H))
println(energies[1])
