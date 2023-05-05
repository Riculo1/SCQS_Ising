using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

#=
Entanglement entropy for smallest hz at g = 0.5.
=#

algorithm = VUMPS()

g = 0.5
hz = 0

Ψ = InfiniteMPS(2, 100)
H = transverse_field_ising(; hx=g, hz=hz)
Ψ_groundstate, envs, δ = find_groundstate(Ψ, H, algorithm)

plt = entanglementplot(Ψ_groundstate)
savefig(plt, "magnetisation_bond.png")