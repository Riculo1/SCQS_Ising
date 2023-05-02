# Created 08/04/2023 by Ian Lateur
# Analysing the difference between analytic disperions relation
# and the one found with MPS with the quasiparticle ansatz

using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

function analytical_disperion(J, g, k)
    # Completely modified with 2's and 1/2's, CHECK THIS!!!
    return J*sqrt(g^2 + 1/4 - g*cos(k))
end

J = 1.0
g_tuple = (0.01, 0.1, 0.5, 1, 2)
hz = 0.0

d = 2  # physical dimension ℂ²
D = 20  # virtual dimension, CHECK IF THIS IS OK
Ψ = InfiniteMPS(d, D)

cur_plot = plot(; title="Disperion Relations: MPS and Analytical")

for (i, g) in enumerate(g_tuple)
    print("Calculating g = $g...\n")

    # ground state of new H
    H = transverse_field_ising(; J=J, hx=g, hz=hz)
    Ψ_groundstate, envs, δ = find_groundstate(Ψ, H, VUMPS(maxiter=200))
    
    # excitations
    momenta = range(-π, π, 24)  # will be multithreaded
    energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), momenta, Ψ_groundstate)

    # plot
    plot!(cur_plot, momenta, real.(energies); label="MPS g = $g", ls=:solid, color=i)  # color as integer :)
    plot!(cur_plot, momenta, analytical_disperion.(J, g, momenta); label="Analytical g = $g", ls=:dash, color=i)
end

@show cur_plot