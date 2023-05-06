using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

#=
This code calculates the excitation energy for different bond dimensions.
It also compares this result with the theoretical equation using the zeros of the Airy function.
=#

VUMPS_alg = VUMPS(verbose=false)  # Ground state algorithm
QuasiparticleAnsatz_alg = QuasiparticleAnsatz()  # Excitation algorithm

g = 0.0125
hz = 0.0025
momenta = [0]
num = 5  # Number of excitations
Ds = range(10, 100, step=10)  # Bond dimensions

# Theoretical solutions with zeros of the Airy function
m0 = (1-2g)/2
A = (hz)^(2/3)*(g/(1-2g))^(1/3)
airy = [2.33811, 4.08795, 5.52056, 6.7867144, 7.94413, 9.02265, 10.0402, 11.0085, 11.936, 12.8288, 13.6915]  # zeros of airy function
theory = (airy .*A .+ 2*m0)[1:num]

to_plot = [[] for _ in range(1,num)]  # Array of energies for different bond dimensions for every excitation

for (i,D) in enumerate(Ds)
    Ψ = InfiniteMPS(2, D)
    H = transverse_field_ising(; hx=g, hz=hz)  # Hamiltionian

    Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
    energies, quasiparticles = excitations(H, QuasiparticleAnsatz_alg, momenta, Ψ_ground, envs, num=num)
    energies = real.(energies)

    for (j,E) in enumerate(energies)
        append!(to_plot[j], E)
    end
    @info "completed D = $D"
    clear_cache()
end

# Plot of energy in function of bond dimension
labels = reshape(["Excitation $n" for n in range(1, num)], (1, num))
plt1 = plot(Ds, to_plot; label=labels, title="Excitation energy", xlabel="Bond dimension", ylabel="E", legend=:outertopright)

savefig(plt1, "airy bond.png")

# Plot of relative error with theory in function of bond dimension
error = [[] for _ in range(1,num)]
for (i, exc) in enumerate(to_plot)
    error[i] = abs.(exc .- theory[i])./theory[i]
end

plt2 = plot(Ds, error; label=labels, title="Excitation energy error", xlabel="Bond dimension", ylabel="|(E-theory)/theory|", legend=:outertopright, yaxis=:log)

savefig(plt2, "airy bond error.png")