using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

VUMPS_alg = VUMPS(verbose=false)

g = 0.25
hz = 0.0002
momenta = [0]
num = 5
Ds = range(10, 100, step=5)

E0 = (1-2g)/2
A = (hz)^(2/3)*(g/(1-2g))^(1/3)
airy = [2.33811, 4.08795, 5.52056, 6.7867144, 7.94413, 9.02265, 10.0402, 11.0085, 11.936, 12.8288, 13.6915]  # zeros of airy function

theory = (airy .*A .+ 2*E0)[1:num]

to_plot = [[] for _ in range(1,num)]

for (i,D) in enumerate(Ds)
    Ψ = InfiniteMPS(2, D)
    H = transverse_field_ising(; hx=g, hz=hz)

    Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
    energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), momenta, Ψ_ground, envs, num=num)
    energies = real.(energies)

    for (j,E) in enumerate(energies)
        append!(to_plot[j], E)
    end
    @info "completed D = $D"
    clear_cache()
end

labels = reshape(["Excitation $n" for n in range(1, num)], (1, num))
plt1 = plot(Ds, to_plot; label=labels, title="Excitation energy", xlabel="Bond dimension", ylabel="E", legend=:outertopright)

savefig(plt1, "airy bond test.png")

error = [[] for _ in range(1,num)]
for (i, exc) in enumerate(to_plot)
    error[i] = abs.(exc .- theory[i])./theory[i]
end

plt2 = plot(Ds, error; label=labels, title="Excitation energy error", xlabel="Bond dimension", ylabel="|(E-theory)/theory|", legend=:outertopright, yaxis=:log)

savefig(plt2, "airy bond error test.png")