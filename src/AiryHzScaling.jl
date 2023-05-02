# Plotting ΔE(hz) to maybe see 3/2 power scaling

# g should NOT be 0 or excitations won't work. (product states don't work)
# but the airy approximation is probably better when g is smaller?

using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

# airy = [2.33811, 4.08795, 5.52056, 6.7867144, 7.94413, 9.02265]  # zeros of airy function

g = 0.01  # g should NOT be 0 exactly, but Airy approximation is probably worse if
hzs = collect(range(0, 2, step=0.02)) 
hzs[1] = 0.01  # don't use 0 exactly, otherwise no symmetry breaking
num = 4

D = 4  # TODO: check
Ψ = InfiniteMPS(2, D)

ΔE = []
for hz in hzs
    @info "Computing hz = $hz"
    H = transverse_field_ising(; hx=g, hz=hz)
    Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS())

    energies, quasiparticles = excitations(H, QuasiparticleAnsatz(), [0], Ψ_ground, envs; num = num)
    push!(ΔE, vec(real.(energies)))
end

plt = plot(title="ΔE(hz) for excited states")
for excitation in range(1, num)
    plot!(plt, hzs, [ΔE[i][excitation] for i in range(1, length(hzs))], label="Excitation $excitation");
end

savefig(plt, "HzScaling D=4 ipv 30.png")
@show plt