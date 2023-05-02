using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

#=
Comparison of the calculated dispersion relation for the first excited state with the the relation found in the sylabus "Strongly Correlated Quantum Systems".
For g < 0.5 a spin up and down ground state got added together to get the actual first excitation and not the second one.
=#

# Integrators
VUMPS_alg = VUMPS(maxiter=200)
QuasiparticleAnsatz_alg = QuasiparticleAnsatz()

gs = [0.01, 0.1, 0.2, 0.5, 1, 2]  # g values
momenta = range(-π, π, 25)
D = 20  # Bond dimension
Ψ = InfiniteMPS(2, D)

function theoreticalDispersion(k; J=1, g=0.01)  # Theoretical dispersion relation for the first excited state
    return J/2*sqrt((cos(k)-2*g)^2 + sin(k)^2)
end

function guessΨ(g, hz=0.1)  # Get a ground state which has a spin aligned according to hz
    H = transverse_field_ising(; hx=g)
    Hguess = transverse_field_ising(; hx=g, hz=hz)

    Ψ = InfiniteMPS(2, D)
    Ψhz, envs, δ = find_groundstate(Ψ, Hguess, VUMPS_alg)
    Ψguess, envsguess, δ = find_groundstate(Ψhz, H, VUMPS_alg)

    while hz*real(sum(expectation_value(Ψguess, σᶻ))) < 0  # hz and σᶻ same sign
        Ψ = InfiniteMPS(2, D)
        Ψhz, envs, δ = find_groundstate(Ψ, Hguess, VUMPS_alg)
        Ψguess, envsguess, δ = find_groundstate(Ψhz, H, VUMPS_alg)
    end
    return Ψguess, envsguess
end

plt = plot(; title="Dispersion relations", xlabel="Momenta", ylabel="ΔE", legend=:outertopright)

for (i,g) in enumerate(gs)
    H = transverse_field_ising(; hx=g)
    
    if g < 0.5
        Ψup, envsup = guessΨ(g)  # Spin up ground state
        Ψdown, envsdown = guessΨ(g, -0.1)  # Spin down ground state
        energies, quasiparticles = excitations(H, QuasiparticleAnsatz_alg, momenta, Ψup, envsup, Ψdown, envsdown)
    else
        Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
        energies, quasiparticles = excitations(H, QuasiparticleAnsatz_alg, momenta, Ψ_ground, envs)
    end

    energies = real.(energies)
    @info "completed g = $g"
    plot!(plt, momenta, energies; label="MPS g=$g", color=i, ls=:solid)
    plot!(plt, momenta, theoreticalDispersion.(momenta; g=g); label="Analytical g=$g", color=i, ls=:dash)
end

savefig(plt, "figuren/dispersion_theory.png")
@show plt