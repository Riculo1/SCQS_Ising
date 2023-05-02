# Created 08/04/2023 by Ian Lateur
# Analysing the magnetization of ground states

using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

J = 1.0
hz_tuple = (0., 0.01, 0.1, 1.0)

N = 20
g_range = range(0, 1, N)

d = 2  # physical dimension ℂ²
D = 20  # virtual dimension, CHECK IF THIS IS OK
Ψ = InfiniteMPS(d, D)

Sz = @mpoham sum(σᶻ{i} for i in vertices(InfiniteChain(1)))

cur_plot = plot(; title="Magnetization", xlabel="g", ylabel="magnetization")

for hz in hz_tuple
    m_array = Array{Float64}(undef, N)

    for (i, g) in enumerate(g_range)
        print("Calculating hz = $hz, g = $g\n")

        # ground state of new H
        H = transverse_field_ising(; J=J, hx=g, hz=hz)
        Ψ_groundstate, envs, δ = find_groundstate(Ψ, H, VUMPS(maxiter=200))
        
        # Magnetization
        m = expectation_value(Ψ_groundstate, Sz)
        m_array[i] = abs(real(m[1]))

    end
    # plot
    plot!(cur_plot, g_range, m_array; label="hz = $hz")
end

@show cur_plot