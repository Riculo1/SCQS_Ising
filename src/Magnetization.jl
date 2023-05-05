using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots

#=
Analysing the magnetization of ground states
=#

# Integrators
VUMPS_alg = VUMPS(verbose=false)

J = 1.0
hz_tuple = (0., 0.01, 0.05, 0.1)

N = 21
g_range = collect(range(0, 1, N))
g_range[1] = 0.001

d = 2  # physical dimension ℂ²
D = 20  # virtual dimension
Ψ = InfiniteMPS(d, D)

Sz = @mpoham sum(σᶻ{i} for i in vertices(InfiniteChain(1)))

cur_plot = plot(; title="Magnetization", xlabel="g", ylabel="magnetization")

for hz in hz_tuple
    m_array = Array{Float64}(undef, N)

    for (i, g) in enumerate(g_range)
        @info "Calculating hz = $hz, g = $g"

        # ground state of new H
        H = transverse_field_ising(; J=J, hx=g, hz=hz)
        Ψ_groundstate, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
        
        # Magnetization
        m = expectation_value(Ψ_groundstate, Sz)
        m_array[i] = abs(real(m[1]))

        clear_cache()

    end
    # plot
    plot!(cur_plot, g_range, m_array; label="hz = $hz")
end

savefig(cur_plot, "magnetisation.png")
@show cur_plot