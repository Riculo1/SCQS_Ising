using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots
using CurveFit

#=
Calculate the value of the critical exponent in function of the bond dimension and compare it with the real value
=#

# Integrators
VUMPS_alg = VUMPS(verbose=false)

N = 20
g_range = collect(range(0.45, 0.499, N))
gc = 0.5  # Critical point

d = 2  # physical dimension ℂ²
bondD = collect(range(0, 25, step=5))  # virtual dimension
bondD[1] = 2

Sz = @mpoham sum(σᶻ{i} for i in vertices(InfiniteChain(1)))

plt = plot(bondD, [1/8 for _ in bondD], color=:red, ls=:dash, title="Critical exponent", xlabel="Bond dimension", ylabel="Exponent", label="β")

βs = []
for D in bondD
    m_array = Array{Float64}(undef, N)
    for (i, g) in enumerate(g_range)
        @info "Calculating D = $D, g = $g"

        # ground state of new H
        Ψ = InfiniteMPS(d, D)
        H = transverse_field_ising(; hx=g)
        Ψ_groundstate, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
        
        # Magnetization
        m = expectation_value(Ψ_groundstate, Sz)
        m_array[i] = abs(real(m[1]))
    end

    a, β = power_fit((gc.-g_range)./gc, m_array)
    @info "fit parameters are a = $a and β = $β for the function m = a*((gc-g)/gc)^β"
    
    mag_fit = plot(g_range, m_array, ls=:dash, title="Magnetisation with D = $D", xlabel="g", ylabel="m", label="calculated")
    plot!(mag_fit, g_range, a*((gc.-g_range)./gc).^β, label="fit")
    savefig(mag_fit, "magnetisation fit for D=$D.png")

    append!(βs, β)
end
plot!(plt, bondD, βs, label="Computed exponent", color=1)

savefig(plt, "critical_exponent.png")
@show plt