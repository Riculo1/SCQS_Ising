using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots
using JLD2
using CurveFit

#=
Plotting ΔE(hz) to maybe see 3/2 power scaling
g should NOT be 0 or excitations won't work. (product states don't work)
but the airy approximation is probably better when g is smaller?
=#

# Integrators
VUMPS_alg = VUMPS(verbose=false)
QuasiparticleAnsatz_alg = QuasiparticleAnsatz()

airy = [2.33811, 4.08795, 5.52056, 6.7867144, 7.94413, 9.02265]  # zeros of airy function

g = 0.0125
hzs = collect(range(0.005, 0.010, length=11))
momenta = 0
num = 5  # Number of excitations
D = 100  # Bond dimension
Ψ = InfiniteMPS(2, D)

data = []
for n in range(1,num)  # Check if data has already been calculated, if so, load the file
    if isfile("data/excitation=$n, D=$D.jld2")
        push!(data, load("data/excitation=$n, D=$D.jld2"))
    end
end

plt = plot(title="ΔE(hz) for excited states", xlabel="hz", ylabel="ΔE")

E0 = (1-2g)/2
A = (hzs).^(2/3).*(g/(1-2g))^(1/3)
for n in range(1, num)
    theory = (airy[n] .* A .+ 2*E0)
    plot!(hzs, theory, label="Theory excitation $n", ls=:dash, color=n)
end

ΔE = []
for hz in hzs
    if length(data) >= num && haskey(data[num], "g=$g, hz=$hz")  # Value has already been calculated, use this value
        @info "D = $D, g = $g and hz = $hz has already been computed"

        energies = []

        for n in range(1,num)
            append!(energies, data[n]["g=$g, hz=$hz"])
        end

        push!(ΔE, vec(real.(energies)))

    else
        @info "computing D = $D, g = $g and hz = $hz"
        H = transverse_field_ising(; hx=g, hz=hz)
        Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
        energies, quasiparticles = excitations(H, QuasiparticleAnsatz_alg, [momenta], Ψ_ground, envs, num=num)
        energies = real.(energies)

        push!(ΔE, vec(energies))

        for (l,E) in enumerate(energies)
            if !isfile("data/excitation=$l, D=$D.jld2")  # File does not exist
                save("data/excitation=$l, D=$D.jld2", Dict("g=$g, hz=$hz" => E))
            else  # Entry does not exist
                jldopen("data/excitation=$l, D=$D.jld2", "r+") do file
                    write(file, "g=$g, hz=$hz", E)
                end
            end
        end
    end
end

for excitation in range(1, num)
    plot!(plt, hzs, [ΔE[i][excitation] for i in range(1, length(hzs))], label="Computed excitation $excitation", color=excitation, legend=:outertopright);
end

E_list = []
hz_list = []
for n in range(1,num)
    append!(E_list, ([ΔE[i][n] for i in range(1, length(hzs))] .- 2*E0) ./airy[n] .*((1-2g)/g)^(1/3))
    append!(hz_list, hzs)
end

a, b = power_fit(hz_list, E_list)
@info "a = $a and b = $b in f(hz) = a*hz^b"

savefig(plt, "HzScaling D = $D.png")
@show plt