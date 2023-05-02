using TensorKit, TensorOperations
using MPSKit, MPSKitModels
using Plots
using Statistics
using LinearAlgebra
using JLD2

# Integrators
VUMPS_alg = VUMPS(verbose=false)
QuasiparticleAnsatz_alg = QuasiparticleAnsatz()

airy = [2.33811, 4.08795, 5.52056, 6.7867144, 7.94413, 9.02265]  # zeros of airy function

# Parameters
bondD = [40]  # Bond dimension
gs = collect(range(0, 0.5, length=41))  # g values
gs[1] = 0.001  # Should not start at 0
gs[end] = 0.499  # Should not end at 0.5, because division by 0
hzs = collect(range(0, 0.1, length=41)) # hz values
append!(hzs, collect(range(0, 1, length=41))[6:end])
hzs[1] = 0.0002  # Should not start at 0
momenta = 0  # Momenta
num = 5  # Number of excitations
devision = 10  # Value to determine when the continuum has arrived

for D in bondD
    diff = [zeros(Float64, length(hzs), length(gs)) for _ in range(1,num)]  # Heatmap for every excitation
    data = []

    for n in range(1,num)  # Check if data has already been calculated, if so, load the file
        if isfile("data/excitation=$n, D=$D.jld2")
            push!(data, load("data/excitation=$n, D=$D.jld2"))
        end
    end

    for (i,g) in enumerate(gs)
        for (j,hz) in enumerate(hzs)

            # Airy
            E0 = (1-2g)/2
            A = (hz)^(2/3)*(g/(1-2g))^(1/3)
            theory = (airy .* A .+ 2*E0)[1:num]

            if length(data) >= num && haskey(data[num], "g=$g, hz=$hz")  # Value has already been calculated, use this value
                @info "D = $D, g = $g and hz = $hz has already been computed"

                energies = []
                for l in range(1,num)
                    append!(energies, data[l]["g=$g, hz=$hz"])
                end

                gap = abs(energies[1] - energies[2])

                for l in range(1,num)
                    if l <= 2  # First two excitations
                        diff[l][j,i] = abs((energies[l]-theory[l])/theory[l])
                    elseif abs(energies[l] - energies[l-1]) > gap/devision  # Are these excitations?
                        diff[l][j,i] = abs((energies[l]-theory[l])/theory[l])
                    else  # The continuum
                        diff[l][j,i] = NaN
                    end
                end

            else  # Calculate the value for every excitation
                @info "computing D = $D, g = $g and hz = $hz"
                # Excitations
                Ψ = InfiniteMPS(2, D)
                H = transverse_field_ising(; hx=g, hz=hz)
                Ψ_ground, envs, δ = find_groundstate(Ψ, H, VUMPS_alg)
                energies, quasiparticles = excitations(H, QuasiparticleAnsatz_alg, [momenta], Ψ_ground, envs, num=num)
                energies = real.(energies)

                # Compare Airy and calculated excitations
                gap = abs(energies[1] - energies[2])

                for (l,E) in enumerate(energies)

                    if l <= 2  # First two excitations
                        diff[l][j,i] = abs((energies[l]-theory[l])/theory[l])
                    elseif abs(energies[l] - energies[l-1]) > gap/devision  # Are these excitations?
                        diff[l][j,i] = abs((energies[l]-theory[l])/theory[l])
                    else  # The continuum
                        diff[l][j,i] = NaN
                    end

                    if !isfile("data/excitation=$l, D=$D.jld2")  # File does not exist
                        save("data/excitation=$l, D=$D.jld2", Dict("g=$g, hz=$hz" => E))
                    else  # Entry does not exist
                        jldopen("data/excitation=$l, D=$D.jld2", "r+") do file
                            write(file, "g=$g, hz=$hz", E)
                        end
                    end

                end
                
                clear_cache()
            end
        end
        println("")
    end

    # Make a nice heatmap from the errors for every excitation
    overlay = zeros(length(hzs), length(gs))
    for (m,to_plot) in enumerate(diff)
        overlay = max.(overlay, to_plot)
        plt = heatmap(gs, hzs, to_plot, clim=(0, 0.1), c=:thermal, xlabel = "g", ylabel = "hz", colorbar_title = "error", title="Relative error for excitation $m\nWith k = $momenta and bond dimension D = $D")
        savefig(plt, "Heatmap Airy D=$D exc=$m HD.png")
    end
    plt2 = heatmap(gs, hzs, overlay, c=:thermal, clim=(0, 0.1), xlabel = "g", ylabel = "hz", colorbar_title = "error", title="Maximum relative error\nWith k = $momenta and bond dimension D = $D")

    savefig(plt2, "Heatmap Airy overlay HD.png")
    println("")
end

@info "DONE :D"