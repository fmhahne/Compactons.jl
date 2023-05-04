using Printf
using JLD2, CodecZlib
using UnPack
using KinkOscillonInteraction

mkpath("data/kink_oscillon_scattering_energies")

let l = 1.0, V = 0.75
    energies = Dict()

    for α in 0.00:0.01:0.99
        for v₀ in 0.00:0.01:0.99
            print("Simulating l=$l, V=$V, α=$α, v₀=$v₀ … ")
            @unpack E₁, E₂, E₃ = simulation(KinkOscillon(; l, V, α, v₀); dx=5e-3, sampling=2)
            energies[(l=l, V=V, α=α, v₀=v₀)] = (E₁=E₁, E₂=E₂, E₃=E₃)
            println("done")
        end
    end

    filename = "data/kink_oscillon_scattering_energies/l=$l,V=$V.jld2"
    jldsave(filename, true; energies)
end

let l = 1.0, α = 0.0
    energies = Dict()

    for V in 0.00:0.01:0.99
        for v₀ in 0.00:0.01:0.99
            print("Simulating l=$l, V=$V, α=$α, v₀=$v₀ … ")
            @unpack E₁, E₂, E₃ = simulation(KinkOscillon(; l, V, α, v₀); dx=5e-3, sampling=2)
            energies[(l=l, V=V, α=α, v₀=v₀)] = (E₁=E₁, E₂=E₂, E₃=E₃)
            println("done")
        end
    end

    filename = "data/kink_oscillon_scattering_energies/l=$l,alpha=$α.jld2"
    jldsave(filename, true; energies)
end

let l = 1.0, v₀ = 0.0
    energies = Dict()

    for V in 0.00:0.01:0.99
        for α in 0.00:0.01:0.99
            print("Simulating l=$l, V=$V, α=$α, v₀=$v₀ … ")
            @unpack E₁, E₂, E₃ = simulation(KinkOscillon(; l, V, α, v₀); dx=5e-3, sampling=2)
            energies[(l=l, V=V, α=α, v₀=v₀)] = (E₁=E₁, E₂=E₂, E₃=E₃)
            println("done")
        end
    end

    filename = "data/kink_oscillon_scattering_energies/l=$l,v0=$v₀.jld2"
    jldsave(filename, true; energies)
end

let V = 0.75, α = 0.0
    energies = Dict()

    for l in 0.50:0.025:3.00
        for v₀ in 0.00:0.01:0.99
            print("Simulating l=$l, V=$V, α=$α, v₀=$v₀ … ")
            @unpack E₁, E₂, E₃ = simulation(KinkOscillon(; l, V, α, v₀); dx=5e-3, sampling=2)
            energies[(l=l, V=V, α=α, v₀=v₀)] = (E₁=E₁, E₂=E₂, E₃=E₃)
            println("done")
        end
    end

    filename = "data/kink_oscillon_scattering_energies/V=$V,alpha=$α.jld2"
    jldsave(filename, true; energies)
end

let V = 0.75, v₀ = 0.0
    energies = Dict()

    for l in 0.50:0.025:3.00
        for α in 0.00:0.01:0.99
            print("Simulating l=$l, V=$V, α=$α, v₀=$v₀ … ")
            @unpack E₁, E₂, E₃ = simulation(KinkOscillon(; l, V, α, v₀); dx=5e-3, sampling=2)
            energies[(l=l, V=V, α=α, v₀=v₀)] = (E₁=E₁, E₂=E₂, E₃=E₃)
            println("done")
        end
    end

    filename = "data/kink_oscillon_scattering_energies/V=$V,v0=$v₀.jld2"
    jldsave(filename, true; energies)
end

let α = 0.00, v₀ = 0.0
    energies = Dict()

    for l in 0.50:0.025:3.00
        for V in 0.00:0.01:0.99
            print("Simulating l=$l, V=$V, α=$α, v₀=$v₀ … ")
            @unpack E₁, E₂, E₃ = simulation(KinkOscillon(; l, V, α, v₀); dx=5e-3, sampling=2)
            energies[(l=l, V=V, α=α, v₀=v₀)] = (E₁=E₁, E₂=E₂, E₃=E₃)
            println("done")
        end
    end

    filename = "data/kink_oscillon_scattering_energies/alpha=$α,v0=$v₀.jld2"
    jldsave(filename, true; energies)
end
