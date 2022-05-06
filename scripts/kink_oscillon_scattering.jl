using Printf
using JLD2, CodecZlib
using UnPack
using Compactons

mkpath("data/kink_oscillon_scattering")

ls = 0.5:0.5:2.0
Vs = 0.0:0.25:0.75
αs = 0.0:0.25:0.75
v₀s = 0.0:0.25:0.75

let V = 0.75, α = 0.0, v₀ = 0.0
    for l ∈ ls
        filename = @sprintf "data/kink_oscillon_scattering/l=%.2f,V=%.2f,alpha=%.2f,v0=%.2f.jld2" l V α v₀
        print("Producing $filename … ")

        @unpack x, t, η, H, E₁, E₂, E₃ = kink_oscillon_scattering(l, V, α, v₀; dx=1e-3, sampling=10)
        jldsave(filename, true; x, t, η, H, E₁, E₂, E₃)
        println("done")
    end
end

let l = 0.5, α = 0.0, v₀ = 0.0
    for V ∈ Vs
        filename = @sprintf "data/kink_oscillon_scattering/l=%.2f,V=%.2f,alpha=%.2f,v0=%.2f.jld2" l V α v₀
        print("Producing $filename … ")

        @unpack x, t, η, H, E₁, E₂, E₃ = kink_oscillon_scattering(l, V, α, v₀; dx=1e-3, sampling=10)
        jldsave(filename, true; x, t, η, H, E₁, E₂, E₃)
        println("done")
    end
end

let l = 1, V = 0.6, v₀ = 0.0
    for α ∈ αs
        filename = @sprintf "data/kink_oscillon_scattering/l=%.2f,V=%.2f,alpha=%.2f,v0=%.2f.jld2" l V α v₀
        print("Producing $filename … ")

        @unpack x, t, η, H, E₁, E₂, E₃ = kink_oscillon_scattering(l, V, α, v₀; dx=1e-3, sampling=10)
        jldsave(filename, true; x, t, η, H, E₁, E₂, E₃)
        println("done")
    end
end

let l = 0.75, V = 0.8, α = 0
    for v₀ ∈ v₀s
        filename = @sprintf "data/kink_oscillon_scattering/l=%.2f,V=%.2f,alpha=%.2f,v0=%.2f.jld2" l V α v₀
        print("Producing $filename … ")

        @unpack x, t, η, H, E₁, E₂, E₃ = kink_oscillon_scattering(l, V, α, v₀; dx=1e-3, sampling=10)
        jldsave(filename, true; x, t, η, H, E₁, E₂, E₃)
        println("done")
    end
end
