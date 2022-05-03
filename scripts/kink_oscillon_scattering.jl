using Printf
using JLD2, CodecZlib
using UnPack
using Compactons

mkpath("data/kink_oscillon_scattering")

for l ∈ 0.5:0.5:3.0
    for V ∈ 0.0:1e-1:0.9
        for α ∈ 0.00:0.25:0.75
            filename = @sprintf "data/kink_oscillon_scattering/l=%.2f,V=%.2f,alpha=%.2f.jld2" l V α
            print("Producing $filename … ")

            @unpack x, t, η, H, E₁, E₂, E₃ = kink_oscillon_scattering(l, V, α)
            jldsave(filename, true; x, t, η, H, E₁, E₂, E₃)
            println("done")
        end
    end
end
