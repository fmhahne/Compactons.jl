using Printf
using JLD2, CodecZlib
using UnPack
using Compactons

mkpath("data/kink_oscillon_scattering_energies")

let l = 1.0
    filename = @sprintf "data/kink_oscillon_scattering_energies/l=%.2f.jld2" l
    energies = Dict()

    for V ∈ 0.00:0.01:0.99
        for α ∈ 0.00:0.01:0.99
            print("Simulating l=$l, V=$V, α=$α … ")
            energies[(V=V, α=α)] = kink_oscillon_scattering_energies(l, V, α)
            println("done")
        end
    end

    jldsave(filename, true; energies)
end
