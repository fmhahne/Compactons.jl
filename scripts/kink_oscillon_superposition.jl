using Printf
using JLD2, CodecZlib
using Compactons

mkpath("data/kink_oscillon_superposition")

for l ∈ 0.5:0.5:3.0
    for α ∈ 0.00:0.25:0.75
        filename = @sprintf "data/kink_oscillon_superposition/l=%.2f,alpha=%.2f.jld2" l α
        print("Producing $filename … ")

        x, t, η, H = kink_oscillon_superposition(l, α)
        jldsave(filename, true; x, t, η, H)
        println("done")
    end
end
