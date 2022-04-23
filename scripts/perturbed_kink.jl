using Printf
using JLD2, CodecZlib
using Compactons

mkpath("data/perturbed_kink")

for ϵ ∈ -0.50:0.05:0.50
    filename = @sprintf "data/perturbed_kink/eps=%.2f.jld2" ϵ
    print("Producing $filename … ")

    x, t, η, H = perturbed_kink(ϵ)
    jldsave(filename, true; x, t, η, H)
    println("done")
end
