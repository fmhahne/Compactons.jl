using Printf
using JLD2, CodecZlib
using Compactons

mkpath("data/kink_antikink_scattering")

for V ∈ 0.00:1e-2:0.99
    filename = @sprintf "data/kink_antikink_scattering/V=%.2f.jld2" V
    print("Producing $filename … ")

    x, t, η, H = kink_antikink_scattering(V)
    jldsave(filename, true; x, t, η, H)
    println("done")
end
